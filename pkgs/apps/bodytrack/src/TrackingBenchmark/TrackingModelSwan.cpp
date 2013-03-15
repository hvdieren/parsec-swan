//-------------------------------------------------------------
//      ____                        _      _
//     / ___|____ _   _ ____   ____| |__  | |
//    | |   / ___| | | |  _  \/ ___|  _  \| |
//    | |___| |  | |_| | | | | |___| | | ||_|
//     \____|_|  \_____|_| |_|\____|_| |_|(_) Media benchmarks
//                  
//	  2006, Intel Corporation, licensed under Apache 2.0 
//
//  file : TrackingModelSwan.cpp
//  author :      Scott Ettinger - scott.m.ettinger@intel.com
//  description : Observation model for kinematic tree body 
//				  tracking threaded with Pipe	
//				  
//  modified : 
//--------------------------------------------------------------

#if defined(HAVE_CONFIG_H)
# include "config.h"
#endif

#include "TrackingModelSwan.h"
#include <vector>
#include <string>
#include <stdio.h>
#include "system.h"
#include "swan/wf_interface.h"
#include <iostream>
#include "swan/tbb_interface.h"
#include "CovarianceMatrix.h"
using namespace obj;
using namespace std;


#if (STACK_FRAME_SIZE!=32768)
	error STACK_FRAME_SIZE!!!
#endif

#define GRAIN_SIZE 8

//-----------------TIMERS---------------------
void stimer_tick(stimer_t *timer)
{
        gettimeofday(&(timer->start), 0);
}

float stimer_tuck(stimer_t *timer, const char *msg)
{
        gettimeofday(&(timer->end), 0);

        timer->diff = (timer->end.tv_sec - timer->start.tv_sec)
                                + (timer->end.tv_usec - timer->start.tv_usec) * 0.000001;

        if (msg)
                printf("%s: %.3f seconds\n", msg, timer->diff);

        return timer->diff;
}
//-------------END-TIMERS---------------------




//templated conversion to string with field width
template<class T>
inline string str(T n, int width = 0, char pad = '0')
{	stringstream ss;
	ss << setw(width) << setfill(pad) << n;
	return ss.str();
}


//Initialize the tracking model parameters
bool TrackingModelPipe::Initialize(const string &path, int cameras, int layers)
{
	mPath = path;
	mNCameras = cameras;
	mFGMaps.resize(cameras);
	mEdgeMaps.resize(cameras);
	vector<string> calibFiles(cameras);											//set camera calibration file paths
	for(int i = 0; i < cameras; i++)
		calibFiles[i] = path + "CALIB" + DIR_SEPARATOR + "Camera" + str(i + 1) + ".cal";
	if(!InitCameras(calibFiles)) return false;									//initialize camera calibration parameters				
	if(!InitGeometry(path + "BodyShapeParameters.txt")) return false;			//initialize geometry parameters
	if(!LoadInitialState(path + "InitialPose.txt")) return false;				//initialize body pose angles and translations
	if(!LoadPoseParameters(path + "PoseParameters.txt")) return false;			//initialize pose statistics
	GenerateStDevMatrices(layers, mPoses[0].Params(), mStdDevs);				//generate annealing rates for particle filter using pose parameters
	return true;
}

void initialize_wrapper(TrackingModelPipe *model, string *path, int cameras, int layers, int *res)
{
	*res = 0;
	if(model->Initialize(*path, cameras, layers))	*res = 1;
	return;
}

//Allocate space for multi-threaded code
void TrackingModelPipe::SetNumThreads(int n)
{	
	mPoses.resize(n);  
	mBodies.resize(n); 
	mProjections.resize(n); 
	mImageMeasurements.resize(n);
}

void SetNumThreads_wrapper(TrackingModelPipe *model, int n)
{
	model->SetNumThreads(n);
}


void SetNumFrames_wrapper(TrackingModelPipe* model, int frames)
{
	model->SetNumFrames(frames);
}

//object used with TBB parallel for to do the gradient magnitude and thresholding operation
template<class T>
class DoGradientMagThresholdPipe{

  const FlexImage<T,1> *p_tmp, *p_src; 
  float threshold; 
  
public:
	void operator()(const blocked_range<int>& r) const
	{
		leaf_call(DoGradientMagThresholdPipe<T>::leaf_func, this, &r);
		//DoGradientMagThresholdPipe<T>::leaf_func( this, &r);
	}

	static void leaf_func(const DoGradientMagThresholdPipe<T> * self, const blocked_range<int> * r)
	{
		const FlexImage<T,1> &tmp = *self->p_tmp;
		const FlexImage<T,1> &src = *self->p_src; 
		//int y = r->begin();
		for(int y = r->begin(); y < r->end(); y++)//y+=r->grainsize())//y++)
		{
			Im8u *p = &src(1,y), *ph = &src(1,y - 1), *pl = &src(1,y + 1), *pr = &tmp(1,y);
			for(int x = 1; x < src.Width() - 1; x++)
			{	
				float xg = -0.125f * ph[-1] + 0.125f * ph[1] - 0.250f * p[-1] + 0.250f * p[1] - 0.125f * pl[-1] + 0.125f * pl[1];	//calc x and y gradients
				float yg = -0.125f * ph[-1] - 0.250f * ph[0] - 0.125f * ph[1] + 0.125f * pl[-1] + 0.250f * pl[0] + 0.125f * pl[1];
				float mag = xg * xg + yg * yg;				//calc magnitude and threshold
				*pr = (mag < self->threshold) ? 0 : 255;
				 p++; ph++; pl++; pr++;
			}
		}
	}
	DoGradientMagThresholdPipe(FlexImage<T,1> *_source, FlexImage<T,1> *_tmp, float _threshold) : 
		p_tmp(_tmp), p_src(_source), threshold(_threshold) {}

};




//object used with TBB parallel for to do row filtering
template<class T, class T2> 
  class DoFlexFilterRowVPipe{
    
    FlexImage<T, 1> *p_src, *p_dst; 
    T2 *kernel;
    int pn;   
 
  public:
    
    void operator() ( const blocked_range<int>& r) const
	{
 		leaf_call(DoFlexFilterRowVPipe<T,T2>::leaf_func, this, &r);
		//DoFlexFilterRowVPipe<T,T2>::leaf_func( this, &r);
	}
	static void leaf_func(const DoFlexFilterRowVPipe<T, T2> * self, const blocked_range<int> *r)
	{
      int n = self->pn;
      FlexImage<T, 1> &src = *self->p_src; 
      FlexImage<T,1> &dst = *self->p_dst; 
      int source_width = src.Width();

      for(int y = r->begin(); y < r->end(); y++)//y+=r->grainsize())// y++)
	  {
		T *psrc = &src(n, y);
		T *pdst = &dst(n, y); 
		for(int x = n; x < source_width - n; x++)
	    {
	      int k = 0; 
	      T2 acc = 0; 
	      for(int i = -n; i <= n; i++)
			acc += (T2)(psrc[i] * self->kernel[k++]);
	      
			*pdst = (T)acc;
			pdst++;
			psrc++;
	    }
	  }
    }

    DoFlexFilterRowVPipe(FlexImage<T,1> *_source, FlexImage<T,1> *_dest, T2 *_kernel, int _n) :
	p_src(_source), p_dst(_dest), kernel(_kernel), pn(_n) {}
};


//object used with TBB parallel for to do column filtering
template<class T, class T2> 
  class DoFlexFilterColumnVPipe{
    
    FlexImage<T, 1> *p_src, *p_dst; 
    T2 *kernel;
    int pn;

  public:
    
    void operator() ( const blocked_range<int>& r) const
	{
		leaf_call(DoFlexFilterColumnVPipe<T,T2>::leaf_func, this, &r);
		//DoFlexFilterColumnVPipe<T,T2>::leaf_func( this, &r);
    }

	static void leaf_func(const DoFlexFilterColumnVPipe<T, T2> * self, const blocked_range<int> *r)
	{
		FlexImage<T,1> &src = *self->p_src; 
		FlexImage<T,1> &dst = *self->p_dst;
		int source_width = src.Width(); 
		int sb = src.StepBytes();
		int n = self->pn;

		for(int y = r->begin(); y < r->end(); y++)//y+=r->grainsize())//y++)
		{
		  T *psrc = &src(0, y);
		  T *pdst = &dst(0, y); 
		  for(int x = 0; x < source_width; x++)
		  {
			int k = 0; 
			T2 acc = 0; 
			for(int i = -n; i <= n; i++)
				acc += (T2)(*(T *)((char *)psrc + sb * i) * self->kernel[k++]);
		      
			*pdst = (T)acc;
			pdst++;
			psrc++;
		   }
		}
	}
    DoFlexFilterColumnVPipe(FlexImage<T,1> *_source, FlexImage<T,1> *_dest, T2 *_kernel, int _n) : 
		p_src(_source), p_dst(_dest), kernel(_kernel), pn(_n) {}
																				
};

//TBB threaded - 1D filter Row wise 1 channel any type data or kernel valid pixels only
template<class T, class T2>
bool FlexFilterRowVPipe(FlexImage<T,1> &src, FlexImage<T,1> &dst, T2 *kernel, int kernelSize, bool allocate = true)
{
	//cout<<"sth sunarthsh FlexFilterRowVPipe!!!"<<endl;
	if(kernelSize % 2 == 0)									//enforce odd length kernels
		return false;
	if(allocate)
		dst.Reallocate(src.Size());
	dst.Set((T)0);
	int n = kernelSize / 2, h = src.Height();
	
	FlexImage<T,1> *psrc = &src;
	FlexImage<T,1> *pdst = &dst;
	
	blocked_range<int> range(0, h, GRAIN_SIZE);
	//Pipe parallel for
	parallel_for(range, DoFlexFilterRowVPipe<T,T2>(psrc, pdst, kernel, n));
	//DoFlexFilterRowVPipe<T,T2>(psrc, pdst, kernel, n)(range);
	
	return true;
}

//TBB threaded - 1D filter Column wise 1 channel any type data or kernel valid pixels only
template<class T, class T2>
inline bool FlexFilterColumnVPipe(FlexImage<T,1> &src, FlexImage<T,1> &dst, T2 *kernel, int kernelSize, bool allocate = true)
{
	//cout<<"sthn sunarthsh FlexFilterColumnVPipe!!!"<<endl;
	if(kernelSize % 2 == 0)									//enforce odd length kernels
		return false;
	if(allocate)
		dst.Reallocate(src.Size());
	dst.Set((T)0);
	int n = kernelSize / 2;
	int h = src.Height() - n;

	FlexImage<T,1> *psrc = &src;
	FlexImage<T,1> *pdst = &dst;
	//Pipe parallel for
	blocked_range<int> range(n, h, GRAIN_SIZE); 
	parallel_for(range, DoFlexFilterColumnVPipe<T,T2>(psrc, pdst, kernel, n));
	//DoFlexFilterColumnVPipe<T,T2>(psrc, pdst, kernel, n)(range);

	return true;
}

void GaussianBlurPipe(FlexImage8u *src, FlexImage8u *dst)
{

	//std::cout<<"sth sunarthsh GaussianBlurPipe!!!"<<endl;
    FlexImage8u tmp(src->Size());
  	float k[] = {0.12149085090552f, 0.14203719483447f, 0.15599734045770f, 0.16094922760463f, 0.15599734045770f, 0.14203719483447f, 						0.12149085090552f};
	//FlexFilterRowVPipe(*src, tmp, k, 7);

	bool res;
	//inlining FlexFilterRowVPipe
	//-------------------------------------------------------
//bool FlexFilterRowVPipe(FlexImage<T,1> &src, FlexImage<T,1> &dst, T2 *kernel, int kernelSize, bool allocate = true)
	bool allocate = true;
	if(7 % 2 == 0)									//enforce odd length kernels
		res = false;
	if(allocate)
		tmp.Reallocate(src->Size());
	tmp.Set((Im8u)0);
	int n = 7 / 2, h = src->Height();
	
	FlexImage8u *psrc = src;
	FlexImage8u *pdst = &tmp;
	
	blocked_range<int> range(0, h, GRAIN_SIZE);
	//Pipe parallel for
	parallel_for(range, DoFlexFilterRowVPipe<Im8u,float>(psrc, pdst, k, n));
	//DoFlexFilterRowVPipe<Im8u,float>(psrc, pdst, k, n)(range);
	
	res = true;
	//-------------------------------------------------------


	FlexFilterColumnVPipe(tmp, *dst, k, 7);
}

#if 0
FlexImage8u* GradientMagThresholdPipe(FlexImage8u *src, float* threshold)
{
	//std::cout<<"sth sunarthsh GradientMagThresholdPipe!!!"<<endl;
	FlexImage8u r(src->Size());
	ZeroBorder(r);
	FlexImage8u * s = src;
	FlexImage8u * rr = &r;
	blocked_range<int> range(1, src->Height() - 1);

	//Pipe parallel for  
	//parallel_for(range, DoGradientMagThresholdPipe<Im8u>(s, rr, *threshold));
		DoGradientMagThresholdPipe<Im8u>(s, rr, *threshold)(range);
	return &r; 
}
#endif


//Generate an edge map from the original camera image
void TrackingModelPipe::CreateEdgeMap(FlexImage8u &src, FlexImage8u &dst)
{

	//std::cout<<"sth sunarthsh CreateEdgeMap!!!"<<std::endl;
	float f = 16.0f;
  //FlexImage8u *gr = GradientMagThresholdPipe(&src, &f);		//calc gradient magnitude and threshold
  //GaussianBlurPipe(gr, &dst);								//Blur to create distance error map

	//inlining GradientMagThresholdPipe
	//---------------------------------------------------------
	FlexImage8u r(src.Size());
	ZeroBorder(r);
	FlexImage8u * s = &src;
	FlexImage8u * rr = &r;
	blocked_range<int> range(1, src.Height() - 1);

	//Pipe parallel for  
	parallel_for(range, DoGradientMagThresholdPipe<Im8u>(s, rr, f));
		//DoGradientMagThresholdPipe<Im8u>(s, rr, f)(range);
	//----------------------------------------------------------
	//GaussianBlurPipe(rr, &dst);								//Blur to create distance error map

	//inlining GaussianBlurPipe
	//----------------------------------------------------------
	 FlexImage8u tmp(rr->Size());
  	float k[] = {0.12149085090552f, 0.14203719483447f, 0.15599734045770f, 0.16094922760463f, 0.15599734045770f, 0.14203719483447f, 						0.12149085090552f};
	//FlexFilterRowVPipe(*rr, tmp, k, 7);
	//inlining FlexFilterRowVPipe
	//-------------------------------------------------------
	//bool FlexFilterRowVPipe(FlexImage<T,1> &src, FlexImage<T,1> &dst, T2 *kernel, int kernelSize, bool allocate = true)
	bool res;
	bool allocate = true;
	//if(7 % 2 == 0)									//enforce odd length kernels
	//	res = false;
	//if(allocate)
		tmp.Reallocate(rr->Size());
	tmp.Set((Im8u)0);
	int n = 3, h = rr->Height(); //7 / 2, h = rr->Height();
	
	FlexImage8u *psrc = rr;
	FlexImage8u *pdst = &tmp;
	
	blocked_range<int> range2(0, h, GRAIN_SIZE);
	//Pipe parallel for
	parallel_for(range, DoFlexFilterRowVPipe<Im8u,float>(psrc, pdst, k, n));
	//DoFlexFilterRowVPipe<Im8u,float>(psrc, pdst, k, n)(range2);
	
	res = true;
	//-------------------------------------------------------
	

	FlexFilterColumnVPipe(tmp, dst, k, 7);
	//----------------------------------------------------------
}

//Generate an edge map from the original camera image
inline void ComputeEdgeMapsPipe(FlexImage8u *src, FlexImage8u *dst)
{
	//cout<<"sthn ComputeEdgeMapsPipe!!!"<<endl;
	FlexImage8u* psrc = src;
	float *f = (float*)malloc(sizeof(float));
	*f = 16.0f;
  //FlexImage8u *gr = GradientMagThresholdPipe( psrc, f);		//calc gradient magnitude and threshold

	//inlining GradientMagThresholdPipe
	//---------------------------------------------------------
	FlexImage8u r(src->Size());
	ZeroBorder(r);
	FlexImage8u * s = src;
	FlexImage8u * rr = &r;
	blocked_range<int> range(1, src->Height() - 1);

	//Pipe parallel for  
	parallel_for(range, DoGradientMagThresholdPipe<Im8u>(s, rr, *f));
		//DoGradientMagThresholdPipe<Im8u>(s, rr, *f)(range);
	//----------------------------------------------------------
  //GaussianBlurPipe(rr, dst);								//Blur to create distance error map

	//inlining GaussianBlurPipe
	//----------------------------------------------------------
	 FlexImage8u tmp(rr->Size());
  	float k[] = {0.12149085090552f, 0.14203719483447f, 0.15599734045770f, 0.16094922760463f, 0.15599734045770f, 0.14203719483447f, 						0.12149085090552f};
	//FlexFilterRowVPipe(*rr, tmp, k, 7);
	//inlining FlexFilterRowVPipe
	//-------------------------------------------------------
	//bool FlexFilterRowVPipe(FlexImage<T,1> &src, FlexImage<T,1> &dst, T2 *kernel, int kernelSize, bool allocate = true)
	bool res;
	bool allocate = true;
	//if(7 % 2 == 0)									//enforce odd length kernels
	//	res = false;
	//if(allocate)
		tmp.Reallocate(rr->Size());
	tmp.Set((Im8u)0);
	int n = 3, h = rr->Height();// 7 / 2, h = rr->Height();
	
	//FlexImage8u *psrc = rr;
	FlexImage8u *pdst = &tmp;
	
	blocked_range<int> range2(0, h, GRAIN_SIZE);
	//Pipe parallel for
	parallel_for(range2, DoFlexFilterRowVPipe<Im8u,float>(rr, pdst, k, n));
	//DoFlexFilterRowVPipe<Im8u,float>(rr, pdst, k, n)(range2);
	
	res = true;
	//-------------------------------------------------------
	
	FlexFilterColumnVPipe(tmp, *dst, k, 7);
	//----------------------------------------------------------
}


//TBB block class to process images in parallel
class DoProcessImages	{

	ImageSet *mEdgeMaps;
	vector<string> *mFGfiles, *mImageFiles;
	BinaryImageSet *mFGmaps;

public:

	DoProcessImages(vector<string> *FGfiles, vector<string> *ImageFiles, ImageSet *edgeMaps, BinaryImageSet *FGMaps)
		: mEdgeMaps(edgeMaps), mFGfiles(FGfiles), mImageFiles(ImageFiles), mFGmaps(FGMaps) {};

	void operator() (const blocked_range<int> &r) const
	{
		//call(DoProcessImages::leaf_func, this, &r);
		DoProcessImages::leaf_func( this, &r);				//to avoid call we made leaf_func inline in this case
	}
	
	static inline void leaf_func(const DoProcessImages * self, const blocked_range<int> *r)
	{
		FlexImage8u im;
		for(int i = r->begin(); i < r->end(); i++)//i+=r->grainsize())//i++)
		{	if(!FlexLoadBMP((*self->mFGfiles)[i].c_str(), im))							//Load foreground maps and raw images
			{	cout << "Unable to load image: " << (*self->mFGfiles)[i].c_str() << endl;
				return;
			}	
			(*self->mFGmaps)[i].ConvertToBinary(im);									//binarize foreground maps to 0 and 1
			if(!FlexLoadBMP((*self->mImageFiles)[i].c_str(), im))
			{	cout << "Unable to load image: " << (*self->mImageFiles)[i].c_str() << endl;
				return;
			}
			//call(ComputeEdgeMapsPipe, &im, &(*self->mEdgeMaps)[i]);									//Create edge maps
			ComputeEdgeMapsPipe( &im, &(*self->mEdgeMaps)[i]);
		}
	}
};


//load and process all images for new observation at a given time(frame)
bool TrackingModelPipe::GetObservation(float timeval)
{
	int frame = (int)timeval;													//generate image filenames
	int n = mCameras.GetCameraCount();
	vector<string> FGfiles(n), ImageFiles(n);
	for(int i = 0; i < n; i++)													
	{	FGfiles[i] = mPath + "FG" + str(i + 1) + DIR_SEPARATOR + "image" + str(frame, 4) + ".bmp";
		ImageFiles[i] = mPath + "CAM" + str(i + 1) + DIR_SEPARATOR + "image" + str(frame, 4) + ".bmp";
	}
	FlexImage8u im;
	for(int i = 0; i < (int)FGfiles.size(); i++)
	{	if(!FlexLoadBMP(FGfiles[i].c_str(), im))								//Load foreground maps and raw images
		{	cout << "Unable to load image: " << FGfiles[i].c_str() << endl;
			return 0;
		}
		mFGMaps[i].ConvertToBinary(im);											//binarize foreground maps to 0 and 1
		if(!FlexLoadBMP(ImageFiles[i].c_str(), im))
		{	cout << "Unable to load image: " << ImageFiles[i].c_str() << endl;
			return 0;
		}
		//call(ComputeEdgeMapsPipe, &im, &(mEdgeMaps[i]));									//Create edge maps
		ComputeEdgeMapsPipe( &im, &(mEdgeMaps[i]));
	}
	return 1;
}

void GetObservation_wrapper(TrackingModelPipe* model, float timeval, int *res)
{
	//int res;
	*res = model->GetObservation( timeval );
}


void TrackingModelPipe::SetObservation(ImageSetToken* token) 
{
	mEdgeMaps = 
	token->edgeMaps; 
	mFGMaps = token->FGmaps;
}

void setObservation_wrapper(TrackingModelPipe *model, ImageSetToken *token)
{
	model->SetObservation(token);
}

/*
void SetMFGMaps_wrapper(TrackingModelPipe* model, std::vector<BinaryImage > FGMaps)
{
	model->SetMFGMaps(FGMaps);
}
*/

//TBB pipeline stage function STAGE 1
//void *TrackingModelPipe::operator ()(void *inToken)
//spawn(stage1, (indep<TrackingModelPipe>)modelObj, (outdep<ImageSetToken*>)token);
//void stage1(TrackingModelPipe * modelObj, ImageSetToken* tokenObj);

/*
int stage1_wrapper(TrackingModelPipe * modelObj, outdep<ImageSetToken*> tokenObj, int currFrame)
{
	//TrackingModelPipe* model = (TrackingModelPipe*) modelObj;
	//tokenObj = (ImageSetToken*) new ImageSetToken;//(ImageSetToken*)calloc(1, sizeof(ImageSetToken));
	ImageSetToken* tk = (ImageSetToken*)tokenObj;
	//leaf_call(stage1, (TrackingModelPipe*) modelObj, (ImageSetToken*)tk);
//	stage1( (TrackingModelPipe*) modelObj, (ImageSetToken*)tk, currFrame);

}
*/

void stage1(/*indep<TrackingModelPipe *>*/TrackingModelPipe * modelObj, outdep<ImageSetToken*> tokenObj, int currFrame, inoutdep<void>)
{
	//stimer_t tmr_stage1;
	//stimer_tick(&tmr_stage1);
	TrackingModelPipe* model = (TrackingModelPipe*) modelObj;
	int n = model->mCameras.GetCameraCount();
	// ImageSetToken *token = (ImageSetToken *)tokenObj;
	ImageSetToken *token = new ImageSetToken;

	token->edgeMaps.resize(n);
	token->FGmaps.resize(n);
	
	cout << "Processing frame : " << currFrame <<endl;

	vector<string> FGfiles(n), ImageFiles(n);
	for(int i = 0; i < n; i++)
	{	FGfiles[i] = model->mPath + "FG" + str(i + 1) + DIR_SEPARATOR + "image" + str(currFrame, 4) + ".bmp";
		ImageFiles[i] = model->mPath + "CAM" + str(i + 1) + DIR_SEPARATOR + "image" + str(currFrame, 4) + ".bmp";
	}

	blocked_range<int> range(0, n);
	//Swan parallel for
	//parallel_for(range, DoProcessImages( pFGfiles, pImageFiles, pedgeMaps, pFGmaps));
	parallel_for(range, DoProcessImages( &FGfiles, &ImageFiles, &(token->edgeMaps), &(token->FGmaps)));
	//DoProcessImages( pFGfiles, pImageFiles, pedgeMaps, pFGmaps)(range);

	tokenObj = token;
	//(TrackingModel*)model->SetMFGMaps(token->FGmaps);
	//stimer_tuck(&tmr_stage1, "Stage 1 took");
	return;
}
