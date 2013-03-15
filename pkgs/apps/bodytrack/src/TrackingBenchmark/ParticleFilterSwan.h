 //------------------------------------------------------------------------
//      ____                        _      _
//     / ___|____ _   _ ____   ____| |__  | |
//    | |   / ___| | | |  _  \/ ___|  _  \| |
//    | |___| |  | |_| | | | | |___| | | ||_|
//     \____|_|  \_____|_| |_|\____|_| |_|(_) Media benchmarks
//                           
//	 © 2006, Intel Corporation, licensed under Apache 2.0 
//
//  file :	 ParticleFilterSwan.h
//  author : 
//
//  description : Intel TBB parallelized version of the particle filter
//                    object dervied from ParticleFilter.h
//  modified : 
//--------------------------------------------------------------------------

#ifndef PARTICLEFILTERPIPE_H
#define PARTICLEFILTERPIPE_H

#if defined(HAVE_CONFIG_H)
# include "config.h"
#endif

#include "ParticleFilter.h"
#include "TrackingModelSwan.h"
//#include "TBBtypes.h"
//#include "tbb/task_scheduler_init.h"
//#include "tbb/parallel_for.h"
//#include "tbb/parallel_reduce.h"
//#include "tbb/blocked_range.h"
//#include "tbb/pipeline.h"
#include "/home/kchronaki/scheduler/scheduler/tbb_interface.h"
#include <fstream>
#include <iostream>

#if (STACK_FRAME_SIZE!=32768)
	error STACK_FRAME_SIZE!!!
#endif

#undef min

#define WORKUNIT_SIZE_NEWPARTICLES 32
#define WORKUNIT_SIZE_CALCWEIGHTS 32
using namespace std;
using namespace obj;
template<class T> 
class ParticleFilterPipe : public ParticleFilter<T> {

  	using ParticleFilter<T>:: mModel;
	using ParticleFilter<T>:: mWeights;
	using ParticleFilter<T>:: mParticles;
	using ParticleFilter<T>:: mNewParticles;
	using ParticleFilter<T>:: mBestParticle;
	using ParticleFilter<T>:: mNParticles;
	using ParticleFilter<T>:: mMinParticles;
	using ParticleFilter<T>:: mBins;  
	using ParticleFilter<T>:: mRnd;
        using ParticleFilter<T>:: mInitialized;
        using ParticleFilter<T>:: mCdf;
        using ParticleFilter<T>:: mSamples;
	typedef typename ParticleFilter<T>::fpType fpType;
	typedef typename ParticleFilter<T>::Vectorf Vectorf;

private:	
    T* getModel() { return mModel; };

protected:
	std::vector<int> mIndex;																//list of particles to regenerate
	//std::ofstream mPoseOutFile;																//output pose file
	bool mOutputBMP;																		//write bitmap output flag
	unsigned int mFrame;																	//current frame being processed

	std::vector<unsigned char> mValid;														//storage for valid particle flags

	//calculate particle weights - threaded version 
	

	


	public:
	//TBB pipeline stage function
	static void stage_operator(ParticleFilterPipe *, ImageSetToken *token);

public:
	void CalcWeights(std::vector<Vectorf > &particles)__attribute__((always_inline));	//calculate particle weights based on model likelihood	
	//Set number of particles to n and generate initial values
	void InitializeParticles(int n);		
	//New particle generation - threaded version 
	void GenerateNewParticles(int k)__attribute__((always_inline));	

	std::ofstream mPoseOutFile;																//output pose file
	ParticleFilterPipe() : mOutputBMP(false), mFrame(0) {};
	virtual ~ParticleFilterPipe() { };

	void WritePose(std::ostream &f, std::vector<float> &pose);
	 T* getmModel() { return mModel; };
	
	//sets
	void setOutputFile(const char *fname) {mPoseOutFile.open(fname); };
	void setOutputBMP(bool flag) {mOutputBMP = flag; };
	
	std::ofstream getmPoseOutFile()	{ return mPoseOutFile; };
	bool getmOutputBMP()	{ return mOutputBMP; };
	unsigned int getmFrame()	{ return mFrame; }
	unsigned int increasemFrame() { return ++mFrame; }
	T GetModel() { return *mModel; };
	//Particle filter update
	bool Update(fpType timeval)__attribute__((always_inline));


protected :

	//----------------------------------- TBB block computing objects ------------------------------------

	//particle generation block computing object
	template<class S>
	class DoGenerateNewParticlesPipe {

	  S *mModel;
	  std::vector<Vectorf> &mNewParticles, &mParticles;
	  std::vector<RandomGenerator> &mRnd; 
	  std::vector<int> &mIndex;
	  
	  int k;
		//distribute particle randomly according to given standard deviations
		inline void AddGaussianNoise(Vectorf &p, const Vectorf &stdDevs, RandomGenerator &rnd) const
		{	for(uint i = 0; i < stdDevs.size(); i++)
				p[i] += (fpType)rnd.RandN() * stdDevs[i];				
		}
	 public: 
		DoGenerateNewParticlesPipe(int _k, S* _model, std::vector<Vectorf> &_mNewParticles,std::vector<Vectorf> &_mParticles, std::vector<RandomGenerator> &_mRnd, std::vector<int> &_mIndex) 
				: mModel(_model), mNewParticles(_mNewParticles), mParticles(_mParticles), mRnd(_mRnd), mIndex(_mIndex), k(_k) {} 
		
		void operator() (const blocked_range<int> &r ) const 
		 {	
			leaf_call(DoGenerateNewParticlesPipe<S>::leaf_func, this, &r);
			//DoGenerateNewParticlesPipe<S>::leaf_func( this, &r);
			
		 }
		 
		 static void leaf_func(const DoGenerateNewParticlesPipe<S> * self, const blocked_range<int> * r)
		 {
			for( int i = r->begin(); i < r->end(); i++ )
			{	
				self->mNewParticles[i] = self->mParticles[self->mIndex[i]];
				self->AddGaussianNoise(self->mNewParticles[i], self->mModel->StdDevs()[self->k], self->mRnd[i]);
			}
		 }
	};

	//likelihood block computing object
	template<class S>
	class DoCalcLikelihoods {
	private:
			S *mModel;
			Vectorf *mParticles;
			fpType *mWeights;
			unsigned char *mValid;
	public:

		DoCalcLikelihoods(S *model, Vectorf *particles, fpType *weights, unsigned char *valid) : mModel(model), mParticles(particles), mWeights(weights), mValid(valid) {};

		void operator()( const blocked_range<int>& r ) const
		{
			leaf_call( DoCalcLikelihoods<S>::leaf_func, this, &r );
			//DoCalcLikelihoods<S>::leaf_func( this, &r );
		}		
		static void leaf_func(const DoCalcLikelihoods<S> * self, const blocked_range<int> * r )		
		{
			for( int i = r->begin(); i<r->end(); ++i)//i+=r->grainsize())//++i )
			{	bool vflag;
				//std::cout<<"in leaf_func DoCalcLikelihoods (calculates mWeights also...) mWeights["<<i<<"] = "<<self->mWeights[i]<<endl;
				self->mWeights[i] = self->mModel->LogLikelihood(self->mParticles[i], vflag, i);						//compute log-likelihood weights for each particle
				//std::cout<<"in leaf_func DoCalcLikelihoods AFTER CALCULATION: mWeights["<<i<<"] = "<<self->mWeights[i]<<endl;
				self->mValid[i] = vflag ? 1 : 0;
			}
		}
	};

};

template<class T> 
void ParticleFilterPipe<T>::CalcWeights(std::vector<Vectorf> &particles){


	//printf("edw in particlefilterpipe, calcweights!\n");
	//std::cout<<"in CALCWEIGHTS -> "<<std::endl; 
	mBestParticle =0; 
	fpType total =0, best =0, minWeight = 1e30f, annealingFactor = 1;
	mWeights.resize(particles.size());

	//p_particles = &particles; 
	//p_weights = &mWeights;//Returns pointer to mWeights
	int np = (int) particles.size(); 
	std::vector<unsigned char> &valid = mValid;
	valid.resize(np); //Resize global valid to num. particles
	uint i = 0;

	unsigned char * pvalid = &valid[0];
	float * pmWeights = &mWeights[0];
	Vectorf * pparticles = &particles[0];

	//TrackingModelPipe * pmModel = &mModel;
	int startpoint = 0;
	int step = WORKUNIT_SIZE_CALCWEIGHTS;
	//parallel code to calculate likelihoods										//&particles[0], &mWeights[0], &valid[0]
	//parallel_for(startpoint, np, step, DoCalcLikelihoods<T>(mModel, pparticles, pmWeights, pvalid));

	blocked_range<int> range(startpoint, np, step);
	parallel_for(range, DoCalcLikelihoods<T>(mModel, pparticles, pmWeights, pvalid));
	//for(int j = startpoint; j<np; j+=step)
		//DoCalcLikelihoods<T>(mModel, pparticles, pmWeights, pvalid)(range);
	i = 0;
	while(i < particles.size())
	{	
		if(!valid[i])//if not valid(model prior), remove the particle from the list
		{	
		  particles[i] = particles[particles.size() - 1];
		  mWeights[i] = mWeights[particles.size() - 1];
		  valid[i] = valid[valid.size() - 1];
		  particles.pop_back(); mWeights.pop_back(); valid.pop_back();
		}
		else{
		  minWeight = std::min(mWeights[i++], minWeight);				//find minimum weight
		}
	//std::cout<<"in CALCWEIGHTS -> mWeights["<<i<<"] = "<<mWeights[i]; 
	}
	  //std::cout<<"in CALCWEIGHTS -> MINWEIGHT = "<<minWeight<<std::endl; 
	//bail out if not enough valid particles
	if((int)particles.size() < mMinParticles) return;					
	mWeights -= minWeight; //shift weights to zero for numerical stability
	if(mModel->StdDevs().size() > 1) 
		annealingFactor = BetaAnnealingFactor(mWeights, 0.5f);			//calculate annealing factor if more than 1 step
	//std::cout<<"in CALCWEIGHTS -> mWeights.size = "<<mWeights.size()<<std::endl; 
	for(i = 0; i < mWeights.size(); i++)
	{	
	  double wa = annealingFactor * mWeights[i];
	  //std::cout<<"in CALCWEIGHTS -> mWeights before calculation = "<<mWeights[i]<<" annealingFactor = "<<annealingFactor<<std::endl; 
	  mWeights[i] = (float)exp(wa);									//exponentiate log-likelihoods scaled by annealing factor
	  total += mWeights[i];											//save sum of all weights
	  //std::cout<<"in CALCWEIGHTS -> mWeights = "<<mWeights[i]<<std::endl; 
	   if(i == 0 || mWeights[i] > best)									//find highest likelihood particle
	    {	
	      best = mWeights[i];
	      mBestParticle = i;
	    }
	}
	mWeights *= fpType(1.0) / total;										//normalize weights
	//std::cout<<"CALCWEIGHTS end mWeights = "<<mWeights<<std::endl;
}
void calcWeights_wrapper(ParticleFilterPipe<TrackingModelPipe> *pf, std::vector<vector<float>> *particles)
{
	pf->CalcWeights(*particles);
}

//Generate a set of inital particles
template<class T>
void ParticleFilterPipe<T>::InitializeParticles(int n)
{	
	mRnd.resize(n);														//initialize random number generators
	for(int i = 0; i < n; i++)
		mRnd[i].Seed(i * 2);
	mNParticles = n;
	Vectorf initVector = mModel->InitialState();						//get inital state vector 
	if(initVector.size() != mModel->StdDevs()[0].size() )
		std::cout << "Warning!! PF::Model initial Vector and stdDev vector are not equal lengths!" << std::endl;
	bool minValid = false;
	while(!minValid)
	{	mParticles.resize(n);
		for(int i = 0; i < n; i++)
		{	Vectorf &p = mParticles[i];									//create new particle with randomized initial value
			p = initVector;
			AddGaussianNoise(p, mModel->StdDevs()[0], mRnd[i]);			//distribute particles randomly about the initial vector
		}
		CalcWeights(mParticles);								//calculate initial weights and remove any particles that are invalid by model prior
		//call(calcWeights_wrapper, this, &mParticles);
		minValid = (int)mParticles.size() >= mMinParticles;				//repeat until minimum number of valid particles is met
		if(!minValid)							
			std::cout << "Warning : initial particle set does not meet minimum number of particles. Resampling.." << std::endl;
	}
	mCdf.resize(n);														//allocate space 
	mSamples.resize(n);
	mBins.resize(n);
	mInitialized = true;
}

void initializeParticles_wrapper(ParticleFilterPipe<TrackingModelPipe> *pf, int n)
{
	pf->InitializeParticles(n);
}


template<class T>
void ParticleFilterPipe<T>::GenerateNewParticles(int k)
{
  int p = 0;
  mNewParticles.resize(mNParticles);
  mIndex.resize(mNParticles);
  for(int i = 0; i < (int)mBins.size(); i++)
    for(uint j = 0; j < mBins[i]; j++)
      mIndex[p++] = i;

  blocked_range<int> range(0, mNParticles, WORKUNIT_SIZE_NEWPARTICLES);
  parallel_for(range, DoGenerateNewParticlesPipe<T>(k, getModel(), mNewParticles, mParticles, mRnd, mIndex));
  //DoGenerateNewParticlesPipe<T>(k, getModel(), mNewParticles, mParticles, mRnd, mIndex)(range);
  
}

void GenerateNewParticles_wrapper(ParticleFilterPipe<TrackingModelPipe> *pf, int k)
{
	pf->GenerateNewParticles(k);
}


//Particle filter update (model and observation updates must be called first)  
template<class T>
bool ParticleFilterPipe<T>::Update(fpType timeval)							//weights have already been computed from previous step or initialization
{						
	if(!mInitialized)														//check for proper initialization
	{	std::cout << "Update Error : Particles not initialized" << std::endl; 
		return false;
	}	
	int res;
	//mModel->GetObservation(timeval, &res);
	//call(GetObservation_wrapper, mModel, timeval, &res);
	//if(!res)
	{//	std::cout << "Update Error : Model observation failed for time : " << timeval << std::endl;
		//return false;
	}
	//cout<< "in Update model->stddevs.size = "<<(int)mModel->StdDevs().size()<<endl;
	for(int k = (int)mModel->StdDevs().size() - 1; k >= 0 ; k--)			//loop over all annealing steps starting with highest
	{	CalcCDF(mWeights, mCdf);						//Monte Carlo re-sampling 
		Resample(mCdf, mBins, mSamples, mNParticles);		
		bool minValid = false;
		while(!minValid)
		{	GenerateNewParticles(k);
			//call(GenerateNewParticles_wrapper, this, k);				//parallel_for inside
			CalcWeights(mNewParticles);									//calculate particle weights and remove any invalid particles
			//call(calcWeights_wrapper, this, &mNewParticles);		//parallel_for inside
			minValid = (int)mNewParticles.size() >= mMinParticles;		//repeat if not enough valid particles
			if(!minValid) 
				std::cout << "Not enough valid particles - Resampling!!!" << std::endl;
		}
		mParticles = mNewParticles;						//save new particle set
	}
	return true;
}

//write a given pose to a stream
template<class T>
void ParticleFilterPipe<T>::WritePose(std::ostream &f, std::vector<float> &pose)
{	
	//cout << "in WritePose pose.size = " << pose.size() << " pose[i] = " << pose[0] << endl;
	for(int i = 0; i < (int)pose.size(); i++)
	{
		mPoseOutFile << pose[i] << " ";
		f << pose[i] << " ";
		//cout<< "pose["<<i<<"] = "<<pose[i]<<endl;
	}
	mPoseOutFile << std::endl;
	f << std::endl;
}

void update_wrapper(ParticleFilterPipe<TrackingModelPipe> *pf, float currFrame)
{
	pf->Update((float)currFrame);
}

template<class T>
static void ParticleFilterPipe<T>::stage_operator(ParticleFilterPipe<T> * pf, ImageSetToken *images)
{
        pf->mModel->SetObservation(images);
        pf->Update(0);                                                                                                                              //call particle filter update
        std::vector<float> estimate;                                                                                    //expected pose from particle distribution
        pf->ParticleFilter<T>::Estimate(estimate);                                                                                                          //get average pose of the particle distribution
        pf->WritePose(pf->mPoseOutFile, estimate);
        if(pf->mOutputBMP)
                pf->mModel->OutputBMP(estimate, pf->mFrame++);                                                          //save output bitmap file

        delete images;
}



void stage2(/*obj::inoutdep<ParticleFilterPipe<TrackingModelPipe>*> */ParticleFilterPipe<TrackingModelPipe>* pf, obj::indep<ImageSetToken*> /* ImageSetToken**/ tokenObj, string *output, int currFrame, obj::inoutdep<void>)
{
    leaf_call( ParticleFilterPipe<TrackingModelPipe>::stage_operator, pf, (ImageSetToken*)tokenObj );
}
#if 0
{
	//stimer_t tmr_stage2;
	//stimer_t update_timer;
	//stimer_tick(&tmr_stage2);
	ofstream outputFileAvg((output)->c_str());
	ParticleFilterPipe<TrackingModelPipe>* PF = (ParticleFilterPipe<TrackingModelPipe>*) pf;
	ImageSetToken *token = (ImageSetToken *)tokenObj;

	TrackingModelPipe* model = (TrackingModelPipe*)PF->getmModel();
	leaf_call(setObservation_wrapper, model, token);
	//model->SetObservation(token);		
	//stimer_tick(&update_timer);
	PF->Update(0);													//**************************************
	//stimer_tuck(&update_timer, "Update took");
	//call(update_wrapper, PF, (float)currFrame);
	//leaf_call(setObservation_wrapper, model, token);
	std::vector<float> estimate;
	PF->Estimate(estimate);									//get average pose of the particle distribution
//	tokenObj = token;
	PF->WritePose(outputFileAvg, estimate);
	
	if(PF->getmOutputBMP())
		model->OutputBMP(estimate, currFrame);
	delete token;
	//stimer_tuck(&tmr_stage2, "Stage 2 took");
	return;
}
#endif


/*
int stage2_wrapper(obj::indep<ParticleFilterPipe<TrackingModelPipe>*> pf, obj::inoutdep<ImageSetToken*> token, string output, int currFrame)
{
	//TrackingModelPipe* model = (TrackingModelPipe*)pf->getmModel();
	//model.setObservation((ImageSetToken)*token);
	ImageSetToken* tk = (ImageSetToken*)token;
	//leaf_call(stage2, (ParticleFilterPipe<TrackingModelPipe>*)pf, (ImageSetToken*)tk);
	//stage2((ParticleFilterPipe<TrackingModelPipe>*)pf, (ImageSetToken*)tk, &output, currFrame);
}
*/

#endif
