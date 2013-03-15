//-------------------------------------------------------------
//      ____                        _      _
//     / ___|____ _   _ ____   ____| |__  | |
//    | |   / ___| | | |  _  \/ ___|  _  \| |
//    | |___| |  | |_| | | | | |___| | | ||_|
//     \____|_|  \_____|_| |_|\____|_| |_|(_) Media benchmarks
//                  
//	  2006, Intel Corporation, licensed under Apache 2.0 
//
//  file : TrackingModelSwan.h
//  author : 
//  description : Observation model for kinematic tree body 
//				  tracking threaded with TBB.
//				  
//  modified : 
//--------------------------------------------------------------

#ifndef TRACKINGMODELTBB_H
#define TRACKINGMODELTBB_H

#if defined(HAVE_CONFIG_H)
# include "config.h"
#endif

#include "TrackingModel.h"
#include "swan/wf_interface.h"
#include "TBBtypes.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include "swan/tbb_interface.h"
using namespace obj;
using namespace std;

//-----------------TIMERS---------------------

#include <sys/time.h>

typedef struct {
        struct  timeval start, end;
        float   diff;
} stimer_t;

void stimer_tick(stimer_t *timer);
float stimer_tuck(stimer_t *timer, const char *msg);

//--------------END-TIMERS---------------------



class TrackingModelPipe : public TrackingModel {

protected:
	
	unsigned int mCurFrame;			//current frame to process
	unsigned int mNumFrames;		//total frames to be processed

	//Generate an edge map from the original camera image - threaded
	virtual void CreateEdgeMap(FlexImage8u &src, FlexImage8u &dst);


public:
	TrackingModelPipe() : mCurFrame(0), mNumFrames(0) {};
	virtual ~TrackingModelPipe() {}

	//sets
	void SetNumFrames(unsigned int frames) {mNumFrames = frames; };

	//load and process images - overloaded for future threading
	bool GetObservation(float timeval); // int* res);

	//give the model object the processed images
	void SetObservation(ImageSetToken* token);
	

	void incrmCurFrame() 	{ mCurFrame++; };
	unsigned int getmCurFrame()	{ return mCurFrame; };
	unsigned int getmNumFrames()	{ return mNumFrames; };
	
	void SetNumThreads(int n);
	bool Initialize(const string &path, int cameras, int layers);
	//generate processed images for pipeline stage - these get passed to the next pipe stage defined by ParticleFilterTBB.h
	//void *operator()(void *inToken);

};
void ComputeEdgeMapsPipe(FlexImage8u *src, FlexImage8u *dst);
FlexImage8u* GradientMagThresholdPipe(FlexImage8u *src, float* threshold)__attribute__((always_inline));
void GaussianBlurPipe(FlexImage8u *src, FlexImage8u *dst);
template<class T, class T2>
bool FlexFilterColumnVPipe(FlexImage<T,1> &src, FlexImage<T,1> &dst, T2 *kernel, int kernelSize, bool allocate = true);
template<class T, class T2>
bool FlexFilterRowVPipe(FlexImage<T,1> &src, FlexImage<T,1> &dst, T2 *kernel, int kernelSize, bool allocate = true);

void setObservation_wrapper(TrackingModelPipe *model, ImageSetToken *token);

void stage1(/*indep<TrackingModelPipe *>*/TrackingModelPipe * modelObj, outdep<ImageSetToken*> tokenObj, int currFrame, inoutdep<void>);
//int stage1_wrapper(TrackingModelPipe * modelObj, outdep<ImageSetToken*> tokenObj, int currFrame);

// void SetMFGMaps_wrapper(TrackingModelPipe* model, std::vector<BinaryImage> FGMaps);
void GetObservation_wrapper(TrackingModelPipe* model, float timeval, int *res)__attribute__((always_inline));
void SetNumThreads_wrapper(TrackingModelPipe *model, int n);
void SetNumFrames_wrapper(TrackingModelPipe* model, int frames);
void initialize_wrapper(TrackingModelPipe *model, string *path, int cameras, int layers, int *res);
#endif //TRACKINGMODELTBB_H
