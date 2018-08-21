#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include <opencv2/highgui/highgui.hpp>
#include <stdio.h>
#include "opencv2/core/core.hpp"  
#include "highgui.h"  
#include "opencv2/imgproc/imgproc.hpp"  
#include "opencv2/features2d/features2d.hpp"  
#include "opencv2/nonfree/nonfree.hpp"  
#include "opencv2/legacy/legacy.hpp"  
#include "opencv2/opencv.hpp"
#include "Thirdparty/DBoW2/DUtils/Random.h"

#include "RTMatrix.h"
#include "MapPoint.h"
#include "Frame.h"

#include <Eigen/StdVector>
#include <Eigen/Dense>

#include "Thirdparty/g2o/g2o/stuff/sampler.h"
#include "Thirdparty/g2o/g2o/core/factory.h"

namespace ORB_SLAM2
{
	class Frame;
	class RTMatrix;
	class Particle;
	//class GaussianSampler;
class Particle
{
public:
	RTMatrix Pose;
	RTMatrix Speedvt;// speed*frame time
	cv::Mat Rcw;
	cv::Mat tcw;
	std::vector<float> ParticleSpeed;
	std::vector<bool> pOutlier;
	std::vector<float> ParticleAcceleration;

	float Likelihood, Weight, WeightSum;
	float SumDiff;
	float SumDiffX;
	float SumDiffY;

	int OutlierCounter;
	int InlierCounter;

	// for time consuming test
	vector<cv::Mat> vMRcw;
	vector<cv::Mat> vMtcw;
	vector<cv::Mat> vMx3Dw;
	vector<cv::Mat> vMx3Dc;
	
	// for time consuming test

	// for rectangle idea
	void Projection(Frame& CurrentpFrame,Frame& pLastFrame, vector<pair<cv::Point2f,cv::Point2f>> &Rectangle);
	vector<cv::Point2f> vp2fProjections;
	void ObservationRect(Frame& CurrentFrame, vector<int> KeysToRect);
	// for rectangle idea

	// for Projection grid idea
	void Projection(vector<float> &Indices, const vector<cv::Point3f> &x3Dw, vector<pair<cv::Point2f, cv::Point2f>> &rectangle, vector<int> &ProjectionIdx, const int &N);
	void Projection2(vector<float> &Indices, const vector<cv::Point3f> &x3Dw, vector<int> &ProjectionIdx, const int &N);
	bool ProjectionOne(vector<float> &Indices, const cv::Point3f *x3Dw, int &ProjectionIdx);

	bool ProjectionOne(vector<float> &Indices, const cv::Point3f *x3Dw, vector<pair<cv::Point2f, cv::Point2f>> &rectangle, int &ProjectionIdx);
	bool ProjectionOneFixedPoint(vector<float> &Indices, const cv::Point3f *x3Dw, vector<pair<cv::Point2f, cv::Point2f>> &rectangle, int &ProjectionIdx);
	bool ProjectionOne1(vector<float> &Indices, const cv::Point3f &x3Dw, vector<pair<cv::Point2f, cv::Point2f>> &rectangle, int &ProjectionIdx);
	bool ProjectionOne2(vector<float> &Indices, const cv::Point3f &x3Dw, vector<pair<cv::Point2f, cv::Point2f>> &rectangle, int &ProjectionIdx);
	bool ProjectionOne3(vector<float> &Indices, const cv::Point3f &x3Dw, vector<pair<cv::Point2f, cv::Point2f>> &rectangle, int &ProjectionIdx);

	vector<cv::Mat> vx3Dc;
	vector<float> vxc, vyc, vinvzc;

	void ObservationRectGrid(Frame* CurrentFrame, vector<int> KeysToRect, vector<int> ProjectionWeight);
	void ObservationRectGrid(Frame* CurrentFrame, vector<int> &KeysToRect);
	// for Projection grid idea

	Particle();

	void ParticleInitialization(RTMatrix &InitializationPose, vector<float> &InitialSpeed);
	void Prediction();
	void PredictionMat();
	void PredictionRot();
	void PredictionTrans();
	void Observation(Frame& CurrentFrame);

	// for HW test
	void PiecewiseDivision(float &input, float& output);
	void FixedPointTrans(float &input, float&output);
	// for HW test


	// for time consuming test
	void Observation1(Frame* CurrentFrame);
	void Observation11(Frame* CurrentFrame);
	void Observation12(Frame* CurrentFrame);
	void Observation13(Frame* CurrentFrame);
	void Observation2(Frame* CurrentFrame);
	void Observation3(Frame* CurrentFrame);
	// for time consuming test

	void Clone(Particle *Clone_Source);
	Eigen::Matrix<double, 3, 3> toMatrix3d(const cv::Mat &cvMat3);

	void Reset();

	//~Particle();

};
}



#endif