#ifndef _PARTICLE_FILTER_H_
#define _PARTICLE_FILTER_H_

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
#include <thread>
#include <mutex>

#include "Particle.h"


namespace ORB_SLAM2
{
#define PROJECTION_GRID_COL 320
#define PROJECTION_GRID_ROW 240
#define PROJECTION_GRID_SIZE 2

	class Frame;
	class RTMatrix;
	class Particle;

	class ProjectionNode
	{
		public:
			ProjectionNode() :bNoMore(false), GlobalSelection(6,0){}
			
			std::vector<cv::Point3f> vMapPointPositions;
			std::vector<cv::Point2f> vProjections;
			std::vector<int> vObs;
			std::vector<int> vMapPointIndices;
			cv::Point2i UL, UR, BL, BR;
			std::list<ProjectionNode>::iterator lit;
			void DivideNode(ProjectionNode &n1, ProjectionNode &n2, ProjectionNode &n3, ProjectionNode &n4);

			vector<int> GlobalSelection;

			bool bNoMore;
	};

class ParticleFilter
{
public:
	int nFrameIdx;
	std::vector<Particle*> My_PF;

	RTMatrix preResultOfFiltering;//output result of last farme
	RTMatrix ResaultOfFiltering;//output result
	bool bError;

	vector<float> Indices;

	double Sum_Weight;
	unsigned int MaxWeightParticleID;
	unsigned int MaxWeightParticleOutlier;
	double MaxWeightParticleTotalDiff;
	vector<vector<Particle*>> vvppGroups;
	vector<RTMatrix> vrtGroupCenters;

	//for clustering/grouping
	vector<vector<Particle*>> vvpCommonGroups;
	vector<RTMatrix> vrtGroupCenter;
	vector<float> GroupLikelihoods;
	vector<float> vAngleRadius;
	vector<float> vTransRadius;
	float MaxGroupL;
	int MaxGroupLIdx;
	//for clustering/grouping

	// projection range
	vector<pair<cv::Point2f, cv::Point2f>> vpairRectangle;
	vector<pair<int, int>> vpairDisAndIdx;
	// projection range

	// Projection grid
	Particle MaxWeightParticle;
	vector<MapPoint*> vmpProjectionGrid[PROJECTION_GRID_COL][PROJECTION_GRID_ROW];
	vector<size_t> vstMPIndices[PROJECTION_GRID_COL][PROJECTION_GRID_ROW];
	vector<MapPoint*> vmpLocalMapPoints;
	vector<size_t> vstLocalMpIndices;
	bool PosInGrid(cv::Point2f &Projection, int& PosX, int& PosY);
	void ProjectionGrid(Frame* CurrentFrame, const Frame* LastFrame);
	void ProjectionGrid2(Frame* CurrentFrame, const Frame* LastFrame);
	void DistributeProjectionOctTree(const Frame* LastFrame, int &MapPointNumber, int& ProjectionTime);
	bool Project(const cv::Mat & Rcw, const cv::Mat & tcw, cv::Mat &MPPosition, vector<float> &Indices, cv::Point2f& Projection);

	vector<cv::Point3f> vx3Dw;
	vector<int> viMP;
	vector<int> viOb;

	vector<bool> vbNecessary;
	// Projection grid

	//for debug
	ofstream OutputGroupL;
	ofstream GroupedParticles;
	ofstream GroupResult;
	ofstream TimeEachFrame;
	ofstream FeaturesPosition;
	ofstream ProMatchNumber;
	ofstream OutputAbnormalMPs;
	//for debug

	// Time consuming
	vector<float> vTimesPrediction;
	vector<float> vTimesObservation;
	vector<float> vTimesReweight;
	vector<float> vTimesResampling;
	vector<float> vTimesEstimation;

	vector<float> vTimesObProjection;
	vector<float> vTimesObMatching;
	vector<float> vTimesObLikelihood;

	vector<float> vTimeObservation1;
	vector<float> vTimeObservation2;
	vector<float> vTimeObservation3;
	vector<float> vTimeObservation4;
	vector<float> vTimeObservation5;


	double ttrack1;
	double ttrack2;
	double ttrack3;
	double ttrack4;
	double ttrack5;

	double tObProjection;
	double tObMatching;
	double tObLikelihood;

	double ttrackO1;
	double ttrackO2;
	double ttrackO3;
	double ttrackO4;
	double ttrackO5;

	ParticleFilter();

	void ParticleFilterInitialization(RTMatrix Initialzation_Pose, vector<float> Initial_Speed);
	void PFTrackCurrentPose(Frame* CurrentFrame, const Frame* LastFrame);
	void Prediction();
	void Prediction2();//predict particles to two groups, one is large rotation, one is large transition
	void Observation(Frame& pFrame);

	void Observation123(Frame* pFrame);//for test the time consuming in each step of observation
	void Reweight();
	void Resampling();
	void Estimation();


	// for rectangle idea
	void PFTrackCurrentPoseRect(Frame* CurrentFrame, Frame* LastFrame);
	void ObservationRect(Frame* CurrentFrame, const Frame* LastFrame);
	void Projection(Frame* CurrentFrame, const Frame* LastFrame);
	int MatchingByProjection(Frame* CurrentFrame, const Frame *LastFrame, int EnlargeTh);
	int DescriptorDistance(const cv::Mat &a, const cv::Mat &b);
	void ComputeThreeMaxima(vector<int>* histo, const int L, int &ind1, int &ind2, int &ind3);
	int nmatches;
	// for rectangle idea

	// for small range matching idea
	void Observation2(Frame* CurrentFrame, const Frame* LastFrame);
	int MatchingByProjection2(Frame* CurrentFrame, const Frame *LastFrame, int EnlargeTh);
	void GetMapPointsInGrid(int& PosX, int& PosY, vector<size_t>& vstMapPointIndices, vector<MapPoint*>& vmpMapPoints);
	// for small range matching idea

	// for debug
	void OutputParticleCloud();
	void OutputOutliers();
	void OutputRectangles(Frame& LastFrame);
	void OutputAllParticles(vector<vector<Particle*>> &vvpOutput);
	void OutputPrediction();
	void OutputMatching(Frame& CurrentFrame);
	void OutputMPs(vector<cv::Point3f> vmMapPoints, vector<int> viObs);
	void OutputProjections(Frame& CurrentFrame, const Frame& LastFrame);
	void OutputProjections(Frame& LastFrame, vector<cv::Point3f> &vp3fNeMapPoints);
	void OutputLastFrameMaxWeightPose();
	void OutputAllMPs(Frame* LastFrame);
	void OutputDescriptors(Frame* CurrentFrame, Frame* LastFrame);
	void OutputLastFrameProjections(Frame* LastFrame);
	vector<int> KeysToRect;
	// for debug

	//for common clustering
	void SelectiveEstimation();
	void CommonGrouping();// from Deng's first proposal, seems does not work :)
	vector<float> RTMatrixDistance(RTMatrix &Matrix1, RTMatrix &Matrix2);
	RTMatrix UpdateGroupCenter(RTMatrix &GroupCenter, unsigned int &GroupIdx, RTMatrix &NewParticlePose);
	RTMatrix CalculateWeigtedCenter(vector<Particle*> Group);
	void CalculateRaidius();
	float CommonGroupingLikelihood(vector<Particle*> Group, RTMatrix WeightedCenter, float AngleRadius, float TransRadius);

	//for common clustering

	void Grouping();

	void Reset();
	//~Tracking();

	std::mutex mPFPos;

};

}//namespace ORB_SLAM2


#endif