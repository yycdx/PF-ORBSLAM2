#include "ParticleFilter.h"

using namespace std;

namespace ORB_SLAM2
{
	const int ParticleNumber = 500;
	const float CommonGroupingAngleTH = 0.03;//thresholds for common grouping
	const float CommonGroupingTransTH = 0.03;
	const int ZeroGroupLTH = 10;
	const int Resolution_x = 640;
	const int Resolution_y = 480;
	const int TH_HIGH = 100;
	const int HISTO_LENGTH = 30;
	const int fps = 762;
	const int GlobalSelctionLayer = 2;
	const int FASTNumber = 1500;
	const float SearchRange = 8.0;
	const float GridRange = 3.0;

	const vector<int> Around9X = { -1, 0, 1, -1, 0, 1, -1, 0, 1 };
	const vector<int> Around9Y = { -1, -1, -1, 0, 0, 0, 1, 1, 1 };

	const float InvGridX = (float)1 / 16;
	const float InvGridY = (float)1 / 16;

	const int ProjectionN = 250;

	ParticleFilter::ParticleFilter()//Number of particles
	{
		Sum_Weight = 0.0;
		MaxWeightParticleID = 0;
		MaxWeightParticleOutlier = 0;
		nFrameIdx = 0;
		MaxWeightParticleTotalDiff = 0;

		bError = false;

		MaxGroupL = 0;
		MaxGroupLIdx = 0;

		OutputGroupL.open("D:\\analysis\\GroupL\\GroupLikelihood.txt");
		GroupResult.open("D:\\analysis\\GroupL\\ClusteringResult.txt");
		TimeEachFrame.open("D:\\analysis\\projection time vs feature number\\TimeVSFeature.txt");
		ProMatchNumber.open("D:\\analysis\\Projection and matching number\\ProjectionAndMatchingNumber.txt");
		OutputAbnormalMPs.open("D:\\analysis\\AbnormalMPS\\abnormalMPS.txt");

		My_PF = vector<Particle*>(ParticleNumber, static_cast<Particle*>(NULL));

		for (int i = 0; i < PROJECTION_GRID_COL; i++)
			for (int j = 0; j < PROJECTION_GRID_ROW; j++)
			{
				vmpProjectionGrid[i][j].reserve(5);
				vstMPIndices[i][j].reserve(5);
			}


		vx3Dw.reserve(ProjectionN * 2);
		viMP.reserve(ProjectionN * 2);
		viOb.reserve(ProjectionN * 2);

		vpairRectangle.reserve(FASTNumber + 100);
		vbNecessary.reserve(FASTNumber + 100);

		vmpLocalMapPoints.reserve(ProjectionN);
		vstLocalMpIndices.reserve(ProjectionN);
	}

	void ParticleFilter::ParticleFilterInitialization(RTMatrix Initialization_Pose, vector<float> Initial_Speed)
	{
		for (size_t i = 0, i_end = ParticleNumber; i < i_end; i++){
			My_PF[i] = new Particle();
			My_PF[i]->ParticleInitialization(Initialization_Pose, Initial_Speed);
			//cout << My_PF[i].Pose.Idx[0] << "\t" << My_PF[i].Pose.Idx[1] << "\t" << My_PF[i].Pose.Idx[2] << "\t" << My_PF[i].Pose.Idx[3] << "\t" << My_PF[i].Pose.Idx[4] << "\t" << My_PF[i].Pose.Idx[5] << endl;
		}
		MaxWeightParticle.ParticleInitialization(Initialization_Pose, Initial_Speed);
		cout << "Projections capacity: " << My_PF[0]->vp2fProjections.capacity() << "\t" << "Projections size: " << My_PF[0]->vp2fProjections.size() << endl;
		Sum_Weight = 0.0;
		preResultOfFiltering.Clone(Initialization_Pose);
	}

	void ParticleFilter::PFTrackCurrentPose(Frame* CurrentFrame, const Frame* LastFrame)
	{
		for (int i = 0; i < PROJECTION_GRID_COL; i++)
			for (int j = 0; j < PROJECTION_GRID_ROW; j++)
			{
				vmpProjectionGrid[i][j].clear();
				vstMPIndices[i][j].clear();
			}

		preResultOfFiltering.Clone(ResaultOfFiltering);
		ResaultOfFiltering.Clear();

		Sum_Weight = 0.0;
		MaxWeightParticleID = 0;
		MaxWeightParticleOutlier = 0;
		MaxWeightParticleTotalDiff = 0.0;

		std::chrono::monotonic_clock::time_point t1 = std::chrono::monotonic_clock::now();
		Prediction();
		std::chrono::monotonic_clock::time_point t2 = std::chrono::monotonic_clock::now();
		Observation2(CurrentFrame, LastFrame);
		std::chrono::monotonic_clock::time_point t3 = std::chrono::monotonic_clock::now();
		Reweight();
		std::chrono::monotonic_clock::time_point t4 = std::chrono::monotonic_clock::now();
		Resampling();
		std::chrono::monotonic_clock::time_point t5 = std::chrono::monotonic_clock::now();

		MaxWeightParticleID = 0;
		MaxWeightParticleOutlier = 0;
		MaxWeightParticleTotalDiff = 0.0;
		Estimation();

		CurrentFrame->SetPose(ResaultOfFiltering.RT);
		CurrentFrame->mvbOutlier = My_PF[MaxWeightParticleID]->pOutlier;
		MaxWeightParticleOutlier = My_PF[MaxWeightParticleID]->OutlierCounter;
		MaxWeightParticle.Pose.Clone(ResaultOfFiltering);// for projection grid

		std::chrono::monotonic_clock::time_point t6 = std::chrono::monotonic_clock::now();

		ttrack1 = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1).count();
		ttrack2 = std::chrono::duration_cast<std::chrono::duration<double> >(t3 - t2).count();
		ttrack3 = std::chrono::duration_cast<std::chrono::duration<double> >(t4 - t3).count();
		ttrack4 = std::chrono::duration_cast<std::chrono::duration<double> >(t5 - t4).count();
		ttrack5 = std::chrono::duration_cast<std::chrono::duration<double> >(t6 - t5).count();

		vTimesPrediction.push_back(ttrack1);
		vTimesObservation.push_back(ttrack2);
		vTimesReweight.push_back(ttrack3);
		vTimesResampling.push_back(ttrack4);
		vTimesEstimation.push_back(ttrack5);

		vTimesObProjection.push_back(tObProjection);
		vTimesObMatching.push_back(tObMatching);
		vTimesObLikelihood.push_back(tObLikelihood);

		nFrameIdx++;
		/*if (nFrameIdx % 100 == 0)
		{
		OutputParticleCloud();
		OutputOutliers();
		}*/


		Reset();
		//cout << "Reset finished" << endl;
	}

	void ParticleFilter::PFTrackCurrentPoseRect(Frame* CurrentFrame, Frame* LastFrame)
	{
		ResaultOfFiltering.Clear();

		Sum_Weight = 0.0;
		MaxWeightParticleID = 0;
		MaxWeightParticleOutlier = 0;
		MaxWeightParticleTotalDiff = 0.0;
		
		std::chrono::monotonic_clock::time_point t1 = std::chrono::monotonic_clock::now();
		Prediction();
		std::chrono::monotonic_clock::time_point t2 = std::chrono::monotonic_clock::now();
		ObservationRect(CurrentFrame, LastFrame);
		std::chrono::monotonic_clock::time_point t3 = std::chrono::monotonic_clock::now();
		Reweight();
		std::chrono::monotonic_clock::time_point t4 = std::chrono::monotonic_clock::now();
		Resampling();
		std::chrono::monotonic_clock::time_point t5 = std::chrono::monotonic_clock::now();

		MaxWeightParticleID = 0;
		MaxWeightParticleOutlier = 0;
		MaxWeightParticleTotalDiff = 0.0;
		Estimation();

		CurrentFrame->SetPose(ResaultOfFiltering.RT);
		CurrentFrame->mvbOutlier = My_PF[MaxWeightParticleID]->pOutlier;
		MaxWeightParticleOutlier = My_PF[MaxWeightParticleID]->OutlierCounter;
		MaxWeightParticle.Pose.Clone(ResaultOfFiltering);// for projection grid

		std::chrono::monotonic_clock::time_point t6 = std::chrono::monotonic_clock::now();

		ttrack1 = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1).count();
		ttrack2 = std::chrono::duration_cast<std::chrono::duration<double> >(t3 - t2).count();
		ttrack3 = std::chrono::duration_cast<std::chrono::duration<double> >(t4 - t3).count();
		ttrack4 = std::chrono::duration_cast<std::chrono::duration<double> >(t5 - t4).count();
		ttrack5 = std::chrono::duration_cast<std::chrono::duration<double> >(t6 - t5).count();

		vTimesPrediction.push_back(ttrack1);
		vTimesObservation.push_back(ttrack2);
		vTimesReweight.push_back(ttrack3);
		vTimesResampling.push_back(ttrack4);
		vTimesEstimation.push_back(ttrack5);

		vTimesObProjection.push_back(tObProjection);
		vTimesObMatching.push_back(tObMatching);
		vTimesObLikelihood.push_back(tObLikelihood);
		preResultOfFiltering.Clone(ResaultOfFiltering);
		nFrameIdx++;
		if (nFrameIdx % 100 == 0)
		{
			OutputParticleCloud();
			OutputLastFrameMaxWeightPose();
			OutputAllMPs(LastFrame);
			OutputDescriptors(CurrentFrame, LastFrame);
			OutputRectangles(*LastFrame);
			OutputLastFrameProjections(LastFrame);
		}


		Reset();
		//cout << "Reset finished" << endl;
	}

	void ParticleFilter::Prediction()
	{
		for (size_t i = 0, i_end = ParticleNumber; i < i_end; i++)
		{
			//cout << "Prediction of particle " << i << " : " << endl;
			My_PF[i]->PredictionMat();
		}
	}

	void ParticleFilter::ObservationRect(Frame* CurrentFrame, const Frame* LastFrame)
	{
		if (Indices.empty())
		{
			const float fx = CurrentFrame->fx;
			const float fy = CurrentFrame->fy;
			const float cx = CurrentFrame->cx;
			const float cy = CurrentFrame->cy;
			const float mnMinX = CurrentFrame->mnMinX;
			const float mnMinY = CurrentFrame->mnMinY;
			const float mnMaxX = CurrentFrame->mnMaxX;
			const float mnMaxY = CurrentFrame->mnMaxY;

			Indices = { fx, fy, cx, cy, mnMinX, mnMinY, mnMaxX, mnMaxY };
		}

		std::chrono::monotonic_clock::time_point tO0 = std::chrono::monotonic_clock::now();
		int N = LastFrame->N;
		vpairRectangle.clear();
		vpairRectangle.resize(N, pair<cv::Point2f, cv::Point2f>(cv::Point2f(Resolution_x, Resolution_y), cv::Point2f(0, 0)));

		ProjectionGrid(CurrentFrame, LastFrame);

		//OutputRectangles(LastFrame);
		//cout << "Projection finished"<<endl;
		std::chrono::monotonic_clock::time_point tO1 = std::chrono::monotonic_clock::now();
		tObProjection = std::chrono::duration_cast<std::chrono::duration<double> >(tO1 - tO0).count();
		TimeEachFrame << tObProjection << "\t";
		//OutputRectangles();
		/*int NoZeroRectanglesCount = 0;
		for (int i = 0; i < vpairRectangle.size(); i++)
		{
			if (vpairRectangle[i].second.x != 0 && vpairRectangle[i].second.y != 0)
			{
				//cout << vpairRectangle[i].first.x << "\t" << vpairRectangle[i].first.y << "\t" << vpairRectangle[i].second.x << "\t" << vpairRectangle[i].second.y << endl;
				NoZeroRectanglesCount++;
			}
		}*/
		//cout << "No zero rectangles number: " << NoZeroRectanglesCount << endl;
		int RangeTh = 0;
		nmatches = MatchingByProjection(CurrentFrame, LastFrame, RangeTh);
		if (nmatches < 20)
		{
			cout << "nmatches is less than 20, enlarge rectangles and try again" << endl;
			fill(CurrentFrame->mvpMapPoints.begin(),CurrentFrame->mvpMapPoints.end(),static_cast<MapPoint*>(NULL));
			nmatches = MatchingByProjection(CurrentFrame, LastFrame, RangeTh*2);
		}
		//cout<<"MatchingByProjection finished"<<endl;
		

		std::chrono::monotonic_clock::time_point tO2 = std::chrono::monotonic_clock::now();
		tObMatching= std::chrono::duration_cast<std::chrono::duration<double> >(tO2 - tO1).count();
		//OutputMatching(CurrentFrame);
		TimeEachFrame << tObMatching << "\t";
		//cout << "ObservationRect: after MatchingByProjection, nmatches = " << nmatches << endl;
		ProMatchNumber << nmatches << "\t" << tObMatching << "\t";

		for (size_t i = 0, i_end = ParticleNumber; i < i_end; i++)
			//My_PF[i].ObservationRect(CurrentFrame, KeysToRect);
			My_PF[i]->ObservationRectGrid(CurrentFrame, KeysToRect);
		//cout<<"Observation finished"<<endl;

		std::chrono::monotonic_clock::time_point tO3 = std::chrono::monotonic_clock::now();

		tObLikelihood = std::chrono::duration_cast<std::chrono::duration<double> >(tO3 - tO2).count();
		ProMatchNumber << tObLikelihood << endl;
		TimeEachFrame << tObLikelihood << endl;
	}

	void ParticleFilter::Observation2(Frame* CurrentFrame, const Frame* LastFrame)
	{
		if (Indices.empty())
		{
			const float fx = CurrentFrame->fx;
			const float fy = CurrentFrame->fy;
			const float cx = CurrentFrame->cx;
			const float cy = CurrentFrame->cy;
			const float mnMinX = CurrentFrame->mnMinX;
			const float mnMinY = CurrentFrame->mnMinY;
			const float mnMaxX = CurrentFrame->mnMaxX;
			const float mnMaxY = CurrentFrame->mnMaxY;

			Indices = { fx, fy, cx, cy, mnMinX, mnMinY, mnMaxX, mnMaxY };
		}

		std::chrono::monotonic_clock::time_point tO0 = std::chrono::monotonic_clock::now();
		int N = LastFrame->N;

		ProjectionGrid2(CurrentFrame, LastFrame);

		//OutputRectangles(LastFrame);
		//cout << "ProjectionGrid2 finished"<<endl;
		std::chrono::monotonic_clock::time_point tO1 = std::chrono::monotonic_clock::now();
		tObProjection = std::chrono::duration_cast<std::chrono::duration<double> >(tO1 - tO0).count();
		TimeEachFrame << tObProjection << "\t";
		//OutputRectangles();
		/*int NoZeroRectanglesCount = 0;
		for (int i = 0; i < vpairRectangle.size(); i++)
		{
		if (vpairRectangle[i].second.x != 0 && vpairRectangle[i].second.y != 0)
		{
		//cout << vpairRectangle[i].first.x << "\t" << vpairRectangle[i].first.y << "\t" << vpairRectangle[i].second.x << "\t" << vpairRectangle[i].second.y << endl;
		NoZeroRectanglesCount++;
		}
		}*/
		//cout << "No zero rectangles number: " << NoZeroRectanglesCount << endl;
		int RangeTh = 3;
		nmatches = MatchingByProjection2(CurrentFrame, LastFrame, RangeTh);
		if (nmatches < 20)
		{
			cout << "nmatches is less than 20, enlarge rectangles and try again" << endl;
			fill(CurrentFrame->mvpMapPoints.begin(), CurrentFrame->mvpMapPoints.end(), static_cast<MapPoint*>(NULL));
			nmatches = MatchingByProjection2(CurrentFrame, LastFrame, RangeTh * 2);
		}
		//cout<<"MatchingByProjection2 finished"<<endl;


		std::chrono::monotonic_clock::time_point tO2 = std::chrono::monotonic_clock::now();
		tObMatching = std::chrono::duration_cast<std::chrono::duration<double> >(tO2 - tO1).count();
		//OutputMatching(CurrentFrame);
		TimeEachFrame << tObMatching << "\t";
		//cout << "ObservationRect: after MatchingByProjection, nmatches = " << nmatches << endl;
		ProMatchNumber << nmatches << "\t" << tObMatching << "\t";

		for (size_t i = 0, i_end = ParticleNumber; i < i_end; i++)
			//My_PF[i].ObservationRect(CurrentFrame, KeysToRect);
			My_PF[i]->ObservationRectGrid(CurrentFrame, KeysToRect);
		//cout<<"Likelihood calculation finished"<<endl;

		std::chrono::monotonic_clock::time_point tO3 = std::chrono::monotonic_clock::now();

		tObLikelihood = std::chrono::duration_cast<std::chrono::duration<double> >(tO3 - tO2).count();
		ProMatchNumber << tObLikelihood << endl;
		TimeEachFrame << tObLikelihood << endl;
	}

	void ParticleFilter::Reweight()
	{
		for (size_t i = 0, i_end = ParticleNumber; i < i_end; i++)
			Sum_Weight += My_PF[i]->Likelihood;
		
		for (size_t i = 0, i_end = ParticleNumber; i < i_end; i++)
		{	
			if (Sum_Weight == 0)
				cout << "Error: Sum of weight is zero!" << endl;
			My_PF[i]->Weight = My_PF[i]->Likelihood / Sum_Weight;
			if (!isfinite(My_PF[i]->Weight))
				cout << "weight is not finite" << endl;
		}


		My_PF[0]->WeightSum = My_PF[0]->Weight;
		for (size_t i = 1, i_end = ParticleNumber; i < i_end; i++)
		{
			My_PF[i]->WeightSum = My_PF[i - 1]->WeightSum + My_PF[i]->Weight;
			//cout << My_PF[i].Pose.Idx[4] << "\t" << My_PF[i].Likelihood << "\t" << My_PF[i].Weight << "\t" << My_PF[i].WeightSum << endl;
		}
		//cout << "Sum weight: " << Sum_Weight << endl;
	}

	void ParticleFilter::Resampling()
	{
		for (size_t i = 0, i_end = ParticleNumber; i < i_end; i++){
			double randw = DUtils::Random::RandomValue(0.0, 1.0);
			//cout << randw << "\t";
			for (size_t j = 0, j_end = ParticleNumber - 1; j < j_end; j++){
				if (randw > My_PF[j]->WeightSum&&randw < My_PF[j + 1]->WeightSum){
					My_PF[i]->Clone(My_PF[j + 1]);
					break;
				}
				else
					continue;
			}
		}
	}

	void ParticleFilter::Estimation()
	{
		//OutputGroupL << nFrameIdx << "----------------------------------------------------" << endl;
		double MaxWeight = 0;
		fill(ResaultOfFiltering.Idx.begin(), ResaultOfFiltering.Idx.end(), 0.0);
		for (size_t i = 0, i_end = ParticleNumber; i < i_end; i++){
			if (My_PF[i]->Weight > MaxWeight){
				MaxWeight = My_PF[i]->Weight;
				MaxWeightParticleID = i;
			}
			for (size_t j = 0, j_end = 6; j < j_end; j++){
				if (!std::isfinite(My_PF[i]->Pose.Idx[j]))
				{
					cout << "PF: tracking result after estimation is not finite" << endl;
					bError = true;
				}

				ResaultOfFiltering.Idx[j] += My_PF[i]->Pose.Idx[j] * My_PF[i]->Weight;
				
			}
			//OutputGroupL << My_PF[i].Likelihood << "\t" << My_PF[i].Weight << endl;
			//cout << My_PF[i].Weight << "\t" << My_PF[i].Likelihood << My_PF[i].OutlierCounter << "\t" << My_PF[i].InlierCounter << "\t" << endl;
			//cout << My_PF[i].Pose.Idx[5] << "\t" << My_PF[i].Weight << "\t";
			//cout << ResaultOfFiltering.Idx[5] << endl;
		}
		//cout << endl;
		ResaultOfFiltering.Idx_to_Mat();
	}

	void ParticleFilter::Reset()
	{
		for (size_t i = 0; i < ParticleNumber; i++)
			My_PF[i]->Reset();
	}

	void ParticleFilter::OutputParticleCloud()
	{
		ofstream OutputParticleCloudDir;
		char OutputDir[80];
		sprintf(OutputDir, "D:\\analysis\\ParticleCloudPoint\\CloudPoint_%06d.txt", nFrameIdx);
		OutputParticleCloudDir.open(OutputDir);
		for (int i = 0; i < ParticleNumber; i++)
		{
			for (int j = 0; j < 3; j++)
				for (int k = 0; k < 4; k++)
					OutputParticleCloudDir << My_PF[i]->Pose.RT.at<float>(j, k) << "\t";
			//OutputParticleCloudDir << My_PF[i]->Weight << "\t";
			OutputParticleCloudDir << endl;
		}
	}

	void ParticleFilter::OutputOutliers()
	{
		ofstream Outliers;
		char OutputDir[80];
		sprintf(OutputDir, "D:\\analysis\\ParticleCloudPoint\\Outliers_%06d.txt", nFrameIdx);
		Outliers.open(OutputDir);
		for (int i = 0; i < ParticleNumber; i++)
			Outliers << My_PF[i]->OutlierCounter << "\t" << My_PF[i]->Weight << endl;
	}

	void ParticleFilter::ProjectionGrid(Frame* CurrentFrame, const Frame* LastFrame)
	{
		const int N = LastFrame->N;

		int NoNULLCount = 0;
		int MapPointNumber = 0;
		int ProjectionTime = 0;

		vx3Dw.clear();
		viMP.clear();
		viOb.clear();

		vbNecessary.clear();
		vbNecessary.resize(N, false);

		std::chrono::monotonic_clock::time_point tP0 = std::chrono::monotonic_clock::now();
		DistributeProjectionOctTree(LastFrame, MapPointNumber, ProjectionTime);
		std::chrono::monotonic_clock::time_point tP1 = std::chrono::monotonic_clock::now();

		std::chrono::monotonic_clock::time_point tP2 = std::chrono::monotonic_clock::now();
		for (int i_Particle = 0; i_Particle < ParticleNumber; ++i_Particle)
		{
			My_PF[i_Particle]->Projection(Indices, vx3Dw, vpairRectangle, viMP, N);
		}
		std::chrono::monotonic_clock::time_point tP3 = std::chrono::monotonic_clock::now();

		double tObOctTree = std::chrono::duration_cast<std::chrono::duration<double>>(tP1 - tP0).count();
		double tObProjection1 = std::chrono::duration_cast<std::chrono::duration<double>>(tP3 - tP2).count();

		ProMatchNumber
			<< N << "\t" << MapPointNumber << "\t" << ProjectionTime << "\t" << tObOctTree << "\t" << tObProjection1 /*<< "\t" << tObProjection2 << "\t" << tObProjection3*/ << "\t";
		TimeEachFrame << MapPointNumber << "\t";

	}

	void ParticleFilter::ProjectionGrid2(Frame* CurrentFrame, const Frame* LastFrame)
	{
		const int N = LastFrame->N;

		int NoNULLCount = 0;
		int MapPointNumber = 0;
		int ProjectionTime = 0;

		vx3Dw.clear();
		viMP.clear();
		viOb.clear();

		vbNecessary.clear();
		vbNecessary.resize(N, false);

		std::chrono::monotonic_clock::time_point tP0 = std::chrono::monotonic_clock::now();
		DistributeProjectionOctTree(LastFrame, MapPointNumber, ProjectionTime);
		std::chrono::monotonic_clock::time_point tP1 = std::chrono::monotonic_clock::now();

		std::chrono::monotonic_clock::time_point tP2 = std::chrono::monotonic_clock::now();
		for (int i_Particle = 0; i_Particle < ParticleNumber; ++i_Particle)
		{
			My_PF[i_Particle]->Projection2(Indices, vx3Dw, viMP, N);
		}
		std::chrono::monotonic_clock::time_point tP3 = std::chrono::monotonic_clock::now();

		double tObOctTree = std::chrono::duration_cast<std::chrono::duration<double>>(tP1 - tP0).count();
		double tObProjection1 = std::chrono::duration_cast<std::chrono::duration<double>>(tP3 - tP2).count();

		ProMatchNumber
			<< N << "\t" << MapPointNumber << "\t" << ProjectionTime << "\t" << tObOctTree << "\t" << tObProjection1 /*<< "\t" << tObProjection2 << "\t" << tObProjection3*/ << "\t";
		TimeEachFrame << MapPointNumber << "\t";


	}

	void ParticleFilter::DistributeProjectionOctTree(const Frame* LastFrame, int &MapPointNumber, int& ProjectionTime)
	{
		const cv::Mat Rcw = MaxWeightParticle.Pose.RT.rowRange(0, 3).colRange(0, 3);
		const cv::Mat tcw = MaxWeightParticle.Pose.RT.rowRange(0, 3).col(3);

		int N = LastFrame->N;
		vector<cv::Point3f> MapPointPositions;
		vector<cv::Point2f> Projections;
		vector<int> vObs;
		vector<int> vMapPointIndices;

		MapPointPositions.reserve(N);
		vObs.reserve(N);
		Projections.reserve(N);
		vMapPointIndices.reserve(N);

		int MapPointNumber0 = 0;
		int ProjectionTime0 = 0;

		{
			unique_lock<mutex> lock(MapPoint::mGlobalMutex);

			for (int i_MP = 0; i_MP < N; i_MP++)
			{
				MapPoint* pMP = LastFrame->mvpMapPoints[i_MP];
				if (pMP)
				{
					if (!LastFrame->mvbOutlier[i_MP])// Current mappoint exist in last frame
					{
						MapPointNumber0++;
						cv::Mat x3Dw = pMP->GetWorldPos();
						cv::Point2f Projection;
						if (Project(Rcw, tcw, x3Dw, Indices, Projection))
						{
							cv::Point3f NewNeMapPoint(x3Dw.at<float>(0), x3Dw.at<float>(1), x3Dw.at<float>(2));

							MapPointPositions.push_back(NewNeMapPoint);
							vObs.push_back(pMP->Observations());
							Projections.push_back(Projection);
							vMapPointIndices.push_back(i_MP);
						}
					}
				}
			}
		}

		const int nIni = round(Resolution_x / Resolution_y);

		const float hX = static_cast<float>(Resolution_x) / nIni;

		list<ProjectionNode> lNodes;
		vector<ProjectionNode*> vpIniNodes;
		vpIniNodes.resize(nIni);
		//initialize initial nodes
		for (int i = 0; i < nIni; i++)
		{
			ProjectionNode ni;
			ni.UL = cv::Point2i(hX*static_cast<float>(i), 0);
			ni.UR = cv::Point2i(hX*static_cast<float>(i + 1), 0);
			ni.BL = cv::Point2i(ni.UL.x, Resolution_y);
			ni.BR = cv::Point2i(ni.UR.x, Resolution_y);
			ni.vProjections.reserve(Projections.size());

			lNodes.push_back(ni);
			vpIniNodes[i] = &lNodes.back();
		}
		//add projections to initial nodes
		for (size_t i = 0; i < Projections.size(); i++)
		{
			const cv::Point2f &kp = Projections[i];
			const int Ob = vObs[i];
			const cv::Point3f &MP = MapPointPositions[i];
			const int MapPointIndex = vMapPointIndices[i];

			vpIniNodes[kp.x / hX]->vProjections.push_back(kp);
			vpIniNodes[kp.x / hX]->vObs.push_back(Ob);
			vpIniNodes[kp.x / hX]->vMapPointPositions.push_back(MP);
			vpIniNodes[kp.x / hX]->vMapPointIndices.push_back(MapPointIndex);
		}

		list<ProjectionNode>::iterator lit = lNodes.begin();

		while (lit != lNodes.end())
		{
			if (lit->vProjections.size() == 1)
			{
				lit->bNoMore = true;
				lit++;
			}
			else if (lit->vProjections.empty())
				lit = lNodes.erase(lit);
			else
				lit++;
		}

		bool bFinish = false;

		int iteration = 0;

		vector<pair<int, ProjectionNode*>> vSizeAndPointerToNode;
		vSizeAndPointerToNode.reserve(lNodes.size() * 4);

		while (!bFinish)
		{
			iteration++;

			int prevSize = lNodes.size();

			lit = lNodes.begin();

			int nToExpand = 0;

			vSizeAndPointerToNode.clear();

			while (lit != lNodes.end())
			{
				if (lit->bNoMore)
				{
					// If node only contains one point do not subdivide and continue
					lit++;
					continue;
				}
				else
				{
					// If more than one point, subdivide
					ProjectionNode n1, n2, n3, n4;
					lit->DivideNode(n1, n2, n3, n4);

					// Add childs if they contain points
					if (n1.vProjections.size()>0)
					{
						lNodes.push_front(n1);
						if (n1.vProjections.size()>1)
						{
							nToExpand++;
							vSizeAndPointerToNode.push_back(make_pair(n1.vProjections.size(), &lNodes.front()));
							lNodes.front().lit = lNodes.begin();
						}
					}
					if (n2.vProjections.size()>0)
					{
						lNodes.push_front(n2);
						if (n2.vProjections.size()>1)
						{
							nToExpand++;
							vSizeAndPointerToNode.push_back(make_pair(n2.vProjections.size(), &lNodes.front()));
							lNodes.front().lit = lNodes.begin();
						}
					}
					if (n3.vProjections.size()>0)
					{
						lNodes.push_front(n3);
						if (n3.vProjections.size()>1)
						{
							nToExpand++;
							vSizeAndPointerToNode.push_back(make_pair(n3.vProjections.size(), &lNodes.front()));
							lNodes.front().lit = lNodes.begin();
						}
					}
					if (n4.vProjections.size()>0)
					{
						lNodes.push_front(n4);
						if (n4.vProjections.size()>1)
						{
							nToExpand++;
							vSizeAndPointerToNode.push_back(make_pair(n4.vProjections.size(), &lNodes.front()));
							lNodes.front().lit = lNodes.begin();
						}
					}

					lit = lNodes.erase(lit);
					continue;
				}
			}

			if ((int)lNodes.size() >= ProjectionN || (int)lNodes.size() == prevSize)
				bFinish = true;
			else if (((int)lNodes.size() + nToExpand * 3) > N)
			{
				while (!bFinish)
				{
					prevSize = lNodes.size();

					vector<pair<int, ProjectionNode*> > vPrevSizeAndPointerToNode = vSizeAndPointerToNode;
					vSizeAndPointerToNode.clear();

					sort(vPrevSizeAndPointerToNode.begin(), vPrevSizeAndPointerToNode.end());
					for (int j = vPrevSizeAndPointerToNode.size() - 1; j >= 0; j--)
					{
						ProjectionNode n1, n2, n3, n4;
						vPrevSizeAndPointerToNode[j].second->DivideNode(n1, n2, n3, n4);

						// Add childs if they contain points
						if (n1.vProjections.size()>0)
						{
							lNodes.push_front(n1);
							if (n1.vProjections.size()>1)
							{
								vSizeAndPointerToNode.push_back(make_pair(n1.vProjections.size(), &lNodes.front()));
								lNodes.front().lit = lNodes.begin();
							}
						}
						if ((int)lNodes.size() >= N)
							break;
						if (n2.vProjections.size()>0)
						{
							lNodes.push_front(n2);
							if (n2.vProjections.size()>1)
							{
								vSizeAndPointerToNode.push_back(make_pair(n2.vProjections.size(), &lNodes.front()));
								lNodes.front().lit = lNodes.begin();
							}
						}
						if ((int)lNodes.size() >= N)
							break;
						if (n3.vProjections.size()>0)
						{
							lNodes.push_front(n3);
							if (n3.vProjections.size()>1)
							{
								vSizeAndPointerToNode.push_back(make_pair(n3.vProjections.size(), &lNodes.front()));
								lNodes.front().lit = lNodes.begin();
							}
						}
						if ((int)lNodes.size() >= N)
							break;
						if (n4.vProjections.size()>0)
						{
							lNodes.push_front(n4);
							if (n4.vProjections.size()>1)
							{
								vSizeAndPointerToNode.push_back(make_pair(n4.vProjections.size(), &lNodes.front()));
								lNodes.front().lit = lNodes.begin();
							}
						}

						lNodes.erase(vPrevSizeAndPointerToNode[j].second->lit);

						if ((int)lNodes.size() >= N)
							break;
					}

					if ((int)lNodes.size() >= N || (int)lNodes.size() == prevSize)
						bFinish = true;
				}
			}
		}

		/*vector<cv::Point3f> vResultMapPoints;
		vector<int> vResultIndices;
		vector<int> vResultObs;
		vResultMapPoints.reserve(2*ProjectionN);
		vResultIndices.reserve(2*ProjectionN);
		vResultObs.reserve(2*ProjectionN);*/
		int ProjectionCount = 0;
		for (list<ProjectionNode>::iterator lit = lNodes.begin(); lit != lNodes.end(); lit++)
		{
			vector<cv::Point3f> &vNodeMapPoints = lit->vMapPointPositions;
			vector<cv::Point2f> &vNodeProjection = lit->vProjections;
			vector<int> &vNodeObs = lit->vObs;
			vector<int> &vNodeIndices = lit->vMapPointIndices;

			cv::Point3f* pMapPoint = &vNodeMapPoints[0];
			cv::Point2f* pProjection = &vNodeProjection[0];
			int maxOb = lit->vObs[0];
			int maxIndex = lit->vMapPointIndices[0];

			for (size_t k = 1; k < vNodeProjection.size(); k++)
			{
				if (vNodeObs[k] > maxOb)
				{
					pMapPoint = &vNodeMapPoints[k];
					pProjection = &vNodeProjection[k];
					maxOb = lit->vObs[k];
					maxIndex = lit->vMapPointIndices[k];
				}
			}
			ProjectionCount++;
			vx3Dw.push_back(*pMapPoint);
			viMP.push_back(maxIndex);
			viOb.push_back(maxOb);
			vbNecessary[maxIndex] = true;
			if (ProjectionCount >= ProjectionN)
				break;
		}
		//vNeMapPoints = vResultMapPoints;
		//vNeMapPointIndices = vResultIndices;
		//vNeObs = vResultObs;
		ProjectionTime0 = vx3Dw.size();
		MapPointNumber = MapPointNumber0;
		ProjectionTime = ProjectionTime0;

		//cout << "Oct tree finish, MapPoin Number: " << MapPointNumber << "\t" << "ProjectionTimes: " << ProjectionTime << endl;


	}

	void ProjectionNode::DivideNode(ProjectionNode & n1, ProjectionNode & n2, ProjectionNode & n3, ProjectionNode & n4)
	{
		const int halfX = ceil(static_cast<float>(UR.x - UL.x) / 2);
		const int halfY = ceil(static_cast<float>(BR.y - UL.y) / 2);

		//Define boundaries of childs
		n1.UL = UL;
		n1.UR = cv::Point2i(UL.x + halfX, UL.y);
		n1.BL = cv::Point2i(UL.x, UL.y + halfY);
		n1.BR = cv::Point2i(UL.x + halfX, UL.y + halfY);
		n1.vProjections.reserve(vProjections.size());

		n2.UL = n1.UR;
		n2.UR = UR;
		n2.BL = n1.BR;
		n2.BR = cv::Point2i(UR.x, UL.y + halfY);
		n2.vProjections.reserve(vProjections.size());

		n3.UL = n1.BL;
		n3.UR = n1.BR;
		n3.BL = BL;
		n3.BR = cv::Point2i(n1.BR.x, BL.y);
		n3.vProjections.reserve(vProjections.size());

		n4.UL = n3.UR;
		n4.UR = n2.BR;
		n4.BL = n3.BR;
		n4.BR = BR;
		n4.vProjections.reserve(vProjections.size());

		n1.GlobalSelection = GlobalSelection;
		n2.GlobalSelection = GlobalSelection;
		n3.GlobalSelection = GlobalSelection;
		n4.GlobalSelection = GlobalSelection;

		//Associate points to childs
		for (size_t i = 0; i<vProjections.size(); i++)
		{
			const cv::Point3f &MP = vMapPointPositions[i];
			const cv::Point2f &kp = vProjections[i];
			const int Ob = vObs[i];
			const int Index = vMapPointIndices[i];
			if (kp.x<n1.UR.x)
			{
				if (kp.y < n1.BR.y)
				{
					n1.vMapPointPositions.push_back(MP);
					n1.vProjections.push_back(kp);
					n1.vObs.push_back(Ob);
					n1.vMapPointIndices.push_back(Index);
				}
				else
				{
					n3.vMapPointPositions.push_back(MP);
					n3.vProjections.push_back(kp);
					n3.vObs.push_back(Ob);
					n3.vMapPointIndices.push_back(Index);
				}

			}
			else if (kp.y < n1.BR.y)
			{
				n2.vMapPointPositions.push_back(MP);
				n2.vProjections.push_back(kp);
				n2.vObs.push_back(Ob);
				n2.vMapPointIndices.push_back(Index);
			}
			else
			{
				n4.vMapPointPositions.push_back(MP);
				n4.vProjections.push_back(kp);
				n4.vObs.push_back(Ob);
				n4.vMapPointIndices.push_back(Index);
			}

		}

		if (n1.vProjections.size() == 1)
			n1.bNoMore = true;
		if (n2.vProjections.size() == 1)
			n2.bNoMore = true;
		if (n3.vProjections.size() == 1)
			n3.bNoMore = true;
		if (n4.vProjections.size() == 1)
			n4.bNoMore = true;

	}

	bool ParticleFilter::Project(const cv::Mat & Rcw, const cv::Mat & tcw, cv::Mat &x3Dw, vector<float> &Indices, cv::Point2f & Projection)
	{
		/*cv::Mat x3Dc = cv::Mat(3, 1, CV_32F, 0.0);
		for (int i = 0; i < 3; i++)
		{
			x3Dc.at<float>(i) =
				Rcw.at<float>(i, 0)*x3Dw.at<float>(0) +
				Rcw.at<float>(i, 1)*x3Dw.at<float>(1) +
				Rcw.at<float>(i, 2)*x3Dw.at<float>(2) +
				tcw.at<float>(i);
		}

		const float xc = x3Dc.at<float>(0);
		const float yc = x3Dc.at<float>(1);
		const float invzc = 1.0 / x3Dc.at<float>(2);*/

		const float xc = Rcw.at<float>(0, 0)*x3Dw.at<float>(0) + Rcw.at<float>(0, 1)*x3Dw.at<float>(1) + Rcw.at<float>(0, 2)*x3Dw.at<float>(2) + tcw.at<float>(0);
		const float yc = Rcw.at<float>(1, 0)*x3Dw.at<float>(0) + Rcw.at<float>(1, 1)*x3Dw.at<float>(1) + Rcw.at<float>(1, 2)*x3Dw.at<float>(2) + tcw.at<float>(1);
		const float invzc = 1.0 / (Rcw.at<float>(2, 0)*x3Dw.at<float>(0) + Rcw.at<float>(2, 1)*x3Dw.at<float>(1) + Rcw.at<float>(2, 2)*x3Dw.at<float>(2) + tcw.at<float>(2));
		//this writing is faster because directly visit memory


		if (invzc < 0)
			return false;

		float u = Indices[0] * xc*invzc + Indices[2];
		float v = Indices[1] * yc*invzc + Indices[3];

		if (u<Indices[4] || u>Indices[6])
			return false;

		if (v<Indices[5] || v>Indices[7])
			return false;

		Projection.x = u;
		Projection.y = v;

		return true;
	}

	void ParticleFilter::OutputProjections(Frame& CurrentFrame, const Frame& LastFrame)
	{
		const cv::Mat Rcw = preResultOfFiltering.RT.rowRange(0, 3).colRange(0, 3);
		const cv::Mat tcw = preResultOfFiltering.RT.rowRange(0, 3).col(3);

		int N = LastFrame.N;

		vector<cv::Point2f> Projections;
		Projections.reserve(N);

		char InputDir[80];
		sprintf(InputDir, "D:\\sequences\\x1\\image_0\\%06d.png", nFrameIdx+15);
		cv::Mat mlm = cv::imread(InputDir, CV_LOAD_IMAGE_COLOR);
		//cv::Mat mlm(480, 640, CV_32FC3, cv::Scalar(255, 255, 255));

		int MapPointCount = 0;
		const int r = 5;
		for (int i_MP = 0; i_MP < N; i_MP++)
		{
			MapPoint* pMP = LastFrame.mvpMapPoints[i_MP];
			if (pMP)
			{
				if (!LastFrame.mvbOutlier[i_MP])// Current mappoint exist in last frame
				{
					cv::Mat x3Dw = pMP->GetWorldPos();
					//cout << x3Dw.t() << endl;
					cv::Mat x3Dc = cv::Mat(3, 1, CV_32F, 0.0);
					for (int i = 0; i < 3; i++)
					{
						x3Dc.at<float>(i) =
							Rcw.at<float>(i, 0)*x3Dw.at<float>(0) +
							Rcw.at<float>(i, 1)*x3Dw.at<float>(1) +
							Rcw.at<float>(i, 2)*x3Dw.at<float>(2) +
							tcw.at<float>(i);
					}

					const float xc = x3Dc.at<float>(0);
					const float yc = x3Dc.at<float>(1);
					const float invzc = 1.0 / x3Dc.at<float>(2);

					if (invzc < 0)
						continue;

					float u = CurrentFrame.fx*xc*invzc + CurrentFrame.cx;
					float v = CurrentFrame.fy*yc*invzc + CurrentFrame.cy;

					cv::Point2f NewProjection(u, v);


					cv::Point2f pt1, pt2;
					pt1.x = NewProjection.x - r;
					pt1.y = NewProjection.y - r;
					pt2.x = NewProjection.x + r;
					pt2.y = NewProjection.y + r;

					cv::rectangle(mlm, pt1, pt2, cv::Scalar(0, 255, 0));
					cv::circle(mlm, NewProjection, 2, cv::Scalar(0, 255, 0), -1);
					MapPointCount++;
				}
			}
		}
		imwrite("D:\\analysis\\ProjectionsImages\\Projections.png", mlm);
		ofstream OutputMP;
		OutputMP.open("D:\\analysis\\ProjectionsImages\\Projections.txt");
		OutputMP << MapPointCount << "\t" << nFrameIdx << endl;
	}

	void ParticleFilter::OutputProjections(Frame& CurrentFrame, vector<cv::Point3f>& vp3fNeMapPoints)
	{
		const cv::Mat Rcw = preResultOfFiltering.RT.rowRange(0, 3).colRange(0, 3);
		const cv::Mat tcw = preResultOfFiltering.RT.rowRange(0, 3).col(3);
		//cout << vp3fNeMapPoints.size()<<endl;

		char InputDir[80];
		sprintf(InputDir, "D:\\sequences\\x1\\image_0\\%06d.png", nFrameIdx+15);
		cv::Mat mlm = cv::imread(InputDir, CV_LOAD_IMAGE_COLOR);

		//cv::Mat mlm(480, 640, CV_32FC3, cv::Scalar(255, 255, 255));

		const int r = 5;
		for (int i = 0, i_end = vp3fNeMapPoints.size(); i < i_end; i++)
		{
			//cout << vp3fNeMapPoints[i].x << "\t" << vp3fNeMapPoints[i].y << "\t" << vp3fNeMapPoints[i].z << endl;
			cv::Mat x3Dc = cv::Mat(3, 1, CV_32F, 0.0);
			for (int j = 0; j < 3; j++)
			{
				x3Dc.at<float>(j) =
					Rcw.at<float>(j, 0)*vp3fNeMapPoints[i].x +
					Rcw.at<float>(j, 1)*vp3fNeMapPoints[i].y +
					Rcw.at<float>(j, 2)*vp3fNeMapPoints[i].z +
					tcw.at<float>(j);
			}
			//cout << x3Dc.t() << "\t";

			const float xc = x3Dc.at<float>(0);
			const float yc = x3Dc.at<float>(1);
			const float invzc = 1.0 / x3Dc.at<float>(2);

			if (invzc < 0)
				continue;

			float u = CurrentFrame.fx*xc*invzc + CurrentFrame.cx;
			float v = CurrentFrame.fy*yc*invzc + CurrentFrame.cy;

			//cout << u << "\t" << v << endl;

			cv::Point2f NewProjection(u, v);

			cv::Point2f pt1, pt2;
			pt1.x = NewProjection.x - r;
			pt1.y = NewProjection.y - r;
			pt2.x = NewProjection.x + r;
			pt2.y = NewProjection.y + r;

			cv::rectangle(mlm, pt1, pt2, cv::Scalar(0, 255, 0));
			cv::circle(mlm, NewProjection, 2, cv::Scalar(0, 255, 0), -1);
		}
		imwrite("D:\\analysis\\ProjectionsImages\\Projections2.png", mlm);
		ofstream OutputMP;
		OutputMP.open("D:\\analysis\\ProjectionsImages\\Projections2.txt");
		OutputMP << vp3fNeMapPoints.size() << "\t" << nFrameIdx << endl;
	}

	void ParticleFilter::OutputLastFrameMaxWeightPose()
	{
		ofstream OutputLastMaxWPoseDir;
		char OutputDir[80];
		sprintf(OutputDir, "D:\\analysis\\LastFrameMaxWeightParticlePose\\LastMaxWPose_%06d.txt", nFrameIdx);
		OutputLastMaxWPoseDir.open(OutputDir);
		for (int i = 0; i < ParticleNumber; i++)
		{
			for (int j = 0; j < 3; j++)
				for (int k = 0; k < 4; k++)
					OutputLastMaxWPoseDir << MaxWeightParticle.Pose.RT.at<float>(j, k) << "\t";
			//OutputParticleCloudDir << My_PF[i]->Weight << "\t";
			OutputLastMaxWPoseDir << endl;
		}

	}

	void ParticleFilter::OutputAllMPs(Frame* LastFrame)
	{
		ofstream OutputAllMapPointsDir;
		char OutputDir[80];
		sprintf(OutputDir, "D:\\analysis\\AllMapPoints\\AllMapPoints_%06d.txt", nFrameIdx);
		OutputAllMapPointsDir.open(OutputDir);
		int N = LastFrame->N;
		for (int i_MP = 0; i_MP < N; i_MP++)
		{
			MapPoint* pMP = LastFrame->mvpMapPoints[i_MP];
			if (pMP)
			{
				if (!LastFrame->mvbOutlier[i_MP])// Current mappoint exist in last frame
				{
					cv::Mat x3Dw = pMP->GetWorldPos();
					cv::Point2f Projection;
					OutputAllMapPointsDir << x3Dw.at<float>(0) << "\t" << x3Dw.at<float>(1) << "\t" << x3Dw.at<float>(2) << endl;
				}
			}
		}
	}

	int ParticleFilter::MatchingByProjection(Frame* CurrentFrame, const Frame* LastFrame, int EnlargeTh)
	{
		//unique_lock<mutex> lock1(MapPoint::mGlobalMutex);
		//FeaturesPosition.open("D:\\analysis\\FeaturePositions\\FeaturesPosition.txt");
		int N = CurrentFrame->N;
		int nmatches = 0;

		// Rotation Histogram (to check rotation consistency)
		vector<int> rotHist[HISTO_LENGTH];
		for (int i = 0; i<HISTO_LENGTH; i++)
			rotHist[i].reserve(500);
		const float factor = 1.0f / HISTO_LENGTH;

		vector<bool> vbMPMatched(LastFrame->N, false);
		vector<int> vnMPMatchedIdx(LastFrame->N, -1);
		vector<int> vnDistance(LastFrame->N, -1);

		KeysToRect.clear();
		KeysToRect.resize(CurrentFrame->N, -1);

		//cout << "MatchingByProjection: KeyPoints in CurrentFrame: " << N << endl;
		int CountTryMatching = 0;
		vector<int> AvgMatchingTimes;
		AvgMatchingTimes.reserve(CurrentFrame->N);

		// Projections in the last frame

		const cv::Mat RcwL = preResultOfFiltering.RT.rowRange(0, 3).colRange(0, 3);
		const cv::Mat tcwL = preResultOfFiltering.RT.rowRange(0, 3).col(3);

		int N_Last = LastFrame->N;
		vector<cv::Point2f> ProjectionInLastFrame;
		ProjectionInLastFrame.reserve(N_Last);
		ProjectionInLastFrame.resize(N_Last, cv::Point2f(0.0, 0.0));
		for (int i_Projection = 0; i_Projection < N_Last; i_Projection++)// assign necessary MapPoints to grids
		{
			if (vbNecessary[i_Projection]){
				MapPoint* pMP = LastFrame->mvpMapPoints[i_Projection];
				if (pMP){
					if (!LastFrame->mvbOutlier[i_Projection]){
						int PosX, PosY;
						cv::Mat x3Dw = pMP->GetWorldPos();
						cv::Point2f ProjectionTry;
						if(Project(RcwL, tcwL, x3Dw, Indices, ProjectionTry))
							ProjectionInLastFrame[i_Projection] = ProjectionTry;
						//cout << ProjectionTry.x << "\t" << ProjectionTry.y << endl;
					}
				}
			}
		}
		// Projections in the last frame

		for (int i = 0; i < CurrentFrame->N; i++)//match from features in current frame to mappoints projection ranges
		{
			const cv::Mat &d = CurrentFrame->mDescriptors.row(i);

			int bestDist = 256;
			int bestIdx2 = -1;

			int nCurrentOctave = CurrentFrame->mvKeys[i].octave;

			//FeaturesPosition << CurrentFrame.mvKeys[i].pt.x << "\t" << CurrentFrame.mvKeys[i].pt.y << endl;

			const bool bCheckLevels = (nCurrentOctave - 1 > 0) || (nCurrentOctave + 1 >= 0);

			int MatchingTime = 0;

			for (int i2 = 0; i2 < LastFrame->N; i2++)//try to match with each map point in last frame
			{
				/*if (abs(CurrentFrame->mvKeysUn[i].pt.x - ProjectionInLastFrame[i2].x)> SearchRange || abs(CurrentFrame->mvKeysUn[i].pt.y - ProjectionInLastFrame[i2].y) > SearchRange)
					continue;*/
				if (abs(CurrentFrame->mvKeysUn[i].pt.x * InvGridX - ProjectionInLastFrame[i2].x * InvGridX) > GridRange ||
					abs(CurrentFrame->mvKeysUn[i].pt.y * InvGridY - ProjectionInLastFrame[i2].y * InvGridY) > GridRange)
					continue;
				/*if (CurrentFrame->mvKeysUn[i].pt.x<vpairRectangle[i2].first.x - EnlargeTh || CurrentFrame->mvKeysUn[i].pt.x>vpairRectangle[i2].second.x + EnlargeTh)
					continue;
				if (CurrentFrame->mvKeysUn[i].pt.y<vpairRectangle[i2].first.y - EnlargeTh || CurrentFrame->mvKeysUn[i].pt.y>vpairRectangle[i2].second.y + EnlargeTh)
					continue;*/

				if (vbNecessary[i2])
				{
					MapPoint* pMP = LastFrame->mvpMapPoints[i2];
					if (pMP)
					{
						if (!LastFrame->mvbOutlier[i2])
						{
							CountTryMatching++;

							int nMapPointLastOctave = LastFrame->mvKeys[i2].octave;//get octave of this map point in lastframe
							if (bCheckLevels)
							{
								if (nMapPointLastOctave < nCurrentOctave - 1)
									continue;
								if (nCurrentOctave + 1 >= 0)
									if (nMapPointLastOctave>nCurrentOctave + 1)
										continue;
							}

							const cv::Mat dMP = pMP->GetDescriptor();
							const int dist = DescriptorDistance(dMP, d);

							MatchingTime++;

							if (dist < bestDist)
							{
								bestDist = dist;
								bestIdx2 = i2;
								KeysToRect[i] = i2;
							}
						}
					}
				}
			}
			AvgMatchingTimes.push_back(MatchingTime);
			if (bestDist <= TH_HIGH)
			{
				if (!vbMPMatched[bestIdx2])//this map point is matched to a feature point for the first time
				{
					vbMPMatched[bestIdx2] = true;
					vnMPMatchedIdx[bestIdx2] = i;
					vnDistance[bestIdx2] = bestDist;
					//CurrentFrame.mvbOutlier[bestIdx2] = false;

					CurrentFrame->mvpMapPoints[i] = LastFrame->mvpMapPoints[bestIdx2];
					nmatches++;
				}
				else//this map point has been matched to a feature point
				{
					if (bestDist < vnDistance[bestIdx2])// A new feature is matched to this map point better
					{
						CurrentFrame->mvpMapPoints[vnMPMatchedIdx[bestIdx2]] = static_cast<MapPoint*>(NULL);
						KeysToRect[vnMPMatchedIdx[bestIdx2]] = -1;
						CurrentFrame->mvpMapPoints[i] = LastFrame->mvpMapPoints[bestIdx2];
						vnMPMatchedIdx[bestIdx2] = i;
						vnDistance[bestIdx2] = bestDist;
					}
				}
			}
		}
		int TotalTime = 0;
		for (int i_Feature = 0; i_Feature < CurrentFrame->N; i_Feature++)
			TotalTime += AvgMatchingTimes[i_Feature];

		ProMatchNumber << TotalTime << "\t";
		
		for (int i = 0; i < LastFrame->N; i++)
		{

			if (vnMPMatchedIdx[i] != -1){
				float rot = LastFrame->mvKeysUn[i].angle - CurrentFrame->mvKeysUn[vnMPMatchedIdx[i]].angle;
				if (rot < 0.0)
					rot += 360.0f;
				int bin = round(rot*factor);
				if (bin == HISTO_LENGTH)
					bin = 0;
				assert(bin >= 0 && bin < HISTO_LENGTH);
				rotHist[bin].push_back(vnMPMatchedIdx[i]);
			}
			else continue;
		}
		
		//cout << "time in making histogram: " << tmakingHist << "\t";

		//cout << "before roration histogram, there are " << nmatches << "nmatches" << endl;
		//Apply rotation consistency
		int ind1 = -1;
		int ind2 = -1;
		int ind3 = -1;

		ComputeThreeMaxima(rotHist, HISTO_LENGTH, ind1, ind2, ind3);

		for (int i = 0; i<HISTO_LENGTH; i++)
		{
			if (i != ind1 && i != ind2 && i != ind3)
			{
				for (size_t j = 0, jend = rotHist[i].size(); j<jend; j++)
				{
					CurrentFrame->mvpMapPoints[rotHist[i][j]] = static_cast<MapPoint*>(NULL);
					nmatches--;
				}
			}
		}
		//cout << "Matching by projection number: " << nmatches << endl;
		return nmatches;
	}

	int ParticleFilter::MatchingByProjection2(Frame* CurrentFrame, const Frame* LastFrame, int EnlargeTh)
	{
		int N = CurrentFrame->N;
		int N_Last = LastFrame->N;
		int nmatches = 0;

		// Rotation Histogram (to check rotation consistency)
		vector<int> rotHist[HISTO_LENGTH];
		for (int i = 0; i<HISTO_LENGTH; i++)
			rotHist[i].reserve(500);
		const float factor = 1.0f / HISTO_LENGTH;

		vector<bool> vbMPMatched(LastFrame->N, false);
		vector<int> vnMPMatchedIdx(LastFrame->N, -1);
		vector<int> vnDistance(LastFrame->N, -1);

		KeysToRect.clear();
		KeysToRect.resize(CurrentFrame->N, -1);

		// grid based matching by projection
		const cv::Mat Rcw = preResultOfFiltering.RT.rowRange(0, 3).colRange(0, 3);
		const cv::Mat tcw = preResultOfFiltering.RT.rowRange(0, 3).col(3);

		for (int i_Projection = 0; i_Projection < N_Last; i_Projection++)// assign necessary MapPoints to grids
		{
			if (vbNecessary[i_Projection]){
				MapPoint* pMP = LastFrame->mvpMapPoints[i_Projection];
				if (pMP){
					if (!LastFrame->mvbOutlier[i_Projection]){
						int PosX, PosY;
						cv::Mat x3Dw = pMP->GetWorldPos();
						cv::Point2f ProjectionTry;
						if (Project(Rcw, tcw, x3Dw, Indices, ProjectionTry)){
							bool bInFrame = PosInGrid(ProjectionTry, PosX,PosY);
							if (bInFrame)
							{
								vmpProjectionGrid[PosX][PosY].push_back(pMP);
								vstMPIndices[PosX][PosY].push_back(i_Projection);
							}

							else
							{
								cout << "Projection: " << PosX << "\t" << PosY << " is not in frame" << endl;
							}

						}
					}
				}
			}
		}

		//cout << "MatchingByProjection: KeyPoints in CurrentFrame: " << N << endl;
		int CountTryMatching = 0;
		int PosX, PosY;
		int PosXtemp, PosYtemp;
		vmpLocalMapPoints.clear();
		vstLocalMpIndices.clear();

		vector<int> AvgMatchingTimes;
		AvgMatchingTimes.reserve(CurrentFrame->N);
		for (int i = 0; i < CurrentFrame->N; i++)//match from features in current frame to mappoints projection ranges
		{
			const cv::Mat &d = CurrentFrame->mDescriptors.row(i);

			int bestDist = 256;
			int bestIdx2 = -1;

			int nCurrentOctave = CurrentFrame->mvKeys[i].octave;

			const bool bCheckLevels = (nCurrentOctave - 1 > 0) || (nCurrentOctave + 1 >= 0);
			// Get map points whose projection in the last frame is in neighbour grids of this feature point
			PosInGrid(CurrentFrame->mvKeys[i].pt, PosX, PosY);

			for (int i_MP = 0; i_MP < 9; ++i_MP)
			{
				PosXtemp = PosX + Around9X[i_MP];
				PosYtemp = PosY + Around9Y[i_MP];
				if (PosXtemp<0 || PosXtemp>(PROJECTION_GRID_COL - 1) || PosYtemp<0 || PosYtemp>(PROJECTION_GRID_ROW - 1))
					continue;
				else
					GetMapPointsInGrid(PosXtemp, PosYtemp, vstLocalMpIndices, vmpLocalMapPoints);
			}
			//AvgMatchingTimes.push_back(vstLocalMpIndices.size());
			// Get map points whose projection int he last frame is in neighbour grids of this feature point
			//cout << "local map point size: " << vstLocalMpIndices.size() << endl;
			int MatchingCount = 0;
			for (int i2 = 0, i2_end = vstLocalMpIndices.size(); i2 < i2_end; ++i2)//try to match with each map point in last frame
			{
				MapPoint* pMP = vmpLocalMapPoints[i2];
				CountTryMatching++;
				int nMapPointLastOctave = LastFrame->mvKeys[vstLocalMpIndices[i2]].octave;//get octave of this map point in lastframe
				if (bCheckLevels)
				{
					if (nMapPointLastOctave < nCurrentOctave - 1)
						continue;
					if (nCurrentOctave + 1 >= 0)
						if (nMapPointLastOctave>nCurrentOctave + 1)
							continue;
				}

				const cv::Mat dMP = pMP->GetDescriptor();
				const int dist = DescriptorDistance(dMP, d);
				MatchingCount++;
				if (dist < bestDist)
				{
					bestDist = dist;
					bestIdx2 = vstLocalMpIndices[i2];
					KeysToRect[i] = vstLocalMpIndices[i2];
				}

			}
			AvgMatchingTimes.push_back(MatchingCount);
			if (bestDist <= TH_HIGH)
			{
				if (!vbMPMatched[bestIdx2])//this map point is matched to a feature point for the first time
				{
					vbMPMatched[bestIdx2] = true;
					vnMPMatchedIdx[bestIdx2] = i;
					vnDistance[bestIdx2] = bestDist;
					//CurrentFrame.mvbOutlier[bestIdx2] = false;

					CurrentFrame->mvpMapPoints[i] = LastFrame->mvpMapPoints[bestIdx2];
					nmatches++;
				}
				else//this map point has been matched to a feature point
				{
					if (bestDist < vnDistance[bestIdx2])// A new feature is matched to this map point better
					{
						CurrentFrame->mvpMapPoints[vnMPMatchedIdx[bestIdx2]] = static_cast<MapPoint*>(NULL);
						KeysToRect[vnMPMatchedIdx[bestIdx2]] = -1;
						CurrentFrame->mvpMapPoints[i] = LastFrame->mvpMapPoints[bestIdx2];
						vnMPMatchedIdx[bestIdx2] = i;
						vnDistance[bestIdx2] = bestDist;
					}
				}
			}
		}
		int totalMatchingTime = 0;
		for (int i_a = 0; i_a < AvgMatchingTimes.size(); i_a++)
			totalMatchingTime += AvgMatchingTimes[i_a];
		ProMatchNumber << totalMatchingTime << "\t";
		//cout << "Matching finish, start rotation check" << endl;
		for (int i = 0; i < LastFrame->N; i++)
		{

			if (vnMPMatchedIdx[i] != -1){
				float rot = LastFrame->mvKeysUn[i].angle - CurrentFrame->mvKeysUn[vnMPMatchedIdx[i]].angle;
				if (rot < 0.0)
					rot += 360.0f;
				int bin = round(rot*factor);
				if (bin == HISTO_LENGTH)
					bin = 0;
				assert(bin >= 0 && bin < HISTO_LENGTH);
				rotHist[bin].push_back(vnMPMatchedIdx[i]);
			}
			else continue;
		}

		//cout << "time in making histogram: " << tmakingHist << "\t";

		//cout << "before roration histogram, there are " << nmatches << "nmatches" << endl;
		//Apply rotation consistency
		int ind1 = -1;
		int ind2 = -1;
		int ind3 = -1;

		ComputeThreeMaxima(rotHist, HISTO_LENGTH, ind1, ind2, ind3);

		for (int i = 0; i<HISTO_LENGTH; i++)
		{
			if (i != ind1 && i != ind2 && i != ind3)
			{
				for (size_t j = 0, jend = rotHist[i].size(); j<jend; j++)
				{
					CurrentFrame->mvpMapPoints[rotHist[i][j]] = static_cast<MapPoint*>(NULL);
					nmatches--;
				}
			}
		}
		//cout << "Matching by projection number: " << nmatches << endl;
		return nmatches;
	}

	bool ParticleFilter::PosInGrid(cv::Point2f &Projection, int& PosX, int& PosY)
	{
		PosX = floor(Projection.x*InvGridX);
		PosY = floor(Projection.y*InvGridY);
		//cout << Projection.x << "\t" << Projection.y << endl;
		//Keypoint's coordinates are undistorted, which could cause to go out of the image
		if (PosX<0 || PosX >= PROJECTION_GRID_COL || PosY<0 || PosY >= PROJECTION_GRID_ROW)
			return false;

		return true;
	}

	void ParticleFilter::GetMapPointsInGrid(int& PosX, int& PosY, vector<size_t>& vstMapPointIndices, vector<MapPoint*>& vmpMapPoints)
	{
		vector<size_t> vstCell = vstMPIndices[PosX][PosY];
		vector<MapPoint*> vmpCell = vmpProjectionGrid[PosX][PosY];
		for (int i = 0, i_end = vstCell.size(); i < i_end; ++i)
		{
			vstMapPointIndices.push_back(vstCell[i]);
			vmpMapPoints.push_back(vmpCell[i]);
		}
	}

	int ParticleFilter::DescriptorDistance(const cv::Mat &a, const cv::Mat &b)
	{
		const int *pa = a.ptr<int32_t>();
		const int *pb = b.ptr<int32_t>();

		int dist = 0;

		for (int i = 0; i<8; i++, pa++, pb++)
		{
			unsigned  int v = *pa ^ *pb;
			v = v - ((v >> 1) & 0x55555555);
			v = (v & 0x33333333) + ((v >> 2) & 0x33333333);
			dist += (((v + (v >> 4)) & 0xF0F0F0F) * 0x1010101) >> 24;
		}

		return dist;
	}

	void ParticleFilter::ComputeThreeMaxima(vector<int>* histo, const int L, int &ind1, int &ind2, int &ind3)
	{
		int max1 = 0;
		int max2 = 0;
		int max3 = 0;

		for (int i = 0; i<L; i++)
		{
			const int s = histo[i].size();
			if (s>max1)
			{
				max3 = max2;
				max2 = max1;
				max1 = s;
				ind3 = ind2;
				ind2 = ind1;
				ind1 = i;
			}
			else if (s>max2)
			{
				max3 = max2;
				max2 = s;
				ind3 = ind2;
				ind2 = i;
			}
			else if (s>max3)
			{
				max3 = s;
				ind3 = i;
			}
		}

		if (max2<0.1f*(float)max1)
		{
			ind2 = -1;
			ind3 = -1;
		}
		else if (max3<0.1f*(float)max1)
		{
			ind3 = -1;
		}
	}

	void ParticleFilter::OutputRectangles(Frame& LastFrame)
	{
		ofstream OutputRectangle;
		char OutputDir[80];
		sprintf(OutputDir, "D:\\analysis\\PR-ORBSLAM rectangle test\\Rectangles_%06d.txt", nFrameIdx);
		OutputRectangle.open(OutputDir);
		//cout << "size of rectangles: " << vpairRectangle.size() << endl;
		for (int i_MP = 0; i_MP < LastFrame.N - 1; i_MP++)
		{
			MapPoint* pMP = LastFrame.mvpMapPoints[i_MP];
			if (pMP)
			{
				if (!LastFrame.mvbOutlier[i_MP] && vpairRectangle[i_MP].second.x!=0)// Current mappoint exist in last frame
				{
					OutputRectangle
						<< vpairRectangle[i_MP].second.x - vpairRectangle[i_MP].first.x << "\t"
						<< vpairRectangle[i_MP].second.y - vpairRectangle[i_MP].first.y << endl;
				}
			}
		}
	}

	void ParticleFilter::OutputPrediction()
	{
		vector<float> MaxIndices(6);
		vector<float> MinIndices(6);

		MaxIndices = My_PF[0]->Pose.Idx;
		MinIndices = My_PF[0]->Pose.Idx;

		for (int i_Particle = 1; i_Particle < ParticleNumber; i_Particle++)
		{	
			for (int i_Idx = 0; i_Idx < 6; i_Idx++)
			{
				if (My_PF[i_Particle]->Pose.Idx[i_Idx]>MaxIndices[i_Idx])
					MaxIndices[i_Idx] = My_PF[i_Particle]->Pose.Idx[i_Idx];
				else
					MinIndices[i_Idx] = My_PF[i_Particle]->Pose.Idx[i_Idx];
			}
		}		
		ofstream OutputPredictionMaxMin;
		char OutputDir[80];
		sprintf(OutputDir, "D:\\analysis\\PredictionResult\\PredictionResult_%06d.txt", nFrameIdx);
		OutputPredictionMaxMin.open(OutputDir);
		for (int i_Idx = 0; i_Idx < 6; i_Idx++)
		{
			OutputPredictionMaxMin << MaxIndices[i_Idx] << "\t";
		}
		OutputPredictionMaxMin << endl;
		for (int i_Idx = 0; i_Idx < 6; i_Idx++)
		{
			OutputPredictionMaxMin << MinIndices[i_Idx] << "\t";
		}
		OutputPredictionMaxMin << endl;
	}

	void ParticleFilter::OutputMatching(Frame& CurrentFrame)
	{
		ofstream OutputMatchings;
		char OutputDir[80];
		sprintf(OutputDir, "D:\\analysis\\matching\\Matchings_%06d.txt", nFrameIdx);
		OutputMatchings.open(OutputDir);
		for (int i_M = 0; i_M < KeysToRect.size(); i_M++)
			if (KeysToRect[i_M]!=-1)
			{
				OutputMatchings << vpairRectangle[KeysToRect[i_M]].first.x << "\t" << vpairRectangle[KeysToRect[i_M]].first.y << "\t" << vpairRectangle[KeysToRect[i_M]].second.x << "\t" << vpairRectangle[KeysToRect[i_M]].second.y << "\t";
				OutputMatchings << CurrentFrame.mvKeysUn[i_M].pt.x << "\t" << CurrentFrame.mvKeysUn[i_M].pt.y << endl;
			}

	}

	void ParticleFilter::OutputMPs(vector<cv::Point3f> MapPoints, vector<int> viObs)
	{
		ofstream OutputMapPoints;
		char OutputDir[80];
		sprintf(OutputDir, "D:\\analysis\\MapPoints\\MapPoints_%06d.txt", nFrameIdx+1);
		OutputMapPoints.open(OutputDir);
		for (int i_M = 0; i_M < MapPoints.size(); i_M++)
			if (MapPoints[i_M].x != 0)
				OutputMapPoints << MapPoints[i_M].x << "\t" << MapPoints[i_M].y << "\t" << MapPoints[i_M].z << "\t" << viObs[i_M] << endl;
	}

	void ParticleFilter::OutputDescriptors(Frame* CurrentFrame, Frame* LastFrame)
	{
		ofstream OutputKeys;
		ofstream OutputMPs;
		char OutputDirKeys[80];
		char OutputMapKeys[80];
		sprintf(OutputDirKeys, "D:\\analysis\\OutputKeyPointsDescriptors\\KeyPointsDiscriptors_%06d.txt", nFrameIdx + 1);
		sprintf(OutputMapKeys, "D:\\analysis\\OutputMapPointsDescriptors\\MapPointsDiscriptors_%06d.txt", nFrameIdx + 1);
		OutputKeys.open(OutputDirKeys);
		OutputMPs.open(OutputMapKeys);
		for (int i = 0; i < CurrentFrame->N; i++)
		{
			const cv::Mat &d = CurrentFrame->mDescriptors.row(i);
			const int *pa = d.ptr<int32_t>();

			OutputKeys << CurrentFrame->mvKeysUn[i].pt.x << "\t" << CurrentFrame->mvKeysUn[i].pt.y << "\t";

			for (int j = 0; j < 8; j++, pa++)
				OutputKeys << hex << *pa << "\t";
			OutputKeys << endl;
		}
		for (int i = 0; i < LastFrame->N; i++)
		{
			MapPoint* pMP = LastFrame->mvpMapPoints[i];
			if (pMP)
			{
				if (!LastFrame->mvbOutlier[i])
				{
					const cv::Mat dMP = pMP->GetDescriptor();
					const int *pa = dMP.ptr<int32_t>();

					for (int j = 0; j < 8; j++, pa++)
						OutputMPs << hex << *pa << "\t";
					OutputMPs << endl;
				}
			}
		}
	}

	void ParticleFilter::OutputLastFrameProjections(Frame* LastFrame)
	{
		ofstream OutputLastProjs;
		char OutputDirProjs[80];
		sprintf(OutputDirProjs, "D:\\analysis\\OutputLastFrameProjs\\LastFrameProjs_%06d.txt", nFrameIdx);
		OutputLastProjs.open(OutputDirProjs);

		int N = LastFrame->N;
		const cv::Mat Rcw = MaxWeightParticle.Pose.RT.rowRange(0, 3).colRange(0, 3);
		const cv::Mat tcw = MaxWeightParticle.Pose.RT.rowRange(0, 3).col(3);

		for (int i_MP = 0; i_MP < N; i_MP++)
		{
			MapPoint* pMP = LastFrame->mvpMapPoints[i_MP];
			if (pMP)
			{
				if (!LastFrame->mvbOutlier[i_MP])// Current mappoint exist in last frame
				{
					cv::Mat x3Dw = pMP->GetWorldPos();
					cv::Point2f Projection;
					if (Project(Rcw, tcw, x3Dw, Indices, Projection))
					{
						OutputLastProjs << Projection.x << "\t" << Projection.y << "\t" << pMP->Observations() << "\n";
					}
				}
			}
		}
	}

	void ParticleFilter::SelectiveEstimation()
	{
		CommonGrouping();
		CalculateRaidius();
		vector<RTMatrix> vrtWeightedGroupCenters(vvpCommonGroups.size());
		GroupLikelihoods.resize(vvpCommonGroups.size(), 0);
		OutputGroupL << nFrameIdx << "----------------------------------------------------" << endl;
		for (size_t i_Group = 0, i_Group_end = vvpCommonGroups.size(); i_Group < i_Group_end; i_Group++)
		{
			vrtWeightedGroupCenters[i_Group].Clone(CalculateWeigtedCenter(vvpCommonGroups[i_Group]));// calculate the weighted center
			GroupLikelihoods[i_Group] = CommonGroupingLikelihood(vvpCommonGroups[i_Group], vrtWeightedGroupCenters[i_Group], vAngleRadius[i_Group], vTransRadius[i_Group]);// calculate the grouplikelihood;
			if (GroupLikelihoods[i_Group] > MaxGroupL)
			{
				MaxGroupL = GroupLikelihoods[i_Group];
				MaxGroupLIdx = i_Group;
			}
		}
		//cout << "Frame          :" << nFrameIdx << endl;
		//cout << "Max GroupL Size: " << vvpCommonGroups[MaxGroupLIdx].size() << endl;
		//cout << "Group size     : " << vvpCommonGroups.size() << endl;

		GroupResult << nFrameIdx << "\t" << vvpCommonGroups[MaxGroupLIdx].size() << "\t" << vvpCommonGroups.size() << endl;
		/*cout << "Group size       : ";
		for (int i_Group = 0; i_Group < vvpCommonGroups.size(); i_Group++)
			cout << vvpCommonGroups[i_Group].size() << "\t";
		cout << endl;

		cout << "Group likelihood : ";
		for (int i_Group = 0; i_Group < vvpCommonGroups.size(); i_Group++)
			cout << GroupLikelihoods[i_Group] << "\t";
		cout << endl;

		cout << "Group AngleRadius: ";
		for (int i_Group = 0; i_Group < vvpCommonGroups.size(); i_Group++)
			cout << vAngleRadius[i_Group] << "\t";
		cout << endl;

		cout << "Group TransRadius: ";
		for (int i_Group = 0; i_Group < vvpCommonGroups.size(); i_Group++)
			cout << vTransRadius[i_Group] << "\t";
		cout << endl;*/

		fill(ResaultOfFiltering.Idx.begin(), ResaultOfFiltering.Idx.end(), 0.0);
		float SumOfWeight = 0;
		for (size_t i = 0, i_end = vvpCommonGroups[MaxGroupLIdx].size(); i < i_end; i++){
			for (size_t j = 0, j_end = 6; j < j_end; j++){
				if (!std::isfinite(My_PF[i]->Pose.Idx[j]))
					cout << "idx is not finite in particle: " << i << endl;
				ResaultOfFiltering.Idx[j] += vvpCommonGroups[MaxGroupLIdx][i]->Pose.Idx[j] * vvpCommonGroups[MaxGroupLIdx][i]->Weight;
				
			}
			SumOfWeight += vvpCommonGroups[MaxGroupLIdx][i]->Weight;
		}

		for (int i_Idx = 0; i_Idx < 6; i_Idx++)
			ResaultOfFiltering.Idx[i_Idx] = ResaultOfFiltering.Idx[i_Idx] / SumOfWeight;
		ResaultOfFiltering.Idx_to_Mat();
		OutputAllParticles(vvpCommonGroups);
	}

	void ParticleFilter::CommonGrouping()
	{
		vector<Particle*> vpNewGroup(1);
		vpNewGroup[0] = My_PF[0];
		vvpCommonGroups.push_back(vpNewGroup);
		vrtGroupCenter.push_back(My_PF[0]->Pose);// set the first particle as the first group

		unsigned int GroupNumber = 1;
		vector<float> vfDistance;
		vector<unsigned int> vuiOKGroupIdx;
		for (unsigned int i = 1; i<ParticleNumber; i++)
		{
			vfDistance.clear();
			vuiOKGroupIdx.clear();
			for (unsigned int Group_n = 0; Group_n < GroupNumber; Group_n++)//for each particle, check the distance between it and all of the group centers
			{
				vfDistance = RTMatrixDistance(My_PF[i]->Pose, vrtGroupCenter[Group_n]);
				if (vfDistance[0] < CommonGroupingAngleTH&&
					vfDistance[1] < CommonGroupingAngleTH&&
					vfDistance[2] < CommonGroupingAngleTH&&
					vfDistance[3] < CommonGroupingTransTH&&
					vfDistance[4] < CommonGroupingTransTH&&
					vfDistance[5] < CommonGroupingTransTH
					)
					vuiOKGroupIdx.push_back(Group_n);
			}
			if (vuiOKGroupIdx.size() == 0)//this particle belongs to no group, create a new group
			{
				vector<Particle*> vpNewGroup(1);
				vpNewGroup[0] = My_PF[i];
				vvpCommonGroups.push_back(vpNewGroup);
				vrtGroupCenter.push_back(My_PF[i]->Pose);
				GroupNumber++;
			}
			else if (vuiOKGroupIdx.size() == 1)//this particle belongs to only one group
			{
				RTMatrix NewGroupCenter;
				vvpCommonGroups[vuiOKGroupIdx[0]].push_back(My_PF[i]);
				NewGroupCenter = UpdateGroupCenter(vrtGroupCenter[vuiOKGroupIdx[0]], vuiOKGroupIdx[0], My_PF[i]->Pose);
			}
			else//this particle belongs to more than one groups, combine these groups
			{
				vector<Particle*> vpCombinedGroup;
				vpCombinedGroup.push_back(My_PF[i]);
				RTMatrix vfCombinedGroupCenter;
				for (size_t Group_Idx = 0, Group_Idx_end = vuiOKGroupIdx.size(); Group_Idx < Group_Idx_end; Group_Idx++)
				{
					for (size_t Particle_Idx = 0, Particle_Idx_end = vvpCommonGroups[vuiOKGroupIdx[Group_Idx]].size(); Particle_Idx < Particle_Idx_end; Particle_Idx++)
					{
						vpCombinedGroup.push_back(vvpCommonGroups[vuiOKGroupIdx[Group_Idx]][Particle_Idx]);//move particles into new group
						for (unsigned int i_Idx = 0; i_Idx < 6; i_Idx++)
							vfCombinedGroupCenter.Idx[i_Idx] += vvpCommonGroups[vuiOKGroupIdx[Group_Idx]][Particle_Idx]->Pose.Idx[i_Idx];
					}
				}
				for (int i_DeleteGroup = vuiOKGroupIdx.size() - 1, i_DeleteGroup_end = 0; i_DeleteGroup >= i_DeleteGroup_end; i_DeleteGroup--)
				{
					vvpCommonGroups.erase(vvpCommonGroups.begin() + vuiOKGroupIdx[i_DeleteGroup]);
					vrtGroupCenter.erase(vrtGroupCenter.begin() + vuiOKGroupIdx[i_DeleteGroup]);
				}
				int CombinedGroupSize = vpCombinedGroup.size();
				for (unsigned int i_Idx = 0; i_Idx < 6; i_Idx++)
					vfCombinedGroupCenter.Idx[i_Idx] = vfCombinedGroupCenter.Idx[i_Idx] / CombinedGroupSize;
				vfCombinedGroupCenter.Idx_to_Mat();
				vvpCommonGroups.push_back(vpCombinedGroup);
				vrtGroupCenter.push_back(vfCombinedGroupCenter);
				GroupNumber = GroupNumber - vuiOKGroupIdx.size() + 1;
			}
		}
	}

	vector<float> ParticleFilter::RTMatrixDistance(RTMatrix &Matrix1, RTMatrix &Matrix2)
	{
		vector<float> vfDistance(6);
		for (unsigned int i = 0; i < 6; i++)
			vfDistance[i] = abs(Matrix1.Idx[i] - Matrix2.Idx[i]);
		return vfDistance;
	}

	RTMatrix ParticleFilter::UpdateGroupCenter(RTMatrix &GroupCenter, unsigned int &GroupIdx, RTMatrix &NewParticlePose)
	{
		unsigned int GroupSize = vvpCommonGroups[GroupIdx].size();
		RTMatrix OldGroupCenter = vrtGroupCenter[GroupIdx];
		RTMatrix NewGroupCenter;
		for (unsigned int i_Idx = 0; i_Idx < 6; i_Idx++)
			NewGroupCenter.Idx[i_Idx] = (OldGroupCenter.Idx[i_Idx] * GroupSize + NewParticlePose.Idx[i_Idx]) / (GroupSize + 1);
		NewGroupCenter.Idx_to_Mat();

		return NewGroupCenter;
	}

	RTMatrix ParticleFilter::CalculateWeigtedCenter(vector<Particle*> Group)
	{
		size_t GroupSize = Group.size();
		float TotalWeight = 0;
		RTMatrix WeightedCenter;
		for (size_t i_Group = 0; i_Group < GroupSize; i_Group++)
		{
			for (int i_Idx = 0; i_Idx < 6; i_Idx++)
				WeightedCenter.Idx[i_Idx] += Group[i_Group]->Pose.Idx[i_Idx] * Group[i_Group]->Weight;
			TotalWeight += Group[i_Group]->Weight;
		}
		for (int i_Idx = 0; i_Idx < 6; i_Idx++)
			WeightedCenter.Idx[i_Idx] = WeightedCenter.Idx[i_Idx] / TotalWeight;
		WeightedCenter.Idx_to_Mat();
		return WeightedCenter;
	}

	void ParticleFilter::CalculateRaidius()
	{
		float MaxAngleRadius = 0;
		float MaxTransRadius = 0;

		vector<float> Distance(6, 0);
		float AngleRadius, TransRadius;
		for (size_t i_Group = 0, i_Group_end = vvpCommonGroups.size(); i_Group < i_Group_end; i_Group++)
		{
			for (size_t i_Particle = 0, i_Particle_end = vvpCommonGroups[i_Group].size(); i_Particle < i_Particle_end; i_Particle++)
			{
				Distance = RTMatrixDistance(vrtGroupCenter[i_Group], vvpCommonGroups[i_Group][i_Particle]->Pose);
				AngleRadius = Distance[0] * Distance[0] + Distance[1] * Distance[1] + Distance[2] * Distance[2];
				TransRadius = Distance[3] * Distance[3] + Distance[4] * Distance[4] + Distance[5] * Distance[5];
				if (AngleRadius > MaxAngleRadius)
					MaxAngleRadius = AngleRadius;
				if (TransRadius > MaxTransRadius)
					MaxTransRadius = TransRadius;
			}
			MaxAngleRadius = sqrt(MaxAngleRadius);
			MaxTransRadius = sqrt(MaxTransRadius);
			vAngleRadius.push_back(exp(-MaxAngleRadius));
			vTransRadius.push_back(exp(-MaxTransRadius));
		}
	}

	float ParticleFilter::CommonGroupingLikelihood(vector<Particle*> Group, RTMatrix WeightedCenter, float AngleRadius, float TransRadius)
	{
		vector<float> Distance;
		float KAngle = 0;
		float KTrans = 0;
		float GroupLikelihood = 0;
		if (Group.size() <= ZeroGroupLTH)
			return 0;
		else
		{
			for (size_t i_Particle = 0, i_Particle_end = Group.size(); i_Particle < i_Particle_end; i_Particle++)
			{
				Distance = RTMatrixDistance(Group[i_Particle]->Pose, WeightedCenter);
				float AngleDistance = sqrt(Distance[0] * Distance[0] + Distance[1] * Distance[1] + Distance[2] * Distance[2]);
				float TransDistance = sqrt(Distance[3] * Distance[3] + Distance[4] * Distance[4] + Distance[5] * Distance[5]);

				if (AngleDistance < 1e-6)
					KAngle = 1.0;
				else
					KAngle = (exp(-AngleDistance) - AngleRadius) / (1 - AngleRadius);
				
				if (TransDistance < 1e-6)
					KTrans = 1.0;
				else
					KTrans = (exp(-TransDistance) - TransRadius) / (1 - TransRadius);
				
				GroupLikelihood += KAngle*KTrans*Group[i_Particle]->Weight;
				OutputGroupL << GroupLikelihood << "\t" << AngleDistance << "\t" << TransDistance << "\t" << KAngle << "\t" << KTrans << "\t" << Group[i_Particle]->Likelihood << "\t" << Group[i_Particle]->Weight << endl;
			}

		}
		for (int i_Particle = 0; i_Particle < ParticleNumber; i_Particle++)
			OutputGroupL << My_PF[i_Particle]->Likelihood << "\t";
		return GroupLikelihood;
	}

	void ParticleFilter::OutputAllParticles(vector<vector<Particle*>> &vvpOutput)
	{
		ofstream GroupedParticles;
		char OutputDir[80];
		sprintf(OutputDir, "D:\\analysis\\ClusteringResult\\ClusteringResult_%06d.txt", nFrameIdx);
		GroupedParticles.open(OutputDir);
		for (size_t i_Particles = 0, i_Particles_end = vvpOutput.size(); i_Particles < i_Particles_end; i_Particles++)
		{
			for (size_t j_Particles = 0, j_Particles_end = vvpOutput[i_Particles].size(); j_Particles < j_Particles_end; j_Particles++)
			{
				for (int i = 0; i < 6; i++)
				{
					GroupedParticles << vvpOutput[i_Particles][j_Particles]->Pose.Idx[i] << "\t";
				}
				GroupedParticles << vvpOutput[i_Particles][j_Particles]->Weight << endl;
			}
			GroupedParticles << endl;
		}
	}
}//namespace ORB_SLAM2
