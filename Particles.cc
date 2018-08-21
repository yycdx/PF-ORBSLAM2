#include "Particle.h"

using namespace std;

namespace ORB_SLAM2
{
	const double GaussianInitial = 0.02;
	const double GaussianEulerAngle = 400.0;
	const double GaussianXYZ = 800.0;
	const double fps = 1000.0;
	const double invfps = 0.001;
	const double sqinvfps = 0.000001;
	const double halfsqinvfps = 0.0000005;
	const pair<unsigned int, unsigned int> Resolution = pair<unsigned int, unsigned int>(640, 480);
	const double DiffXTreshold = 15.0;
	const double DiffYTreshold = 15.0;
	const bool randomSeed = false;
	const int resolution_x = 640;
	const int resolution_y = 480;
	const int FixedPointBitNumber = 1000;

	/*const float PiecewiseInterval[21] = { 0.001,0.101,	0.201,	0.301,	0.401,	0.501,	0.601,	0.701,	0.801,	0.901,	1.001,	1.109,	1.248,	1.426,	1.663,	1.996,	2.493,	3.322,	4.975,	9.900, 1000};
	const float Piecewisek[20] = {-9900, -49.259, - 16.529, - 8.285, - 4.978, - 3.321, - 2.374, - 1.781, - 1.386, - 1.109, - 0.900, - 0.722, - 0.562, - 0.421, - 0.301, - 0.201, - 0.121, - 0.061, - 0.020, 0.000 };
	const float Piecewiseb[20] = {1009, 14.876, 8.297, 	5.816, 	4.490, 	3.660, 	3.090, 	2.675, 	2.358, 	2.109, 	1.900, 	1.702, 	1.502, 	1.302, 	1.102, 	0.902 ,	0.702 ,	0.502, 0.302, 0.102 };*/
	
	const float PiecewiseInterval[25] = { 0.001,	0.1,	0.2,	0.3,	0.4,	0.5,	0.6,	0.7,	0.8,	0.85,	0.9,	0.95,	1,	1.05,	1.1,	1.15,	1.2,	1.3,	1.4,	1.6,	2,	2.5,
		3.3,	4.9,	10};
	const float Piecewisek[24] = {
		-10000, -50, -16.66666667, -8.333333333, -5, -3.333333333, -2.380952381, -1.785714286, -1.470588235, -1.307189542, -1.169590643, -1.052631579, -0.952380952, -0.865800866, -0.790513834,
		-0.724637681, -0.641025641, -0.549450549, -0.446428571, -0.3125, -0.2, -0.121212121, -0.061842919, -0.020408163 };
	const float Piecewiseb[24] = {
		1010,	15,	8.333333333,	5.833333333,	4.5,	3.666666667,	3.095238095,	2.678571429,	2.426470588,	2.287581699,	2.16374269,	2.052631579,	1.952380952,	1.861471861,
		1.778656126,	1.702898551,	1.602564103, 1.483516484,	1.339285714, 1.125,	0.9,	0.703030303,	0.507111936,	0.304081633 };


	/*const float PiecewiseInterval[17] = { 0.001, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 };
	const float Piecewisek[16] = { -1000,-0.5,-0.1667,-0.08333,-0.05,-0.03333,-0.02381,-0.01786,-0.01389,-0.01111,-0.00909,-0.00758,-0.00641,-0.00579,-0.00476,-0.00417 };
	const float Piecewiseb[16] = {1001,1.5,0.8333,0.5833,0.45,0.3667,0.3095,0.2678,0.2361,0.2111,0.1909,0.1742,0.1602,0.1483,0.1381,0.1292};*/


	Particle::Particle() :ParticleAcceleration(6, 0.0), ParticleSpeed(6, 0.0)
	{
		Likelihood = 0.0;
		Weight = 0.0;
		WeightSum = 0.0;
		SumDiff = 0.0;
		SumDiffX = 0.0;
		SumDiffY = 0.0;
		OutlierCounter = 0;
		InlierCounter = 0;
	}

	void Particle::ParticleInitialization(RTMatrix &InitializationPose, vector<float> &InitialSpeed)
	{
		DUtils::Random::SeedRandOnce(0);
		for (size_t i = 0; i < 6; i++){
			Pose.Idx[i] = DUtils::Random::RandomGaussianValue(0.0, GaussianInitial)*InitializationPose.Idx[i] + InitializationPose.Idx[i];
			//while (!std::isfinite(Pose.Idx[i]))
			//	Pose.Idx[i] = DUtils::Random::RandomGaussianValue(0.0, GaussianInitial)*InitializationPose.Idx[i] + InitializationPose.Idx[i];
			//cout << Pose.Idx[i] << "\t";
		}
		//cout << endl;
		Pose.Idx_to_Mat();
		Rcw = Pose.RT.rowRange(0, 3).colRange(0, 3);
		tcw = Pose.RT.rowRange(0, 3).col(3);
		vp2fProjections.reserve(1600);
		ParticleSpeed = InitialSpeed;
		for (int i = 0; i < 6;i++)
			Speedvt.Idx[i] = InitialSpeed[i]*invfps;
		Speedvt.Idx_to_Mat();
		fill(ParticleAcceleration.begin(), ParticleAcceleration.end(), 0.0);
	}

	void Particle::Prediction()
	{
		fill(ParticleAcceleration.begin(), ParticleAcceleration.end(), 0.0);
		DUtils::Random::SeedRandOnce(0);
		for (size_t i = 0; i < 3; i++)//predict euler angles accelerations
			ParticleAcceleration[i] = DUtils::Random::RandomGaussianValue(0.0, GaussianEulerAngle);

		for (size_t i = 3; i < 6; i++)//predict xyz accelerations
			ParticleAcceleration[i] = DUtils::Random::RandomGaussianValue(0.0, GaussianXYZ);

		for (size_t i = 0; i < 6; i++){
			ParticleSpeed[i] += ParticleAcceleration[i] / fps;
			Pose.Idx[i] += ParticleSpeed[i] / fps + ParticleAcceleration[i] / (fps * fps * 2);//apply acceleration and spped to idx
		}
		//cout << Pose.Idx[5] << "\t";
		this->Pose.Idx_to_Mat();

		Rcw = Pose.RT.rowRange(0, 3).colRange(0, 3);
		tcw = Pose.RT.rowRange(0, 3).col(3);

		for (int i = 0; i < 6; i++)
			if (!isfinite(Pose.Idx[i]))
				cout << "Idx " << i << " is not finite in prediction" << endl;
	}

	void Particle::PredictionMat()
	{
		fill(ParticleAcceleration.begin(), ParticleAcceleration.end(), 0.0);
		DUtils::Random::SeedRandOnce(0);
		for (size_t i = 0; i < 3; i++)//predict euler angles accelerations
			ParticleAcceleration[i] = DUtils::Random::RandomGaussianValue(0.0, GaussianEulerAngle);

		for (size_t i = 3; i < 6; i++)//predict xyz accelerations
			ParticleAcceleration[i] = DUtils::Random::RandomGaussianValue(0.0, GaussianXYZ);

		for (size_t i = 3; i < 6; i++){
			ParticleSpeed[i] += ParticleAcceleration[i] / fps;
			Pose.Idx[i] += ParticleSpeed[i] / fps + ParticleAcceleration[i] / (fps * fps * 2);//apply acceleration and spped to idx
		}

		cv::Mat R_at2 = (cv::Mat_<float>(3, 3) <<
			cos(ParticleAcceleration[0] * sqinvfps) * cos(ParticleAcceleration[1] * sqinvfps),
			cos(ParticleAcceleration[0] * sqinvfps) * sin(ParticleAcceleration[1] * sqinvfps) * sin(ParticleAcceleration[2] * sqinvfps) - sin(ParticleAcceleration[0] * sqinvfps) * cos(ParticleAcceleration[2] * sqinvfps),
			cos(ParticleAcceleration[0] * sqinvfps) * sin(ParticleAcceleration[1] * sqinvfps) * cos(ParticleAcceleration[2] * sqinvfps) + sin(ParticleAcceleration[0] * sqinvfps) * sin(ParticleAcceleration[2] * sqinvfps),
			sin(ParticleAcceleration[0] * sqinvfps) * cos(ParticleAcceleration[1] * sqinvfps),
			sin(ParticleAcceleration[0] * sqinvfps) * sin(ParticleAcceleration[1] * sqinvfps) * sin(ParticleAcceleration[2] * sqinvfps) + cos(ParticleAcceleration[0] * sqinvfps) * cos(ParticleAcceleration[2] * sqinvfps),
			sin(ParticleAcceleration[0] * sqinvfps) * sin(ParticleAcceleration[1] * sqinvfps) * cos(ParticleAcceleration[2] * sqinvfps) - cos(ParticleAcceleration[0] * sqinvfps) * sin(ParticleAcceleration[2] * sqinvfps),
			-sin(ParticleAcceleration[1] * sqinvfps),
			cos(ParticleAcceleration[1] * sqinvfps) * sin(ParticleAcceleration[2] * sqinvfps),
			cos(ParticleAcceleration[1] * sqinvfps) * cos(ParticleAcceleration[2] * sqinvfps)
			);

		//Speedvt.RT.colRange(0, 3).rowRange(0, 3) = Speedvt.RT.colRange(0, 3).rowRange(0, 3) * R_at2;
		cv::Mat preSpeedvt = Speedvt.RT.clone();
		for (int i = 0; i < 3; i++)// this writing is faster than using operator * 
			for (int j = 0; j < 3; j++)
				Speedvt.RT.at<float>(i, j) =
				preSpeedvt.at<float>(i, 0)*R_at2.at<float>(0, j) +
				preSpeedvt.at<float>(i, 1)*R_at2.at<float>(1, j) +
				preSpeedvt.at<float>(i, 2)*R_at2.at<float>(2, j);

		Speedvt.RT.at<float>(0, 3) += ParticleAcceleration[3] * sqinvfps;
		Speedvt.RT.at<float>(1, 3) += ParticleAcceleration[4] * sqinvfps;
		Speedvt.RT.at<float>(2, 3) += ParticleAcceleration[5] * sqinvfps;
		Speedvt.Mat_to_Idx();

		cv::Mat R_hat2 = (cv::Mat_<float>(3, 3) <<
			cos(ParticleAcceleration[0] * halfsqinvfps) * cos(ParticleAcceleration[1] * halfsqinvfps),
			cos(ParticleAcceleration[0] * halfsqinvfps) * sin(ParticleAcceleration[1] * halfsqinvfps) * sin(ParticleAcceleration[2] * halfsqinvfps) - sin(ParticleAcceleration[0] * halfsqinvfps) * cos(ParticleAcceleration[2] * halfsqinvfps),
			cos(ParticleAcceleration[0] * halfsqinvfps) * sin(ParticleAcceleration[1] * halfsqinvfps) * cos(ParticleAcceleration[2] * halfsqinvfps) + sin(ParticleAcceleration[0] * halfsqinvfps) * sin(ParticleAcceleration[2] * halfsqinvfps),
			sin(ParticleAcceleration[0] * halfsqinvfps) * cos(ParticleAcceleration[1] * halfsqinvfps),
			sin(ParticleAcceleration[0] * halfsqinvfps) * sin(ParticleAcceleration[1] * halfsqinvfps) * sin(ParticleAcceleration[2] * halfsqinvfps) + cos(ParticleAcceleration[0] * halfsqinvfps) * cos(ParticleAcceleration[2] * halfsqinvfps),
			sin(ParticleAcceleration[0] * halfsqinvfps) * sin(ParticleAcceleration[1] * halfsqinvfps) * cos(ParticleAcceleration[2] * halfsqinvfps) - cos(ParticleAcceleration[0] * halfsqinvfps) * sin(ParticleAcceleration[2] * halfsqinvfps),
			-sin(ParticleAcceleration[1] * halfsqinvfps),
			cos(ParticleAcceleration[1] * halfsqinvfps) * sin(ParticleAcceleration[2] * halfsqinvfps),
			cos(ParticleAcceleration[1] * halfsqinvfps) * cos(ParticleAcceleration[2] * halfsqinvfps)
			);

		//cv::Mat R_vt = Speedvt.RT.colRange(0, 3).rowRange(0, 3)*R_at2;
		cv::Mat R_vt = cv::Mat(3, 3, CV_32F, 0.0);
		for (int i = 0; i < 3; i++)// this writing is faster than using operator * 
			for (int j = 0; j < 3; j++)
				R_vt.at<float>(i, j) =
				Speedvt.RT.at<float>(i, 0)*R_at2.at<float>(0, j) +
				Speedvt.RT.at<float>(i, 1)*R_at2.at<float>(1, j) +
				Speedvt.RT.at<float>(i, 2)*R_at2.at<float>(2, j);

		/*Pose.RT.colRange(0, 3).rowRange(0, 3) =
			Pose.RT.colRange(0, 3).rowRange(0, 3) *
			R_vt * 
			R_hat2;*/

		cv::Mat Rotation1 = cv::Mat(3, 3, CV_32F, 0.0);

		for (int i = 0; i < 3; i++)// this writing is faster than using operator * 
			for (int j = 0; j < 3; j++)
				Rotation1.at<float>(i, j) =
				Pose.RT.at<float>(i, 0)*R_vt.at<float>(0, j) +
				Pose.RT.at<float>(i, 1)*R_vt.at<float>(1, j) +
				Pose.RT.at<float>(i, 2)*R_vt.at<float>(2, j);

		for (int i = 0; i < 3; i++)// this writing is faster than using operator * 
			for (int j = 0; j < 3; j++)
				Pose.RT.at<float>(i, j) =
				Rotation1.at<float>(i, 0)*R_hat2.at<float>(0, j) +
				Rotation1.at<float>(i, 1)*R_hat2.at<float>(1, j) +
				Rotation1.at<float>(i, 2)*R_hat2.at<float>(2, j);

		Pose.RT.at<float>(0, 3) = Pose.Idx[3];
		Pose.RT.at<float>(1, 3) = Pose.Idx[4];
		Pose.RT.at<float>(2, 3) = Pose.Idx[5];

		Pose.Mat_to_Idx();

		Rcw = Pose.RT.rowRange(0, 3).colRange(0, 3);
		tcw = Pose.RT.rowRange(0, 3).col(3);

		for (int i = 0; i < 6; i++)
			if (!isfinite(Pose.Idx[i]))
				cout << "Idx " << i << " is not finite in prediction" << endl;
	}

	void Particle::PredictionTrans()
	{
		fill(ParticleAcceleration.begin(), ParticleAcceleration.end(), 0.0);
		DUtils::Random::SeedRandOnce(0);
		for (size_t i = 0; i < 3; i++)//predict euler angles accelerations
			ParticleAcceleration[i] = DUtils::Random::RandomGaussianValue(0.0, GaussianEulerAngle);

		for (size_t i = 3; i < 6; i++)//predict xyz accelerations
			ParticleAcceleration[i] = DUtils::Random::RandomGaussianValue(0.0, GaussianXYZ);

		for (size_t i = 0; i < 6; i++){
			ParticleSpeed[i] += ParticleAcceleration[i] / fps;
			Pose.Idx[i] += ParticleSpeed[i] / fps + ParticleAcceleration[i] / (fps * fps * 2);//apply acceleration and spped to idx
		}
		//cout << Pose.Idx[5] << "\t";
		for (int i = 0; i < 6; i++)
			if (!isfinite(Pose.Idx[i]))
				cout << "Idx " << i << " is not finite" << endl;
		this->Pose.Idx_to_Mat();
	}

	void Particle::PredictionRot()
	{
		fill(ParticleAcceleration.begin(), ParticleAcceleration.end(), 0.0);
		DUtils::Random::SeedRandOnce(0);
		for (size_t i = 0; i < 3; i++)//predict euler angles accelerations
			ParticleAcceleration[i] = DUtils::Random::RandomGaussianValue(0.0, GaussianXYZ);

		for (size_t i = 3; i < 6; i++)//predict xyz accelerations
			ParticleAcceleration[i] = DUtils::Random::RandomGaussianValue(0.0, GaussianEulerAngle);

		for (size_t i = 0; i < 6; i++){
			ParticleSpeed[i] += ParticleAcceleration[i] / fps;
			Pose.Idx[i] += ParticleSpeed[i] / fps + ParticleAcceleration[i] / (fps * fps * 2);//apply acceleration and spped to idx
		}
		//cout << Pose.Idx[5] << "\t";
		this->Pose.Idx_to_Mat();
	}

	void Particle::Observation(Frame& CurrentFrame)
	{
		int N = CurrentFrame.N;
		vector<bool> pOutlier1(CurrentFrame.N, false);
		this->pOutlier.clear();
		this->pOutlier = pOutlier1;
		for (size_t i = 0; i < N; i++){
			MapPoint* pMP = CurrentFrame.mvpMapPoints[i];
			if (pMP){
				//step 1
				cv::Point2f p2fProjection;

				const cv::Mat Rcw = this->Pose.RT.rowRange(0, 3).colRange(0, 3);
				const cv::Mat tcw = this->Pose.RT.rowRange(0, 3).col(3);

				cv::Mat x3Dw = pMP->GetWorldPos();
				cv::Mat x3Dc = Rcw*x3Dw + tcw;
				//step 2
				const float xc = x3Dc.at<float>(0);
				const float yc = x3Dc.at<float>(1);
				const float invzc = 1.0 / x3Dc.at<float>(2);

				if (invzc < 0)
				{
					pOutlier[i] = true;
					OutlierCounter++;
					continue;
				}


				float u = CurrentFrame.fx*xc*invzc + CurrentFrame.cx;
				float v = CurrentFrame.fy*yc*invzc + CurrentFrame.cy;

				p2fProjection.x = u;//Project from local 3D coordinate to 2D coordinate
				p2fProjection.y = v;
				//step 3
				float DiffX = abs(p2fProjection.x - CurrentFrame.mvKeys[i].pt.x);
				float DiffY = abs(p2fProjection.y - CurrentFrame.mvKeys[i].pt.y);
				//cout << DiffX << "\t" << DiffY << endl;
				if (DiffX > DiffXTreshold || DiffY > DiffYTreshold)
				{
					pOutlier[i] = true;
					OutlierCounter++;
					continue;
				}
				else
				{
					InlierCounter++;
					//SumDiff += DiffX*DiffX + DiffY*DiffY;
					SumDiffX += DiffX*DiffX;
					SumDiffY += DiffY*DiffY;
				}

			}
		}
		float AverageSumDiffX = SumDiffX / float(InlierCounter);
		float AverageSumDiffY = SumDiffY / float(InlierCounter);

		Likelihood = InlierCounter * (1 / AverageSumDiffX)*(1 / AverageSumDiffY);
		//Likelihood = 1 / (SumDiff + 1);
		//Likelihood = InlierCounter*(1 / SumDiffX)*(1 / SumDiffX);
		Likelihood = Likelihood*Likelihood*Likelihood;


		//Likelihood = (1 / SumDiffX)*(1 / SumDiffY);
		//cout << SumDiff << "\t" << Pose.Idx[5] << endl;
	}

	void Particle::ObservationRect(Frame& CurrentFrame, vector<int> KeysToRect)
	{
		int N = CurrentFrame.N;
		vector<bool> pOutlier1(CurrentFrame.N, false);
		this->pOutlier.clear();
		this->pOutlier = pOutlier1;
		for (size_t i = 0; i < N; i++){
			MapPoint* pMP = CurrentFrame.mvpMapPoints[i];
			if (pMP){
				float DiffX = abs(vp2fProjections[KeysToRect[i]].x - CurrentFrame.mvKeys[i].pt.x);
				float DiffY = abs(vp2fProjections[KeysToRect[i]].y - CurrentFrame.mvKeys[i].pt.y);
				//cout << DiffX << "\t" << DiffY << endl;
				if (DiffX > DiffXTreshold || DiffY > DiffYTreshold)
				{
					pOutlier[i] = true;
					OutlierCounter++;
					continue;
				}
				else
				{
					InlierCounter++;
					//SumDiff += DiffX*DiffX + DiffY*DiffY;
					SumDiffX += DiffX*DiffX;
					SumDiffY += DiffY*DiffY;
				}

			}
		}
		float AverageSumDiffX = SumDiffX / float(InlierCounter);
		float AverageSumDiffY = SumDiffY / float(InlierCounter);

		Likelihood = InlierCounter * (1 / AverageSumDiffX)*(1 / AverageSumDiffY);
		Likelihood = Likelihood*Likelihood*Likelihood;
		if (InlierCounter < 10)
			Likelihood = 0;
	}

	void Particle::ObservationRectGrid(Frame* CurrentFrame, vector<int> KeysToRect, vector<int> ProjectionWeight)
	{
		int N = CurrentFrame->N;
		vector<bool> pOutlier1(CurrentFrame->N, false);
		this->pOutlier.clear();
		this->pOutlier = pOutlier1;
		for (size_t i = 0; i < N; i++){
			MapPoint* pMP = CurrentFrame->mvpMapPoints[i];
			if (pMP){
				float DiffX = abs(vp2fProjections[KeysToRect[i]].x - CurrentFrame->mvKeys[i].pt.x)*ProjectionWeight[KeysToRect[i]];
				float DiffY = abs(vp2fProjections[KeysToRect[i]].y - CurrentFrame->mvKeys[i].pt.y)*ProjectionWeight[KeysToRect[i]];
				//cout << DiffX << "\t" << DiffY << endl;
				if (DiffX > DiffXTreshold || DiffY > DiffYTreshold)
				{
					pOutlier[i] = true;
					OutlierCounter++;
					continue;
				}
				else
				{
					InlierCounter++;
					//SumDiff += DiffX*DiffX + DiffY*DiffY;
					SumDiffX += DiffX*DiffX;
					SumDiffY += DiffY*DiffY;
				}

			}
		}
		if (InlierCounter < 5 || SumDiffX == 0 || SumDiffY == 0)
			Likelihood = 0;
		else
		{
			float AverageSumDiffX = SumDiffX / float(InlierCounter);
			float AverageSumDiffY = SumDiffY / float(InlierCounter);

			Likelihood = InlierCounter * (1 / AverageSumDiffX)*(1 / AverageSumDiffY);
			Likelihood = Likelihood*Likelihood*Likelihood;
		}
	}

	void Particle::ObservationRectGrid(Frame* CurrentFrame, vector<int> &KeysToRect)
	{
		int N = CurrentFrame->N;
		vector<bool> pOutlier1(CurrentFrame->N, false);
		this->pOutlier.clear();
		this->pOutlier = pOutlier1;
		for (size_t i = 0; i < N; i++){
			MapPoint* pMP = CurrentFrame->mvpMapPoints[i];
			if (pMP){
				float DiffX = abs(vp2fProjections[KeysToRect[i]].x - CurrentFrame->mvKeys[i].pt.x);
				float DiffY = abs(vp2fProjections[KeysToRect[i]].y - CurrentFrame->mvKeys[i].pt.y);
				//cout << DiffX << "\t" << DiffY << endl;
				if (DiffX > DiffXTreshold || DiffY > DiffYTreshold)
				{
					pOutlier[i] = true;
					OutlierCounter++;
					continue;
				}
				else
				{
					InlierCounter++;
					//SumDiff += DiffX*DiffX + DiffY*DiffY;
					SumDiffX += DiffX*DiffX;
					SumDiffY += DiffY*DiffY;
				}

			}
		}
		if (InlierCounter < 5 || SumDiffX == 0 || SumDiffY == 0)
			Likelihood = 0;
		else
		{
			float AverageSumDiffX = SumDiffX / float(InlierCounter);
			float AverageSumDiffY = SumDiffY / float(InlierCounter);

			Likelihood = InlierCounter * (1 / AverageSumDiffX)*(1 / AverageSumDiffY);
			Likelihood = Likelihood*Likelihood*Likelihood;
		}
	}

	void Particle::Projection(Frame& CurrentFrame,Frame& LastFrame, vector<pair<cv::Point2f, cv::Point2f>> &Rectangle)
	{
		int N = LastFrame.N;
		vp2fProjections.clear();
		vp2fProjections.resize(N, cv::Point2f(0, 0));
		const cv::Mat Rcw = Pose.RT.rowRange(0, 3).colRange(0, 3);
		const cv::Mat tcw = Pose.RT.rowRange(0, 3).col(3);

		//cout << "Particle: start projection" << endl;
		for (int i_MP = 0; i_MP < LastFrame.N; i_MP++)
		{
			MapPoint* pMP = LastFrame.mvpMapPoints[i_MP];
			if (pMP)
			{
				if (!LastFrame.mvbOutlier[i_MP])
				{
					cv::Mat x3Dw = pMP->GetWorldPos();
					cv::Mat x3Dc = Rcw*x3Dw + tcw;

					const float xc = x3Dc.at<float>(0);
					const float yc = x3Dc.at<float>(1);
					const float invzc = 1.0 / x3Dc.at<float>(2);

					if (invzc < 0)
					{
						//cout << "invzc < 0" << endl;
						continue;	
					}


					float u = CurrentFrame.fx*xc*invzc + CurrentFrame.cx;
					float v = CurrentFrame.fy*yc*invzc + CurrentFrame.cy;

					if (u<CurrentFrame.mnMinX || u>CurrentFrame.mnMaxX)
					{
						//cout << "Projection out of range x" << endl;
						continue;
					}

					if (v<CurrentFrame.mnMinY || v>CurrentFrame.mnMaxY)
					{
						//cout << "Projection out of range y" << endl;
						continue;
					}

					vp2fProjections[i_MP].x = u;
					vp2fProjections[i_MP].y = v;

					if (vp2fProjections[i_MP].x < Rectangle[i_MP].first.x)
						Rectangle[i_MP].first.x = vp2fProjections[i_MP].x;
					if (vp2fProjections[i_MP].y < Rectangle[i_MP].first.y)
						Rectangle[i_MP].first.y = vp2fProjections[i_MP].y;
					if (vp2fProjections[i_MP].x > Rectangle[i_MP].second.x)
						Rectangle[i_MP].second.x = vp2fProjections[i_MP].x;
					if (vp2fProjections[i_MP].y > Rectangle[i_MP].second.y)
						Rectangle[i_MP].second.y = vp2fProjections[i_MP].y;
				}
			}
		}
	}

	void Particle::Projection(vector<float> &Indices, const vector<cv::Point3f> &x3Dw, vector<pair<cv::Point2f, cv::Point2f>> &rectangle, vector<int> &ProjectionIdx, const int& N)
	{
		vp2fProjections.clear();
		vp2fProjections.resize(N, cv::Point2f(0, 0));
		for (int i_MP = 0, i_MP_end = ProjectionIdx.size(); i_MP <i_MP_end; ++i_MP)
			ProjectionOneFixedPoint(Indices, &x3Dw[i_MP], rectangle, ProjectionIdx[i_MP]);
			//ProjectionOne(Indices, &x3Dw[i_MP], rectangle, ProjectionIdx[i_MP]);
	}

	void Particle::Projection2(vector<float> &Indices, const vector<cv::Point3f> &x3Dw, vector<int> &ProjectionIdx, const int& N)
	{
		vp2fProjections.clear();
		vp2fProjections.resize(N, cv::Point2f(0, 0));
		for (int i_MP = 0, i_MP_end = ProjectionIdx.size(); i_MP < i_MP_end; ++i_MP)
			ProjectionOne(Indices, &x3Dw[i_MP], ProjectionIdx[i_MP]);
	}

	inline bool Particle::ProjectionOne(vector<float> & Indices, const cv::Point3f *x3Dw, vector<pair<cv::Point2f, cv::Point2f>> &Rectangle, int &ProjectionIdx)
	{
		// step 1
		const float xc = Rcw.at<float>(0, 0)*x3Dw->x + Rcw.at<float>(0, 1)*x3Dw->y + Rcw.at<float>(0, 2)*x3Dw->z + tcw.at<float>(0);
		const float yc = Rcw.at<float>(1, 0)*x3Dw->x + Rcw.at<float>(1, 1)*x3Dw->y + Rcw.at<float>(1, 2)*x3Dw->z + tcw.at<float>(1);
		const float invzc = 1.0 / (Rcw.at<float>(2, 0)*x3Dw->x + Rcw.at<float>(2, 1)*x3Dw->y + Rcw.at<float>(2, 2)*x3Dw->z + tcw.at<float>(2));//this writing is faster because directly visit memory

		if (!isfinite(invzc))
			cout << "invzc is too not finite!" << endl;
		// step 2
		if (invzc < 0)
		{
			//cout << "invzc < 0" << endl;
			return false;
		}


		float u = Indices[0] * xc*invzc + Indices[2];
		float v = Indices[1] * yc*invzc + Indices[3];

		if (u<Indices[4] || u>Indices[6])
		{
			//cout << "Projection out of range x" << endl;
			return false;
		}

		if (v<Indices[5] || v>Indices[7])
		{
			//cout << "Projection out of range y" << endl;
			return false;
		}


		// step 3
		vp2fProjections[ProjectionIdx].x = u;
		vp2fProjections[ProjectionIdx].y = v;

		if (vp2fProjections[ProjectionIdx].x < Rectangle[ProjectionIdx].first.x)
			Rectangle[ProjectionIdx].first.x = vp2fProjections[ProjectionIdx].x;
		if (vp2fProjections[ProjectionIdx].y < Rectangle[ProjectionIdx].first.y)
			Rectangle[ProjectionIdx].first.y = vp2fProjections[ProjectionIdx].y;
		if (vp2fProjections[ProjectionIdx].x > Rectangle[ProjectionIdx].second.x)
			Rectangle[ProjectionIdx].second.x = vp2fProjections[ProjectionIdx].x;
		if (vp2fProjections[ProjectionIdx].y > Rectangle[ProjectionIdx].second.y)
			Rectangle[ProjectionIdx].second.y = vp2fProjections[ProjectionIdx].y;
	}

	inline bool Particle::ProjectionOne(vector<float> & Indices, const cv::Point3f *x3Dw, int &ProjectionIdx)
	{
		// step 1
		const float xc = Rcw.at<float>(0, 0)*x3Dw->x + Rcw.at<float>(0, 1)*x3Dw->y + Rcw.at<float>(0, 2)*x3Dw->z + tcw.at<float>(0);
		const float yc = Rcw.at<float>(1, 0)*x3Dw->x + Rcw.at<float>(1, 1)*x3Dw->y + Rcw.at<float>(1, 2)*x3Dw->z + tcw.at<float>(1);
		float temp_invzc = Rcw.at<float>(2, 0)*x3Dw->x + Rcw.at<float>(2, 1)*x3Dw->y + Rcw.at<float>(2, 2)*x3Dw->z + tcw.at<float>(2);
		float invzc;
		PiecewiseDivision(temp_invzc, invzc);
		//const float invzc = 1.0 / (Rcw.at<float>(2, 0)*x3Dw->x + Rcw.at<float>(2, 1)*x3Dw->y + Rcw.at<float>(2, 2)*x3Dw->z + tcw.at<float>(2));//this writing is faster because directly visit memory

		if (!isfinite(invzc))
			cout << "invzc is too not finite!" << endl;
		// step 2
		if (invzc < 0)
		{
			//cout << "invzc < 0" << endl;
			return false;
		}


		float u = Indices[0] * xc*invzc + Indices[2];
		float v = Indices[1] * yc*invzc + Indices[3];

		if (u<Indices[4] || u>Indices[6])
		{
			//cout << "Projection out of range x" << endl;
			return false;
		}

		if (v<Indices[5] || v>Indices[7])
		{
			//cout << "Projection out of range y" << endl;
			return false;
		}


		// step 3
		vp2fProjections[ProjectionIdx].x = u;
		vp2fProjections[ProjectionIdx].y = v;

	}

	inline bool Particle::ProjectionOneFixedPoint(vector<float> & Indices, const cv::Point3f *x3Dw, vector<pair<cv::Point2f, cv::Point2f>> &Rectangle, int &ProjectionIdx)
	{
		cv::Mat RcwFixed = Rcw.clone();
		for (int i = 0; i < RcwFixed.cols; i++)
			for (int j = 0; j < RcwFixed.rows; j++)
				FixedPointTrans(Rcw.at<float>(i, j), RcwFixed.at<float>(i, j));
		cv::Mat tcwFixed = tcw.clone();
		FixedPointTrans(tcw.at<float>(0), tcwFixed.at<float>(0));
		FixedPointTrans(tcw.at<float>(1), tcwFixed.at<float>(1));
		FixedPointTrans(tcw.at<float>(2), tcwFixed.at<float>(2));

		cv::Point3f x3DwFixed;
		float x, y, z;
		x = x3Dw->x;
		y = x3Dw->y;
		z = x3Dw->z;

		FixedPointTrans(x, x3DwFixed.x);
		FixedPointTrans(y, x3DwFixed.y);
		FixedPointTrans(z, x3DwFixed.z);
		// step 1
		const float xc = RcwFixed.at<float>(0, 0)*x3DwFixed.x + RcwFixed.at<float>(0, 1)*x3DwFixed.y + RcwFixed.at<float>(0, 2)*x3DwFixed.z + tcwFixed.at<float>(0);
		const float yc = RcwFixed.at<float>(1, 0)*x3DwFixed.x + RcwFixed.at<float>(1, 1)*x3DwFixed.y + RcwFixed.at<float>(1, 2)*x3DwFixed.z + tcwFixed.at<float>(1);
		float invzctemp = (RcwFixed.at<float>(2, 0)*x3DwFixed.x + RcwFixed.at<float>(2, 1)*x3DwFixed.y + RcwFixed.at<float>(2, 2)*x3DwFixed.z + tcwFixed.at<float>(2));//this writing is faster because directly visit memory
		float invzc;
		PiecewiseDivision(invzctemp, invzc);
		//cout << "linear division test: "<< z << "\t" << invzctemp << "\t" << invzc << "\t" << (double)(1 / invzctemp) << endl;
		if (!isfinite(invzc))
			cout << "invzc is not finite!" << endl;
		// step 2
		if (invzc < 0)
		{
			//cout << "invzc < 0" << endl;
			return false;
		}


		float u = Indices[0] * xc*invzc + Indices[2];
		float v = Indices[1] * yc*invzc + Indices[3];

		if (u<Indices[4] || u>Indices[6])
		{
			//cout << "Projection out of range x" << endl;
			return false;
		}

		if (v<Indices[5] || v>Indices[7])
		{
			//cout << "Projection out of range y" << endl;
			return false;
		}


		// step 3
		vp2fProjections[ProjectionIdx].x = u;
		vp2fProjections[ProjectionIdx].y = v;

		if (vp2fProjections[ProjectionIdx].x < Rectangle[ProjectionIdx].first.x)
			Rectangle[ProjectionIdx].first.x = vp2fProjections[ProjectionIdx].x;
		if (vp2fProjections[ProjectionIdx].y < Rectangle[ProjectionIdx].first.y)
			Rectangle[ProjectionIdx].first.y = vp2fProjections[ProjectionIdx].y;
		if (vp2fProjections[ProjectionIdx].x > Rectangle[ProjectionIdx].second.x)
			Rectangle[ProjectionIdx].second.x = vp2fProjections[ProjectionIdx].x;
		if (vp2fProjections[ProjectionIdx].y > Rectangle[ProjectionIdx].second.y)
			Rectangle[ProjectionIdx].second.y = vp2fProjections[ProjectionIdx].y;
	}

	bool Particle::ProjectionOne1(vector<float> & Indices, const cv::Point3f &x3Dw, vector<pair<cv::Point2f, cv::Point2f>> &Rectangle, int &ProjectionIdx)
	{
		
		//vx3Dc[ProjectionIdx] = Rcw*x3Dw + tcw;
		//cv::Mat x3Dc = cv::Mat(3, 1, CV_32F, 0.0);
		/*for (int i = 0; i < 3; i++)
			vx3Dc[ProjectionIdx].at<float>(i) =
			Rcw.at<float>(i, 0)*x3Dw.x + 
			Rcw.at<float>(i, 1)*x3Dw.y + 
			Rcw.at<float>(i, 2)*x3Dw.z + 
			tcw.at<float>(i);*/
		vx3Dc[ProjectionIdx].at<float>(0) = Rcw.at<float>(0, 0)*x3Dw.x + Rcw.at<float>(0, 1)*x3Dw.y + Rcw.at<float>(0, 2)*x3Dw.z + tcw.at<float>(0);
		vx3Dc[ProjectionIdx].at<float>(1) = Rcw.at<float>(1, 0)*x3Dw.x + Rcw.at<float>(1, 1)*x3Dw.y + Rcw.at<float>(1, 2)*x3Dw.z + tcw.at<float>(1);
		vx3Dc[ProjectionIdx].at<float>(2) = Rcw.at<float>(2, 0)*x3Dw.x + Rcw.at<float>(2, 1)*x3Dw.y + Rcw.at<float>(2, 2)*x3Dw.z + tcw.at<float>(2);

		vxc[ProjectionIdx] = vx3Dc[ProjectionIdx].at<float>(0);
		vyc[ProjectionIdx] = vx3Dc[ProjectionIdx].at<float>(1);
		vinvzc[ProjectionIdx] = 1.0 / vx3Dc[ProjectionIdx].at<float>(2);

		return true;
	}

	bool Particle::ProjectionOne2(vector<float> & Indices, const cv::Point3f &x3Dw, vector<pair<cv::Point2f, cv::Point2f>> &Rectangle, int &ProjectionIdx)
	{

		/*vxc[ProjectionIdx] = vx3Dc[ProjectionIdx].at<float>(0);
		vyc[ProjectionIdx] = vx3Dc[ProjectionIdx].at<float>(1);
		vinvzc[ProjectionIdx] = 1.0 / vx3Dc[ProjectionIdx].at<float>(2);*/

		if (vinvzc[ProjectionIdx] < 0)
		{
			//cout << "invzc < 0" << endl;
			return false;
		}


		float u = Indices[0] * vxc[ProjectionIdx] * vinvzc[ProjectionIdx] + Indices[2];
		float v = Indices[1] * vyc[ProjectionIdx] * vinvzc[ProjectionIdx] + Indices[3];

		if (u<Indices[4] || u>Indices[6])
		{
			//cout << "Projection out of range x" << endl;
			return false;
		}

		if (v<Indices[5] || v>Indices[7])
		{
			//cout << "Projection out of range y" << endl;
			return false;
		}

		vp2fProjections[ProjectionIdx].x = u;
		vp2fProjections[ProjectionIdx].y = v;
	}

	bool Particle::ProjectionOne3(vector<float> & Indices, const cv::Point3f &x3Dw, vector<pair<cv::Point2f, cv::Point2f>> &Rectangle, int &ProjectionIdx)
	{
		if (vp2fProjections[ProjectionIdx].x < Rectangle[ProjectionIdx].first.x)
			Rectangle[ProjectionIdx].first.x = vp2fProjections[ProjectionIdx].x;
		if (vp2fProjections[ProjectionIdx].y < Rectangle[ProjectionIdx].first.y)
			Rectangle[ProjectionIdx].first.y = vp2fProjections[ProjectionIdx].y;
		if (vp2fProjections[ProjectionIdx].x > Rectangle[ProjectionIdx].second.x)
			Rectangle[ProjectionIdx].second.x = vp2fProjections[ProjectionIdx].x;
		if (vp2fProjections[ProjectionIdx].y > Rectangle[ProjectionIdx].second.y)
			Rectangle[ProjectionIdx].second.y = vp2fProjections[ProjectionIdx].y;

		return true;
	}

	inline void Particle::PiecewiseDivision(float &input, float& output)
	{
		int piecewise = 0;
		if (input < PiecewiseInterval[0])
		{
			output = 0;
			return;
		}
		if (input > PiecewiseInterval[24])
		{
			output = 0;
			return;
		}
		for (int i = 1; i < 24; i++)
			if (input >= PiecewiseInterval[i] && input < PiecewiseInterval[i + 1])
				piecewise = i;
		//cout << piecewise << "\t";
		output = input * Piecewisek[piecewise] + Piecewiseb[piecewise];
	}

	inline void Particle::FixedPointTrans(float &input, float& output)
	{
		int a = input;
		output = (double)((int)(input * FixedPointBitNumber) % FixedPointBitNumber) / FixedPointBitNumber;
		output += a;
	}

	void Particle::Observation1(Frame* CurrentFrame)
	{
		int N = CurrentFrame->N;
		vector<bool> pOutlier1(CurrentFrame->N, false);
		this->pOutlier.clear();
		this->pOutlier = pOutlier1;
		vMx3Dc.clear();
		vMx3Dc.resize(N);
		for (size_t i = 0; i < N; i++){
			MapPoint* pMP = CurrentFrame->mvpMapPoints[i];
			if (pMP){
				//step 1
				const cv::Mat Rcw = this->Pose.RT.rowRange(0, 3).colRange(0, 3);
				const cv::Mat tcw = this->Pose.RT.rowRange(0, 3).col(3);

				cv::Mat x3Dw = pMP->GetWorldPos();
				cv::Mat x3Dc = Rcw*x3Dw + tcw;

				vMx3Dc[i] = x3Dc.clone();
			}
		}
	}

	void Particle::Observation11(Frame* CurrentFrame)
	{
		int N = CurrentFrame->N;
		vector<bool> pOutlier1(CurrentFrame->N, false);
		this->pOutlier.clear();
		this->pOutlier = pOutlier1;
		vMRcw.clear();
		vMRcw.resize(N);
		vMtcw.clear();
		vMtcw.resize(N);
		for (size_t i = 0; i < N; i++){
			MapPoint* pMP = CurrentFrame->mvpMapPoints[i];
			if (pMP){
				//step 1-1
				const cv::Mat Rcw = this->Pose.RT.rowRange(0, 3).colRange(0, 3);
				const cv::Mat tcw = this->Pose.RT.rowRange(0, 3).col(3);

				vMRcw[i] = Rcw.clone();
				vMtcw[i] = tcw.clone();
			}
		}
	}

	void Particle::Observation12(Frame* CurrentFrame)
	{
		int N = CurrentFrame->N;
		vector<bool> pOutlier1(CurrentFrame->N, false);
		this->pOutlier.clear();
		this->pOutlier = pOutlier1;
		vMx3Dw.clear();
		vMx3Dw.resize(N);
		for (size_t i = 0; i < N; i++){
			MapPoint* pMP = CurrentFrame->mvpMapPoints[i];
			if (pMP){
				//step 1-2
				cv::Mat x3Dw = pMP->GetWorldPos();
				vMx3Dw[i] = x3Dw.clone();
			}
		}
	}

	void Particle::Observation13(Frame* CurrentFrame)
	{
		int N = CurrentFrame->N;
		vector<bool> pOutlier1(CurrentFrame->N, false);
		this->pOutlier.clear();
		this->pOutlier = pOutlier1;
		vMx3Dc.clear();
		vMx3Dc.resize(N);
		for (size_t i = 0; i < N; i++){
			MapPoint* pMP = CurrentFrame->mvpMapPoints[i];
			if (pMP){
				//step 1-3
				//cv::Mat x3Dc = (cv::Mat_<float>(3, 1) << 0, 0, 0);
				//x3Dc.at<float>(0, 0) = vMRcw[i].at<float>(0, 0)*vMx3Dw[i].at<float>(0, 0) + vMRcw[i].at<float>(0, 1)*vMx3Dw[i].at<float>(1, 0) + vMRcw[i].at<float>(0, 2)*vMx3Dw[i].at<float>(2, 0) + vMtcw[i].at<float>(0, 0);
				//x3Dc.at<float>(1, 0) = vMRcw[i].at<float>(1, 0)*vMx3Dw[i].at<float>(1, 0) + vMRcw[i].at<float>(1, 1)*vMx3Dw[i].at<float>(1, 0) + vMRcw[i].at<float>(1, 2)*vMx3Dw[i].at<float>(2, 0) + vMtcw[i].at<float>(1, 0);
				//x3Dc.at<float>(2, 0) = vMRcw[i].at<float>(2, 0)*vMx3Dw[i].at<float>(2, 0) + vMRcw[i].at<float>(2, 1)*vMx3Dw[i].at<float>(1, 0) + vMRcw[i].at<float>(2, 2)*vMx3Dw[i].at<float>(2, 0) + vMtcw[i].at<float>(2, 0);
				cv::Mat x3Dc = vMRcw[i] * vMx3Dw[i] + vMtcw[i];
				vMx3Dc[i] = x3Dc.clone();
			}
		}
	}

	void Particle::Observation2(Frame* CurrentFrame)
	{
		int N = CurrentFrame->N;
		vp2fProjections.clear();
		vp2fProjections.resize(N, cv::Point2f(0, 0));
		for (size_t i = 0; i < N; i++){
			MapPoint* pMP = CurrentFrame->mvpMapPoints[i];
			if (pMP){
				cv::Point2f p2fProjection;
				const float xc = vMx3Dc[i].at<float>(0);
				const float yc = vMx3Dc[i].at<float>(1);
				const float invzc = 1.0 / vMx3Dc[i].at<float>(2);

				if (invzc < 0)
				{
					pOutlier[i] = true;
					OutlierCounter++;
					continue;
				}


				float u = CurrentFrame->fx*xc*invzc + CurrentFrame->cx;
				float v = CurrentFrame->fy*yc*invzc + CurrentFrame->cy;

				p2fProjection.x = u;//Project from local 3D coordinate to 2D coordinate
				p2fProjection.y = v;

				vp2fProjections[i] = p2fProjection;
			}
		}
	}

	void Particle::Observation3(Frame* CurrentFrame)
	{
		int N = CurrentFrame->N;
		for (size_t i = 0; i < N; i++){
			MapPoint* pMP = CurrentFrame->mvpMapPoints[i];
			if (pMP){
				float DiffX = abs(vp2fProjections[i].x - CurrentFrame->mvKeys[i].pt.x);
				float DiffY = abs(vp2fProjections[i].y - CurrentFrame->mvKeys[i].pt.y);
				//cout << DiffX << "\t" << DiffY << endl;
				if (DiffX > DiffXTreshold || DiffY > DiffYTreshold)
				{
					pOutlier[i] = true;
					OutlierCounter++;
					continue;
				}
				else
				{
					InlierCounter++;
					//SumDiff += DiffX*DiffX + DiffY*DiffY;
					SumDiffX += DiffX*DiffX;
					SumDiffY += DiffY*DiffY;
				}
			}
		}
		float AverageSumDiffX = SumDiffX / float(InlierCounter);
		float AverageSumDiffY = SumDiffY / float(InlierCounter);

		Likelihood = InlierCounter * (1 / AverageSumDiffX)*(1 / AverageSumDiffY);
		//Likelihood = 1 / (SumDiff + 1);
		//Likelihood = InlierCounter*(1 / SumDiffX)*(1 / SumDiffX);
		Likelihood = Likelihood*Likelihood*Likelihood;

	}

	void Particle::Clone(Particle *CloneSource){
		Pose.Clone(CloneSource->Pose);
		Speedvt.Clone(CloneSource->Speedvt);
		ParticleSpeed = CloneSource->ParticleSpeed;
	}

	Eigen::Matrix<double, 3, 3> Particle::toMatrix3d(const cv::Mat &cvMat3)
	{
		Eigen::Matrix<double, 3, 3> M;

		M << cvMat3.at<float>(0, 0), cvMat3.at<float>(0, 1), cvMat3.at<float>(0, 2),
			cvMat3.at<float>(1, 0), cvMat3.at<float>(1, 1), cvMat3.at<float>(1, 2),
			cvMat3.at<float>(2, 0), cvMat3.at<float>(2, 1), cvMat3.at<float>(2, 2);

		return M;
	}

	void Particle::Reset(){
		fill(ParticleAcceleration.begin(),ParticleAcceleration.end(), 0.0);
		Likelihood = 0.0;
		Weight = 0.0;
		WeightSum = 0.0;
		SumDiff = 0.0;
		SumDiffX = 0.0;
		SumDiffY = 0.0;
		OutlierCounter = 0;
		InlierCounter = 0;
	}
}//namespace ORB_SLAM2