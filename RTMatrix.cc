#include "RTMatrix.h"

using namespace std;

namespace ORB_SLAM2
{
	RTMatrix::RTMatrix() :Idx(6, 0)
	{
		RT = (cv::Mat_<float>(4, 4) <<
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1);
	}


	void RTMatrix::Mat_to_Idx()
	{
		float sy = sqrt(RT.at<float>(0, 0)*RT.at<float>(0, 0) + RT.at<float>(1, 0)*RT.at<float>(1, 0));//这种欧拉角的取法，是围绕世界坐标系的zyx轴的旋转，并不是通常的百度百科里欧拉角的取法
		bool singular = sy < 1e-6;
		float euler_x, euler_y, euler_z;
		if (!singular)
		{
			euler_x = atan2(RT.at<float>(2, 1), RT.at<float>(2, 2));
			euler_y = atan2(-RT.at<float>(2, 0), sy);
			euler_z = atan2(RT.at<float>(1, 0), RT.at<float>(0, 0));
		}
		else
		{
			euler_x = atan2(-RT.at<float>(1, 2), RT.at<float>(1, 1));
			euler_y = atan2(-RT.at<float>(2, 0), sy);
			euler_z = 0;
		}
		Idx[0] = euler_x;
		if (!isfinite(euler_x))
			cout << "euler_x is not finite" << endl;
		Idx[1] = euler_y;
		if (!isfinite(euler_x))
			cout << "euler_y is not finite" << endl;
		Idx[2] = euler_z;
		if (!isfinite(euler_x))
			cout << "euler_z is not finite" << endl;
		Idx[3] = RT.at<float>(0, 3);
		if (!isfinite(euler_x))
			cout << "x is not finite" << endl;
		Idx[4] = RT.at<float>(1, 3);
		if (!isfinite(euler_x))
			cout << "y is not finite" << endl;
		Idx[5] = RT.at<float>(2, 3);
		if (!isfinite(euler_x))
			cout << "z is not finite" << endl;
	}

	void RTMatrix::Idx_to_Mat()
	{
		cv::Mat R_x = (cv::Mat_<float>(3, 3) <<
			1, 0, 0,
			0, cos(Idx[0]), -sin(Idx[0]),
			0, sin(Idx[0]), cos(Idx[0])
			);

		// Calculate rotation about y axis
		cv::Mat R_y = (cv::Mat_<float>(3, 3) <<
			cos(Idx[1]), 0, sin(Idx[1]),
			0, 1, 0,
			-sin(Idx[1]), 0, cos(Idx[1])
			);

		// Calculate rotation about z axis
		cv::Mat R_z = (cv::Mat_<float>(3, 3) <<
			cos(Idx[2]), -sin(Idx[2]), 0,
			sin(Idx[2]), cos(Idx[2]), 0,
			0, 0, 1);

		cv::Mat R = R_z*R_y*R_x;

		R.copyTo(RT.rowRange(0, 3).colRange(0, 3));
		RT.at<float>(0, 3) = Idx[3];
		RT.at<float>(1, 3) = Idx[4];
		RT.at<float>(2, 3) = Idx[5];


	}

	void RTMatrix::Clone(RTMatrix CloneSource)
	{
		RT = CloneSource.RT.clone();
		Idx = CloneSource.Idx;

	}

	void RTMatrix::Clear()
	{
		RT = (cv::Mat_<float>(4, 4) <<
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1);

		fill(Idx.begin(), Idx.end(), 0.0);
	}

	RTMatrix::~RTMatrix(){}

}