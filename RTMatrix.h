#ifndef _RT_Matrix_H_
#define _RT_Matrix_H_

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/objdetect/objdetect.hpp>
#include <iostream>
#include <vector>
#include <string>


namespace ORB_SLAM2
{
class RTMatrix
{
public:
	cv::Mat RT;// projection matrix, including rotation and translation
	std::vector<float> Idx;

	RTMatrix();

	void Mat_to_Idx();
	void Idx_to_Mat();
	void Clone(RTMatrix CloneSource);

	void Clear();

	~RTMatrix();
};
}//namespace ORB_SLAM2


#endif