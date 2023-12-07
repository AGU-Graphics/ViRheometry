#ifndef FluidMat_h
#define FluidMat_h

#include <stdio.h>
#include <opencv2/opencv.hpp>

struct FluidMat final
{
	FluidMat();
	cv::Mat fluid_mat_cal;
	cv::Mat fluid_mat_ref;
	cv::Mat Image;
	
	double cal_ref_diff;
	int sum_count;
	std::vector<double> loss_values;
};


#endif /* FluidMat_h */
