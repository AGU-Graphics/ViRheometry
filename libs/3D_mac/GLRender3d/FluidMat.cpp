#include "FluidMat.h"
#include "main.h"

FluidMat::FluidMat()
{
	fluid_mat_cal = cv::Mat::zeros( g_UI.window_height, g_UI.window_width, CV_8UC3 );
	fluid_mat_ref = cv::Mat::zeros( g_UI.window_height, g_UI.window_width, CV_8UC3 );
	cal_ref_diff = 0;
	sum_count = 0;
	loss_values = {};
}

