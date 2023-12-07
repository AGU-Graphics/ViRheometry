#ifndef main_h
#define main_h

#include <stdio.h>
#include <Eigen/Core>
#include <Eigen/Dense>

struct UI
{
	Eigen::Vector2d center;
	Eigen::Vector2d width;
	int window_width;
	int window_height;
	double extend_ratio;
	
	int mx, my;
};
extern UI g_UI;

#endif /* main_h */
