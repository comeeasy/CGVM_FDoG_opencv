#pragma once

#ifndef __UTILS_H__
#define __UTILS_H__

#include <iostream>
#include <string>
#include <opencv2/opencv.hpp>

#include "FlowPath.h"

namespace utils {

enum IMG_IDX {ORIGIN, INFODRAW, GRAD, TANGENT, ETF, CL, FBL};

cv::Mat     g_monitor_image;
cv::Mat*    g_images_for_monitor;
FlowPath**  g_fpath; 
cv::Point   g_ptOld;
cv::Point   g_ptNew;
int         g_toggleFill = 1;
int         g_image_idx = 0;

void draw_RGB_image(cv::Mat RGB_image, std::string title);

void interactive_monitor(cv::Mat* imgs, FlowPath** fpath);
void on_mouse(int event, int x, int y, int flags, void*);

void draw_fpath(int x, int y);
}



#endif