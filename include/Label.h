#pragma once

#ifndef __LABEL_H__
#define __LABEL_H__

#include "mtoolkits.h"

//라벨링 구조체
struct LABEL{
	int num_pixel;	//픽셀 개수
	V3DF clr;		//대표색(Luv)
	V2DI pos;		//대표위치(임의)
public:
	LABEL& operator = (const LABEL& label){
		num_pixel = label.num_pixel;
		clr[0] = label.clr[0];	clr[1] = label.clr[1];	clr[2] = label.clr[2];
		pos[0] = label.pos[0];	pos[1] = label.pos[1];
		return *this;
	}
};

#endif