#pragma once

#ifndef __FLOWPATH_H__
#define __FLOWPATH_H__

#include "mtoolkits.h"


class FlowPath{
public:
	int sn;	// s ???? ????
	V3DF *Alpha;
	int tn;	// t ???? ????
	V3DF *Beta;
	V3DF *AlphaMax;	// Alpha -> BetaMax?? ??? ?????

	bool is_silhouette;	// silhouette ????
	float importance;	// ??? ????? stream ????

	// for vectorization
	V3DF BBox[2];
	float skltn, cl;
	float length;
	float nlength;	// normalized length
	float tcurv;	//	threshold curvature: negative ????? positive ???????? tolerable?? curvature ????
	float *curv;	//	??? pixel???? ????? curvature
	V2DI *crange;	//	negative ???????? tolerable?? curvature?? ???? pixel?? index ~ positive ?????????? index
	int selected_segment;	// selected index; ?? index?? ???????? +,- ???? crange[selected_segment][] ????

	FlowPath(){
		sn = tn = 0;
		Alpha = NULL;
		Beta = NULL;
		AlphaMax = NULL;
		is_silhouette = false;	importance = 0.0f;

		skltn = cl = length = nlength = 0.0f;
		tcurv = 0;	curv = NULL;	crange = NULL;	selected_segment = 0;
	}
	~FlowPath(){
		if (sn>0)	delete [] Alpha;
		if (tn>0)	delete [] Beta;
		if (!AlphaMax)	delete [] AlphaMax;
		if (curv)	delete [] curv;
		if (crange)	delete [] crange;
	}
	void Init(){
		if ( !Alpha )	delete [] Alpha;
		if ( !Beta )	delete [] Beta;
		if ( !AlphaMax)	delete [] AlphaMax;
		Alpha = NULL;		Beta = NULL;	AlphaMax = NULL;
		sn = tn = 0;
		is_silhouette = false;	importance = 0.0f;

		skltn = cl = length = nlength = 0.0f;
		tcurv = 0;
		if (curv) {
			delete[] curv;
		}		
		curv = NULL;
		if (crange)	{delete[] crange;}	crange = NULL;
		selected_segment = 0;
	}
	void InitShallow(){
		sn = tn = 0;	
		is_silhouette = false;	importance = 0.0f;
		skltn = cl = length = nlength = 0.0f;
		tcurv = 0;	selected_segment = 0;
	}
	void get_length ( ){
		int i;
		this->length = 0.0f;

		for ( i = 0; i < this->sn - 1; i++ )
			this->length += pdist ( this->Alpha[i], this->Alpha[i+1] );
	}
	void Copy(FlowPath *p){
		Init();

		int i;
		sn = p->sn;
		if ( sn>0 ){
			Alpha = new V3DF[sn];
			for (i=0; i<sn; i++)	vcopy(Alpha[i], p->Alpha[i]);
		}

		tn = p->tn;
		if ( tn>0 ){
			Beta = new V3DF[tn];
			for (i=0; i<tn; i++)	vcopy(Beta[i], p->Beta[i]);
		}
	}
};

#endif