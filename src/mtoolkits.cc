#include "mtoolkits.h"



double max_d ( int n, V2DD *list )
{
	double maxval;
	int i;

	for ( i = 1, maxval = list[0][1]; i < n; i++ ) {
		if ( list[i][1] > maxval )
			maxval = list[i][1];
	}

	return maxval;
}

double min_d ( int n, V2DD *list )
{
	double minval;
	int i;

	for ( i = 1, minval = list[0][1]; i < n; i++ ) {
		if ( list[i][1] < minval )
			minval = list[i][1];
	}

	return minval;
}

double max_d ( int n, V3DD *list )
{
	double maxval;
	int i;

	for ( i = 1, maxval = list[0][1]; i < n; i++ ) {
		if ( list[i][1] > maxval )
			maxval = list[i][1];
	}

	return maxval;
}

double min_d ( int n, V3DD *list )
{
	double minval;
	int i;

	for ( i = 1, minval = list[0][1]; i < n; i++ ) {
		if ( list[i][1] < minval )
			minval = list[i][1];
	}

	return minval;
}

float avg3f ( float a, float b, float c )
{
	return (a+b+c)/3.0f;
}

float avg2f ( float a, float b ) 
{	
	return (a+b)/2.0f;	
}

double avg2d ( double a, double b ) 
{	
	return (a+b)/2.0;	
}

int in_a_set ( int v, int nSet, int *Set )
{
	int i;

	for ( i = 0; i < nSet; i++ ) {
		if ( v == Set[i] )
			return 1;
	}

	return 0;
}

int sign ( float v )
{
	if ( v > 0.0 )
		return 1;
	else if ( v == 0 )
		return 0;
	else
		return -1;
}

void swap2f ( V2DF v1, V2DF v2 )
{
	V2DF tp;
	vcopy2 ( tp, v1 );
	vcopy2 ( v1, v2 );
	vcopy2 ( v2, tp );
}

void swap3f ( V3DF v1, V3DF v2 )
{
	V3DF tp;
	vcopy ( tp, v1 );
	vcopy ( v1, v2 );
	vcopy ( v2, tp );
}

void swap2i ( V2DI v1, V2DI v2 )
{
	V2DI tp;
	vcopy2i ( tp, v1 );
	vcopy2i ( v1, v2 );
	vcopy2i ( v2, tp );
}

void swap ( int *l1, int *l2 )
{
	int t;
	t = *l1;
	*l1 = *l2;
	*l2 = t;
}


void swap ( float *l1, float *l2 )
{
	float t;
	t = *l1;
	*l1 = *l2;
	*l2 = t;
}

void swap ( double *l1, double *l2 )
{
	double t;
	t = *l1;
	*l1 = *l2;
	*l2 = t;
}

void vscale2 ( V2DF v, float div )
{
  v[0] *= div;
  v[1] *= div;
}

void vscale2 ( V2DD v, double div )
{
  v[0] *= div;
  v[1] *= div;
}

void vscale ( V3DF v, float t )
{
  v[0] *= t;
  v[1] *= t;
  v[2] *= t;
}

void vscale ( V3DD tar, double t, V3DD src )
{
  tar[0] = t * src[0];
  tar[1] = t * src[1];
  tar[2] = t * src[2];
}

void vscale ( V3DF tar, float t, V3DF src )
{
  tar[0] = t * src[0];
  tar[1] = t * src[1];
  tar[2] = t * src[2];
}

void vscale ( V3DD v, double t )
{
  v[0] *= t;
  v[1] *= t;
  v[2] *= t;
}

inline void vcopy ( V3DF des, V3DF src )
{
	des[0] = src[0];
	des[1] = src[1];
	des[2] = src[2];
}

// void vzero ( float *v )
// {
//   v[0] = 0.0f;
//   v[1] = 0.0f;
//   v[2] = 0.0f;
// }

// void vset ( float *v, float x, float y, float z )
// {
//   v[0] = x;
//   v[1] = y;
//   v[2] = z;
// }

// void vsub ( float *dst, float *src1, float *src2 )
// {
//   dst[0] = src1[0] - src2[0];
//   dst[1] = src1[1] - src2[1];
//   dst[2] = src1[2] - src2[2];
// }

// void vcopy ( float *v1, float *v2 )
// {
// 	int i;
//   for (i = 0; i < 3; i++)
//       v1[i] = v2[i];
// }

// void vcross ( float *cross, float *v1, float *v2 )
// {
//   float temp[3];

//   temp[0] = (v1[1] * v2[2]) - (v1[2] * v2[1]);
//   temp[1] = (v1[2] * v2[0]) - (v1[0] * v2[2]);
//   temp[2] = (v1[0] * v2[1]) - (v1[1] * v2[0]);
//   vcopy(cross, temp);
// }

// float vlength ( const float *v )
// {
//   return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
// }

// void vscale ( float *v, float div )
// {
//   v[0] *= div;
//   v[1] *= div;
//   v[2] *= div;
// }

// void vnormal ( float *v )
// {
//   vscale(v, 1.0 / vlength(v));
// }

// float vdot ( const float *v1, const float *v2 )
// {
//   return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
// }

// void vadd ( float *dst, float *src1, float *src2 )
// {
//   dst[0] = src1[0] + src2[0];
//   dst[1] = src1[1] + src2[1];
//   dst[2] = src1[2] + src2[2];
// }

void vcopy ( V3DD des, V3DD src )
{
	des[0] = src[0];
	des[1] = src[1];
	des[2] = src[2];
}

void vcopy2 ( V2DF des, V2DF src )
{
	des[0] = src[0];
	des[1] = src[1];
}

void vcopy2 ( V2DD des, V2DD src )
{
	des[0] = src[0];
	des[1] = src[1];
}

void vcopy4 ( V4DF des, V4DF src )
{
	des[0] = src[0];
	des[1] = src[1];
	des[2] = src[2];
	des[3] = src[3];
}

void vcopy3i ( V3DI des, V3DI src )
{
	des[0] = src[0];
	des[1] = src[1];
	des[2] = src[2];
}

void vcopy2i ( V2DI des, V2DI src )
{
	des[0] = src[0];
	des[1] = src[1];
}

void vcopy4i ( V4DI des, V4DI src )
{
	des[0] = src[0];
	des[1] = src[1];
	des[2] = src[2];
	des[3] = src[3];
}

void vminus ( V3DF v, V3DF v1 )
{
	v[0] -= v1[0];
	v[1] -= v1[1];
	v[2] -= v1[2];
}

void vsum ( V3DF v, V3DF v1 )
{
	v[0] += v1[0];
	v[1] += v1[1];
	v[2] += v1[2];
}

inline void vadd ( V3DF v, V3DF v1, V3DF v2 )
{
	v[0] = v1[0] + v2[0];
	v[1] = v1[1] + v2[1];
	v[2] = v1[2] + v2[2];
}

void vadd ( V3DD v, V3DD v1, V3DD v2 )
{
	v[0] = v1[0] + v2[0];
	v[1] = v1[1] + v2[1];
	v[2] = v1[2] + v2[2];
}

inline void vsub ( V3DF v, V3DF v1, V3DF v2 )
{
	v[0] = v1[0] - v2[0];
	v[1] = v1[1] - v2[1];
	v[2] = v1[2] - v2[2];
}

void vsub ( V3DD v, V3DD v1, V3DD v2 )
{
	v[0] = v1[0] - v2[0];
	v[1] = v1[1] - v2[1];
	v[2] = v1[2] - v2[2];
}

void vzero ( V3DI v )
{
	v[0] = 0;
	v[1] = 0;
	v[2] = 0;
}

inline void vzero4i ( V4DI v )
{
	v[0] = 0;
	v[1] = 0;
	v[2] = 0;
	v[3] = 0;
}

inline void vzero ( V3DD v )
{
	v[0] = 0.0;
	v[1] = 0.0;
	v[2] = 0.0;
}

inline void vzero ( V3DF v )
{
	v[0] = 0.0f;
	v[1] = 0.0f;
	v[2] = 0.0f;
}

void vnegate ( V3DF v )
{
	v[0] *= -1.0f;
	v[1] *= -1.0f;
	v[2] *= -1.0f;
}

void vnegate ( V3DD v )
{
	v[0] *= -1.0;
	v[1] *= -1.0;
	v[2] *= -1.0;
}

void vnegative ( V3DF v1, V3DF v2 )
{
	v1[0] = -1.0f*v2[0];
	v1[1] = -1.0f*v2[1];
	v1[2] = -1.0f*v2[2];
}

void vnegative ( V3DD v1, V3DD v2 )
{
	v1[0] = -1.0*v2[0];
	v1[1] = -1.0*v2[1];
	v1[2] = -1.0*v2[2];
}
/*	Normalize given vector v	*/
void vnorm ( V3DF v )
{
	int i;
	float norm = (float) sqrt ( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );

	if ( norm != 0.0f )
		for ( i = 0; i < 3; i++ )
			v[i] = v[i]/norm;

}

void vnorm ( V3DD v )
{
	int i;
	double norm = sqrt ( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );

	if ( norm != 0 )
		for ( i = 0; i < 3; i++ )
			v[i] = v[i]/norm;

}

void vectori ( V2DI v, int p1, int p2 )
{
	v[0] = p1;
	v[1] = p2;
}

void vectori ( V3DI v, int p1, int p2, int p3 )
{
	v[0] = p1;
	v[1] = p2;
	v[2] = p3;
}

void vectori ( V3DI v, int p1, int p2, int p3, int p4, int p5 )
{
	v[0] = p1;
	v[1] = p2;
	v[2] = p3;
	v[3] = p4;
	v[4] = p5;
}

void vectori ( V3DI v, V3DI src, V3DI des )
{
	v[0] = src[0] - des[0];
	v[1] = src[1] - des[1];
	v[2] = src[2] - des[2];
}

void nvectori ( V3DF v, V3DI src, V3DI des )
{
	v[0] = (float) (src[0] - des[0]);
	v[1] = (float) (src[1] - des[1]);
	v[2] = (float) (src[2] - des[2]);

	vnorm ( v );
}

void vscaled_sum3 ( V3DF w, float wt, V3DF v )
{
	w[0] += wt * v[0];
	w[1] += wt * v[1];
	w[2] += wt * v[2];
}

void vscaled_sum4 ( V4DF w, float wt, V4DF v )
{
	w[0] += wt * v[0];
	w[1] += wt * v[1];
	w[2] += wt * v[2];
	w[3] += wt * v[3];
}

float vleng ( V3DF v )
{
	return sqrt ( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
}

float vleng2 (V3DF v )
{
	return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

void vector2 ( V2DI v, int f1, int f2 )
{
	v[0] = f1;
	v[1] = f2;
}

void vector2 ( V2DF v, float f1, float f2 )
{
	v[0] = f1;
	v[1] = f2;
}

void vector2 ( V2DD v, double f1, double f2 )
{
	v[0] = f1;
	v[1] = f2;
}

void vector ( V3DF v, float f1, float f2, float f3 )
{
	v[0] = f1;
	v[1] = f2;
	v[2] = f3;
}

void vector ( V4DF v, float f1, float f2, float f3, float f4 )
{
	v[0] = f1;
	v[1] = f2;
	v[2] = f3;
	v[3] = f4;
}

void vector ( V3DD v, double f1, double f2, double f3 )
{
	v[0] = f1;
	v[1] = f2;
	v[2] = f3;
}

/*	v = p1 - p2		*/
void vector2 ( V2DF v, V2DF p1, V2DF p2 )
{
	v[0] = p1[0] - p2[0];
	v[1] = p1[1] - p2[1];
}

void vector2 ( V2DD v, V2DD p1, V2DD p2 )
{
	v[0] = p1[0] - p2[0];
	v[1] = p1[1] - p2[1];
}

/*	v = p1 - p2		*/
void vector ( V3DF v, V3DF p1, V3DF p2 )
{
	v[0] = p1[0] - p2[0];
	v[1] = p1[1] - p2[1];
	v[2] = p1[2] - p2[2];
}

void vector ( V3DD v, V3DD p1, V3DD p2 )
{
	v[0] = p1[0] - p2[0];
	v[1] = p1[1] - p2[1];
	v[2] = p1[2] - p2[2];
}

void vector4 ( V4DF v, float f1, float f2, float f3, float f4 )
{
	v[0] = f1;
	v[1] = f2;
	v[2] = f3;
	v[3] = f4;
}

void vector4 ( V4DD v, double f1, double f2, double f3, double f4 )
{
	v[0] = f1;
	v[1] = f2;
	v[2] = f3;
	v[3] = f4;
}

void nvector ( V3DF v, float f1, float f2, float f3 )
{
	v[0] = f1;
	v[1] = f2;
	v[2] = f3;

	vnorm ( v );
}

void nvector ( V3DD v, double f1, double f2, double f3 )
{
	v[0] = f1;
	v[1] = f2;
	v[2] = f3;

	vnorm ( v );
}

/*	v = |p1 - p2|		*/
void nvector ( V3DF v, V3DF p1, V3DF p2 )
{
	int i;

	for ( i = 0; i < 3; i++ )
		v[i] = p1[i] - p2[i];

	vnorm ( v );
}

void nvector ( V3DD v, V3DD p1, V3DD p2 )
{
	int i;

	for ( i = 0; i < 3; i++ )
		v[i] = p1[i] - p2[i];

	vnorm ( v );
}

void nvector ( V3DF v, V3DI p1, V3DI p2 )
{
	int i;

	for ( i = 0; i < 3; i++ )
		v[i] = (float) (p1[i] - p2[i]);

	vnorm ( v );
}

void nvector ( V3DD v, V3DF p1, V3DF p2 )
{
	int i;

	for ( i = 0; i < 3; i++ )
		v[i] = p1[i] - p2[i];

	vnorm ( v );
}

void nvector ( V2DD v, double f1, double f2 )
{
	double norm = sqrt(v[0]*v[0] + v[1]*v[1]);

	v[0] = f1/norm;
	v[1] = f2/norm;
}

/*	vc = v1 X v2	*/
void nvcross ( V3DF vc, V3DF v1, V3DF v2 )
{
	vc[0] = v1[1]*v2[2] - v1[2]*v2[1];
	vc[1] = v1[2]*v2[0] - v1[0]*v2[2];
	vc[2] = v1[0]*v2[1] - v1[1]*v2[0];

	vnorm ( vc );
}

void nvcross ( V3DD vc, V3DD v1, V3DD v2 )
{
	vc[0] = v1[1]*v2[2] - v1[2]*v2[1];
	vc[1] = v1[2]*v2[0] - v1[0]*v2[2];
	vc[2] = v1[0]*v2[1] - v1[1]*v2[0];

	vnorm ( vc );
}

void nvcross ( V3DF vc, V3DD v1, V3DF v2 )
{
	vc[0] = (float) v1[1]*v2[2] - (float) v1[2]*v2[1];
	vc[1] = (float) v1[2]*v2[0] - (float) v1[0]*v2[2];
	vc[2] = (float) v1[0]*v2[1] - (float) v1[1]*v2[0];

	vnorm ( vc );
}

void nvcross ( V3DF vc, V3DF v1, V3DD v2 )
{
	vc[0] = v1[1]* (float) v2[2] - v1[2]* (float) v2[1];
	vc[1] = v1[2]* (float) v2[0] - v1[0]* (float) v2[2];
	vc[2] = v1[0]* (float) v2[1] - v1[1]* (float) v2[0];

	vnorm ( vc );
}

/*	vc = v1 X v2	*/
inline void vcross ( V3DF vc, V3DF v1, V3DF v2 )
{
	vc[0] = v1[1]*v2[2] - v1[2]*v2[1];
	vc[1] = v1[2]*v2[0] - v1[0]*v2[2];
	vc[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

void vcross ( cv::Vec3f &vc, V3DF v1, V3DF v2 )
{
	vc[0] = v1[1]*v2[2] - v1[2]*v2[1];
	vc[1] = v1[2]*v2[0] - v1[0]*v2[2];
	vc[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

void vcross ( V3DD vc, V3DD v1, V3DD v2 )
{
	vc[0] = v1[1]*v2[2] - v1[2]*v2[1];
	vc[1] = v1[2]*v2[0] - v1[0]*v2[2];
	vc[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

float vdot4 ( V4DF v1, V4DF v2 )
{
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2] + v1[3]*v2[3];
}

double vdot4 ( V4DD v1, V4DD v2 )
{
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2] + v1[3]*v2[3];
}

float vdot ( V3DF v1, V3DF v2 )
{
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

double vdot ( V3DD v1, V3DD v2 )
{
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

float vdot_scaled ( V3DF v1, V3DF v2 )
{
	if ( vleng (v1) == 0.0f || vleng(v2) == 0.0f )
		return 0.0f;

	return (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])/(vleng(v1) * vleng(v2));
}

void nvavg ( V3DF vt, V3DF v1, V3DF v2 )
{
	int i;

	for ( i = 0; i < 3; i++ )
		vt[i] = (v1[i] + v2[i])/2;

	vnorm ( vt );
}

int match_V4DI ( V4DI v1, V4DI v2 )
{
	int i;

	for ( i = 0; i < 4; i++ )
		if ( v1[i] != v2[i] )
			return 0;

	return 1;
}

//	we assume vt & nm are normalized
void get_projected_v ( V3DF pvt, V3DF vt, V3DF pt, V3DF nm )
{
	float dot = vdot ( vt, nm );
	nvector ( pvt, vt[0] - dot * nm[0], vt[1] - dot * nm[1], vt[2] - dot * nm[2] );
}

/*	ppt is a projected V3DF of pt on the line which is parallel
	to a unit vector v and passes through pt0		*/
float get_projected ( V3DF ppt, V3DF pt, V3DF v, V3DF pt0 )
{
	int i;
	float t;

	for ( i = 0, t = 0.0; i < 3; i++ )
		t -= ( pt0[i] - pt[i] )*v[i];

	for ( i = 0; i < 3; i++ )
		ppt[i] = v[i]*t + pt0[i];

	return t;

}

double get_projected ( V3DD ppt, V3DD pt, V3DD v, V3DD pt0 )
{
	int i;
	double t;

	for ( i = 0, t = 0.0; i < 3; i++ )
		t -= ( pt0[i] - pt[i] )*v[i];

	for ( i = 0; i < 3; i++ )
		ppt[i] = v[i]*t + pt0[i];

	return t;

}

double get_projected_2 ( V3DD pt, V3DD pt0, V3DD pt1 )
{
	return ((pt[0]-pt0[0])*(pt1[0]-pt0[0])+(pt[1]-pt0[1])*(pt1[1]-pt0[1])+(pt[2]-pt0[2])*(pt1[2]-pt0[2]))/((pt1[0]-pt0[0])*(pt1[0]-pt0[0])+(pt1[1]-pt0[1])*(pt1[1]-pt0[1])+(pt1[2]-pt0[2])*(pt1[2]-pt0[2]));
}

/*	The Euclidean distance between pt1 and pt2	*/
float pdist ( V3DI pt1, V3DI pt2 )
{
	int i;
	float tmp;

	for ( i = 0, tmp = 0.0; i < 3; i++ )
		tmp += (pt1[i]-pt2[i])*(pt1[i]-pt2[i]);
	tmp = (float) sqrt ( tmp );

	return tmp;
}

float pdist ( V3DF pt1, V3DF pt2 )
{
	int i;
	float tmp;

	for ( i = 0, tmp = 0.0; i < 3; i++ )
		tmp += (pt1[i]-pt2[i])*(pt1[i]-pt2[i]);

	return (float) sqrt ( tmp );
}

double pdist ( V3DD pt1, V3DD pt2 )
{
	int i;
	double tmp;

	for ( i = 0, tmp = 0.0; i < 3; i++ )
		tmp += (pt1[i]-pt2[i])*(pt1[i]-pt2[i]);

	return sqrt ( tmp );
}


/*	The Euclidean distance between pt1 and pt2	*/
float pdist2 ( V3DF pt1, V3DF pt2 )
{
	int i;
	float tmp;

	for ( i = 0, tmp = 0.0f; i < 3; i++ )
		tmp += (pt1[i]-pt2[i])*(pt1[i]-pt2[i]);

	return tmp;
}

double pdist2 ( V3DD pt1, V3DD pt2 )
{
	int i;
	double tmp;

	for ( i = 0, tmp = 0.0; i < 3; i++ )
		tmp += (pt1[i]-pt2[i])*(pt1[i]-pt2[i]);

	return tmp;
}


/*	The Ellipsoidal distance function	*/
float edist ( V3DF pt1, V3DF pt2, V3DF elp )
{
	int i;
	V3DF x;

	for ( i = 0; i < 3; i++ )
		x[i] = (pt1[i] - pt2[i])/elp[i];

	return (float) sqrt ( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] );
}

float angle_d_2PIF ( V3DF v1, V3DF v2, V3DF nm )
{
	return angle_r_2PIF ( v1, v2, nm ) * (180.0f/PIF);
}

/*	The angle ( degree ) between v1 and v2	*/
//	0 < angle < 180
float angle_d ( V3DF v1, V3DF v2 )
{
	int i;
	float temp, temp2;

	for ( i = 0, temp = 0.0f; i < 3; i++ ) {
		temp += v1[i]*v2[i];
	}
	temp2 = (float) sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2])*(float) sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]);	
	if ( temp2 == 0.0f )
		return 0.0f;
	temp = temp/temp2;

	temp = ( temp <= 1.0f ) ? temp : 1.0f;

	if ( temp >= 0.999999999f )
		return 0.0f;
	if ( temp <= -0.999999999f )
		return 180.0f;

	return (float) acos(temp)*(180.0f/PIF);
}

/*	The angle ( degree ) between (p1 - p2) and (p3 - p2)	*/
//	0 < angle < 180
float angle_d ( V3DF p1, V3DF p2, V3DF p3 )
{
	V3DF v1, v2;

	nvector ( v1, p1, p2 );
	nvector ( v2, p3, p2 );

	return angle_d ( v1, v2 );
}

float angle_d_2PIF ( V3DF p1, V3DF p2, V3DF p3, V3DF nm )
{
	V3DF v1, v2;

	nvector ( v1, p1, p2 );
	nvector ( v2, p3, p2 );

	return angle_d_2PIF ( v1, v2, nm );
}

/*	The angle ( degree ) between v1 and v2	*/
//	0 < angle < 180
double angle_d ( V3DD v1, V3DD v2 )
{
	int i;
	double temp, temp2;

	for ( i = 0, temp = 0.0; i < 3; i++ ) {
		temp += v1[i]*v2[i];
	}
	temp2 = sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2])* sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]);	
	if ( temp2 == 0.0 )
		return 0.0;
	temp = temp/temp2;

	temp = ( temp <= 1.0f ) ? temp : 1.0;
	return acos(temp)*(180.0/PI);
}

/*	The angle ( radian ) between v1 and v2	*/
//	0 < angle < PIF
float angle_r ( V3DF v1, V3DF v2 )
{
	float temp, temp2;

	temp =vdot ( v1, v2 );
	temp2 = (float) sqrt(vdot (v1, v1))*(float) sqrt(vdot (v2, v2));	
	if ( temp2 == 0.0f )
		return 0.0f;
	temp = temp/temp2;

	temp = ( temp <= 1.0f ) ? temp : 1.0f;
	return (float) acos(temp);
}

/*	The angle ( radian ) between v1 and v2	*/
//	0 < angle < 2PIF
//	from v1 to v2
float angle_r_2PIF ( V3DF v1, V3DF v2, V3DF nm )
{
	float temp, temp2;

	temp =vdot ( v1, v2 );
	temp2 = (float) sqrt(vdot (v1, v1))*(float) sqrt(vdot (v2, v2));	
	if ( temp2 == 0.0f )
		return 0.0f;
	temp = temp/temp2;

	if ( temp > 0.0f ) 
		temp = min2f (temp, 1.0f);
	else
		temp = max2f (temp, -1.0f);

	V3DF vnm;
	nvcross ( vnm, v1, v2 );
	float ac = (float) acos(temp);
	if ( vdot ( vnm, nm ) > 0.0f ) {
		return ac;
	}
	else {
		return 2.0f * PIF - ac;
	}
}

//	from v1 to v2
double angle_r_2PI ( V3DD v1, V3DD v2, V3DD nm )
{
	double temp, temp2;

	temp =vdot ( v1, v2 );
	temp2 = sqrt(vdot (v1, v1))* sqrt(vdot (v2, v2));	
	if ( temp2 == 0.0f )
		return 0.0;
	temp = temp/temp2;

	if ( temp > 0 ) 
		temp = min2d (temp, 1.0);
	else
		temp = max2d (temp, -1.0f);

	V3DD vnm;
	nvcross ( vnm, v1, v2 );
	if ( vdot ( vnm, nm ) > 0.0f ) {
		return  acos(temp);
	}
	else {
		return 2.0 * PI - acos(temp);
	}
}

/*	The angle ( radian ) between v1 and v2	*/
//	0 < angle < PIF
double angle_r ( V3DD v1, V3DD v2 )
{
	int i;
	double temp, temp2;

	for ( i = 0, temp = 0.0; i < 3; i++ ) {
		temp += v1[i]*v2[i];
	}
	temp2 = sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2])* sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]);	
	if ( temp2 == 0.0 )
		return 0.0;
	temp = temp/temp2;

	temp = ( temp <= 1.0f ) ? temp : 1.0;
	return acos(temp);
}

/*	The angle ( radian ) between (p1 - p2) and (p3 - p2)	*/
//	0 < angle < PIF
float angle_r ( V3DF p1, V3DF p2, V3DF p3 )
{
	V3DF v1, v2;

	nvector ( v1, p1, p2 );
	nvector ( v2, p3, p2 );

	return angle_r ( v1, v2 );
}

double angle_r ( V3DD p1, V3DD p2, V3DD p3 )
{
	V3DD v1, v2;

	nvector ( v1, p1, p2 );
	nvector ( v2, p3, p2 );

	return angle_r ( v1, v2 );
}


int is_zero_vector ( V3DF v )
{
	return ( (v[0] == 0.0f) && (v[1] == 0.0f) && (v[2] == 0.0f) );
}

int is_zero_vector ( cv::Vec3f &v )
{
	return ( (v[0] == 0.0f) && (v[1] == 0.0f) && (v[2] == 0.0f) );
}

int is_zero_vector ( V3DD v )
{
	return ( (v[0] == 0.0) && (v[1] == 0.0) && (v[2] == 0.0) );
}

int is_same_vector ( V3DF v1, V3DF v2 )
{
	if ( (v1[0] == v2[0]) && (v1[1] == v2[1]) && (v1[2] == v2[2]) )
		return 1;
	else
		return 0;
}

int is_same_vector ( V3DD v1, V3DD v2 )
{
	if ( (v1[0] == v2[0]) && (v1[1] == v2[1]) && (v1[2] == v2[2]) )
		return 1;
	else
		return 0;
}

int is_same_vector ( V3DI v1, V3DI v2 )
{
	if ( (v1[0] == v2[0]) && (v1[1] == v2[1]) && (v1[2] == v2[2]) )
		return 1;
	else
		return 0;
}

int is_same_vector4i ( V4DI v1, V4DI v2 )
{
	if ( (v1[0] == v2[0]) && (v1[1] == v2[1]) && (v1[2] == v2[2]) && (v1[3] == v2[3]) )
		return 1;
	else
		return 0;
}

void vavg3 ( V3DF mpt, V3DF v1, V3DF v2, V3DF v3 )
{
	mpt[0] = (v1[0] + v2[0] + v3[0])/3.0f;
	mpt[1] = (v1[1] + v2[1] + v3[1])/3.0f;
	mpt[2] = (v1[2] + v2[2] + v3[2])/3.0f;
}

void vavg2 ( V2DF mpt, V2DF v1, V2DF v2 )
{
	int i;

	for ( i = 0; i < 2; i++ )
		mpt[i] = (v1[i] + v2[i])/2.0f;
}

void vavg2 ( V2DD mpt, V2DD v1, V2DD v2 )
{
	int i;

	for ( i = 0; i < 2; i++ )
		mpt[i] = (v1[i] + v2[i])/2.0;
}

void vavg ( V3DF mpt, V3DF v1, V3DF v2 )
{
	int i;

	for ( i = 0; i < 3; i++ )
		mpt[i] = (v1[i] + v2[i])/2.0f;
}

void vavg ( V3DD mpt, V3DD v1, V3DD v2 )
{
	int i;

	for ( i = 0; i < 3; i++ )
		mpt[i] = (v1[i] + v2[i])/2.0;
}

//	nml = (p0 - p1)X(p2 - p1)
void get_plane ( V3DF p0, V3DF p1, V3DF p2, V3DF nml )
{
	V3DF v1, v2;

	nvector ( v1, p0, p1 );
	nvector ( v2, p2, p1 );
	nvcross ( nml, v1, v2 );
}


int is_near_zero ( float val, float thr )
{
	if ( fabs ( val ) <= thr )
		return true;
	else
		return false;
}

int is_near_zero ( double val, double thr )
{
	if ( fabs ( val ) <= thr )
		return true;
	else
		return false;
}

int min2i ( int a, int b )
{
	if ( a < b )
		return a;
	else
		return b;
}

int min4i ( int a, int b, int c, int d )
{
	return min2i ( min2i ( a, b ), min2i ( c, d ) );
}

float min2f ( float a, float b )
{
	if ( a < b )
		return a;
	else
		return b;
}

double min2d ( double a, double b )
{
	if ( a < b )
		return a;
	else
		return b;
}

float min2fp ( float a, float b )
{
	if ( a < 0.0f && b < 0.0f )
		return -1.0f;
	if ( a < 0.0f )
		return b;
	if ( b < 0.0f )
		return a;

	if ( a < b )
		return a;
	else
		return b;
}

float min3f ( float a, float b, float c )
{
	return min2f ( min2f(a, b), c );
}


float min4f ( float a, float b, float c, float d )
{
	return min2f ( min2f(a, b), min2f(c, d) );
}

float min5f ( float a, float b, float c, float d, float e )
{
	return min3f ( min2f(a, b), min2f(c, d), e );
}

float max2f ( float a, float b )
{
	if ( a > b )
		return a;
	else
		return b;
}

double max2d ( double a, double b )
{
	if ( a > b )
		return a;
	else
		return b;
}

int max2i ( int a, int b )
{
	if ( a > b )
		return a;
	else
		return b;
}

float max3f ( float a, float b, float c )
{
	return max2f ( max2f(a, b), c );
}

float max4f ( float a, float b, float c, float d )
{
	return max2f ( max2f(a, b), max2f(c, d) );
}

float max5f ( float a, float b, float c, float d, float e )
{
	return max3f ( max2f(a, b), max2f(c, d), e );
}

double max4d ( double a, double b, double c, double d )
{
	return max2d ( max2d(a, b), max2d(c, d) );
}

int max4fi ( float v[4] )
{
	float max = max2f ( max2f(v[0], v[1]), max2f(v[2], v[3]) );
	for ( int i = 0; i < 4; i++ )
		if ( max == v[i] )
			return i;
	return -1;
}

int min4fip ( float v[4] )
{
	float min = min2fp ( min2fp(v[0], v[1]), min2fp(v[2], v[3]) );
	for ( int i = 0; i < 4; i++ )
		if ( min == v[i] )
			return i;
	return -1;
}

int project_round ( float v )
{
	int t;

	t = (int) floor ( v );

	if ( v - t >= 0.5f )
		return t+1;
	else
		return t;
}

int is_similar_vector ( V3DF p1, V3DF p2, float thr )
{
//	if ( is_near_zero ( p1[0]-p2[0], thr ) && is_near_zero ( p1[1]-p2[1], thr ) && is_near_zero ( p1[2]-p2[2], thr ) )
	if ( is_near_zero ( pdist(p1,p2), thr ) )
		return true;
	else
		return false;
}

int is_similar_vector ( V3DD p1, V3DD p2, double thr )
{
//	if ( is_near_zero ( p1[0]-p2[0], thr ) && is_near_zero ( p1[1]-p2[1], thr ) && is_near_zero ( p1[2]-p2[2], thr ) )
	if ( is_near_zero ( pdist(p1,p2), thr ) )
		return true;
	else
		return false;
}

int is_similar_reverse_vector ( V3DF p1, V3DF p2, float thr )
{
	if ( is_near_zero ( p1[0]+p2[0], thr ) && is_near_zero ( p1[1]+p2[1], thr ) && is_near_zero ( p1[2]+p2[2], thr ) )
		return true;
	else
		return false;
}
/*
void get_normal ( V3DF v1, V3DF v2, V3DF v3, V4DF nm )
{
	V3DF nml, p1, p2;

	nvector ( p1, v1, v2 );
	nvector ( p2, v3, v2 );
	nvcross ( nml, p2, p1 );

	nm[0] = nml[0];		nm[1] = nml[1];		nm[2] = nml[2];
	nm[3] = -(nm[0]*v1[0] + nm[1]*v1[1] + nm[2]*v1[2]);

}
*/
void get_normal ( V3DF nm, V3DF p1, V3DF p2, V3DF p3 )
{
	V3DF v1, v2;

	nvector ( v1, p1, p2 );
	nvector ( v2, p3, p2 );
	nvcross ( nm, v2, v1 );
}

void get_normal ( V3DD nm, V3DD p1, V3DD p2, V3DD p3 )
{
	V3DD v1, v2;

	nvector ( v1, p1, p2 );
	nvector ( v2, p3, p2 );
	nvcross ( nm, v2, v1 );
}

float ldist( V3DF pt, V3DF lp, V3DF rp)   
{
	float ab2, bc2, ca2, at ; //d_ab, d_bc, d_ca;
	float cosA, cosB;

	ab2 = pdist(lp , rp);
	ab2 = ab2*ab2;
	bc2 = pdist(rp , pt);
	bc2 = bc2*bc2;
	ca2 = pdist(pt , lp);
	ca2 = ca2*ca2;

	if ( (bc2!=0)&&(ca2!=0) ) {
		cosA = (bc2-ca2-ab2)/(-2.0f*(float)sqrt(ca2)*(float)sqrt(ab2));
		cosB = (ca2-ab2-bc2)/(-2.0f*(float)sqrt(ab2)*(float)sqrt(bc2));
	}
	else  {
		cosA=0; cosB=0;
	}
	
	if ((cosA==1) || (cosB==1))
		return 0;

	else if ((cosA>0) && (cosB>0)) {    // between line and point
		at = (bc2-ca2-ab2)/(-2.0f*(float)sqrt(ab2));
		return  (float) sqrt(ca2-(at*at));
	}

	else if (cosA<0)
		return (float) sqrt(ca2);   // between A point and pt point

	else if (cosB<0)
		return (float) sqrt(bc2);   // between B point and pt point

	else 
		return 0.0f;
}

//	distance from line segment: outside vertices will get 0
//	returns -1, if p is left to lp
//	returns -2, if p is right to rp
float lsdist ( V3DF pt, V3DF lp, V3DF rp)   
{
	float t;
	V3DF v, vt, ppt;
	float dist = pdist ( lp, rp );

	nvector ( v, rp, lp );
	vector ( vt, pt, lp );
	t = vdot ( vt, v );
	if ( t < 0.0f )
		return -1.0f;
	if (t > dist )
		return -2.0f;

	get_point ( ppt, lp, t, v );
	return pdist ( pt, ppt );

}

float ldist2 ( V3DF pt, V3DF lp, V3DF rp )   
{
	float t;
	V3DF v, vt, ppt;

	nvector ( v, rp, lp );
	vector ( vt, pt, lp );
	t = vdot ( vt, v );
//	if ( t < 0.0f || t > pdist(lp,rp) )
//		return 10000.0f;

	get_point ( ppt, lp, t, v );
	return pdist ( pt, ppt );

}

float ldist3 ( V3DF pt, V3DF p0, V3DF vec )   
{
	V3DF v1, v2;

	nvector ( v1, pt, p0 );
	vcopy ( v2, vec );
	vnorm ( v2 );

	float ang = angle_r ( v1, v2 );

	return pdist ( pt, p0 ) * sin (ang);
}

int In_Out_Test( V3DF pt, V3DF p1, V3DF p2, V3DF p3, V3DF pnorm ) 
{	
	int test1, test2, test3;
	V3DF paxab, pbxbc, pcxca;
	V3DF vab, vbc, vca, vac;
	V3DF vpa, vpb, vpc;

	vector ( vab, p2, p1 );
	vector ( vbc, p3, p2 );
	vector ( vca, p1, p3 );

	vector ( vpa, p1, pt );
	vector ( vpb, p2, pt );
	vector ( vpc, p3, pt );

//	vnorm(vector_pa);
//	vnorm(vector_pb);
//	vnorm(vector_pc);

	vector ( vac, p3, p1 );

	nvcross( pnorm, vab, vac );

	vcross( paxab, vpa, vab );
	vcross( pbxbc, vpb, vbc );
	vcross( pcxca, vpc, vca );

	if ( vdot (paxab, pnorm) >= 0 )
		test1 = 1;
	else 
		test1 = -1;

	if ( vdot (pbxbc, pnorm) >= 0 )
		test2 = 1;
	else 
		test2 = -1;

	if ( vdot (pcxca, pnorm) >= 0 )
		test3 = 1;
	else 
		test3 = -1;

	if ( (test1>0) && (test2>0) && (test3>0) )
		return 1;
	else if ( (test1<0) && (test2>0) && (test3<0) )
		return 2;
	else if ( (test1<0) && (test2<0) && (test3>0) )
		return 3;
	else if ( (test1>0) && (test2<0) && (test3<0) )
		return 4;
	else if ( (test1<0) && (test2>0) && (test3>0) )
		return 5;
	else if ( (test1>0) && (test2<0) && (test3>0) )
		return 6;
	else if ( (test1>0) && (test2>0) && (test3<0) )
		return 7;
	else         // all zero -> error !
		return 0;
	
}

/// distance between point and triangle
float tdist(V3DF pt, V3DF p1, V3DF p2, V3DF p3)
{
	V3DF pnorm;	
	float dist, d;  // ax+by+cz+d : a=plane_normal[0]
					//				b=plane_normal[1]
					//				c=plane_normal[2]
					//				d=-A point dot plane_normal

	int result = In_Out_Test( pt, p1, p2, p3, pnorm );

	switch ( result ) {
	case 1 :         //  point is in the triangle
			d = pnorm[0]*(-p1[0])+pnorm[1]*(-p1[1])+pnorm[2]*(-p1[2]);
			dist = (pnorm[0]*pt[0]+pnorm[1]*pt[1]+pnorm[2]*pt[2]+d);				
			if (dist<0)
				dist = -dist;

			// D = | ax+by+cz+d | / sqrt(a^2+b^2+c^2)
			return dist;///sqrt(vdot(plane_normal, plane_normal));
		break;
	case 2 :
		return pdist(pt, p1);
		break;
	case 3 :
		return pdist(pt, p2);
		break;
	case 4 :
		return pdist(pt, p3);
		break;
	case 5 :
		return ldist(pt, p1, p2);
		break;
	case 6 :
		return ldist(pt, p2, p3);
		break;
	case 7 :
		return ldist(pt, p3, p1);
		break;
	default :
		printf("error -- distance between triangle and point\n");
		printf("pt - %f %f %f\n", pt[0], pt[1], pt[2]);
		exit(0);
		return 0.0;
		break;
	}
}

void rotate_vtx_x ( V3DF vtx, float ang )
{
	int i, j;
	float rad = ang*(PIF/180.0f);
	float rotmat[3][3];
	V3DF rvtx;

	for ( i = 0; i < 3; i++ )
		for ( j = 0; j < 3; j++ )
			rotmat[i][j] = 0.0;
	rotmat[1][1] = (float) cos ( rad );
	rotmat[1][2] = (float) -sin ( rad );
	rotmat[2][1] = (float) sin ( rad );
	rotmat[2][2] = (float) cos ( rad );
	rotmat[0][0] = 1.0f;

	for ( i = 0; i < 3; i++ )
		for ( j = 0, rvtx[i] = 0; j < 3; j++ )
			rvtx[i] += rotmat[i][j]*vtx[j];

//	printf("%f, %f, %f: %f, %f, %f\n", vtx[0], vtx[1], vtx[2], rvtx[0], rvtx[1], rvtx[2] );

	vcopy ( vtx, rvtx );
}

void rotate_vtx_y ( V3DF vtx, float ang )
{
	int i, j;
	float rad = ang*(PIF/180.0f);
	float rotmat[3][3];
	V3DF rvtx;

	for ( i = 0; i < 3; i++ )
		for ( j = 0; j < 3; j++ )
			rotmat[i][j] = 0.0;
	rotmat[0][0] = (float) cos ( rad );
	rotmat[0][2] = (float) -sin ( rad );
	rotmat[2][0] = (float) sin ( rad );
	rotmat[2][2] = (float) cos ( rad );
	rotmat[1][1] = 1.0f;

	for ( i = 0; i < 3; i++ )
		for ( j = 0, rvtx[i] = 0; j < 3; j++ )
			rvtx[i] += rotmat[i][j]*vtx[j];

//	printf("%f, %f, %f: %f, %f, %f\n", vtx[0], vtx[1], vtx[2], rvtx[0], rvtx[1], rvtx[2] );

	vcopy ( vtx, rvtx );
}

void rotate_vtx_z ( V3DF vtx, float ang )
{
	int i, j;
	float rad = ang*(PIF/180.0f);
	float rotmat[3][3];
	V3DF rvtx;

	for ( i = 0; i < 3; i++ )
		for ( j = 0; j < 3; j++ )
			rotmat[i][j] = 0.0;
	rotmat[0][0] = (float) cos ( rad );
	rotmat[0][1] = (float) -sin ( rad );
	rotmat[1][0] = (float) sin ( rad );
	rotmat[1][1] = (float) cos ( rad );
	rotmat[2][2] = 1.0f;

	for ( i = 0; i < 3; i++ )
		for ( j = 0, rvtx[i] = 0; j < 3; j++ )
			rvtx[i] += rotmat[i][j]*vtx[j];

//	printf("%f, %f, %f: %f, %f, %f\n", vtx[0], vtx[1], vtx[2], rvtx[0], rvtx[1], rvtx[2] );

	vcopy ( vtx, rvtx );
}

void rotate_vtx ( V3DF vtx/*���*/, V3DF pt, V3DF ang, V3DF org/*�߽�*/ )
{
	V3DF newpt;

	vcopy ( newpt, pt );
	vsub ( newpt, newpt, org );

	if ( ang[0] != 0.0 )
		rotate_vtx_x ( newpt, -ang[0] );
	if ( ang[1] != 0.0 )
		rotate_vtx_y ( newpt, -ang[1] );
	if ( ang[2] != 0.0 )
		rotate_vtx_z ( newpt, -ang[2] );

	vadd ( vtx, newpt, org );
}

void rotate_pt ( V3DF pts, V3DF ang, V3DF org )
{
	V3DF newpt;

	vcopy ( newpt, pts );
	vsub ( newpt, newpt, org );
	if ( ang[0] != 0.0 )
		rotate_vtx_x ( newpt, ang[0] );
	if ( ang[1] != 0.0 )
		rotate_vtx_y ( newpt, ang[1] );
	if ( ang[2] != 0.0 )
		rotate_vtx_z ( newpt, ang[2] );

	vadd ( pts, newpt, org );
}

void rotate_by_matrix ( V3DF rpt, Matrix4DF rot_mat, V3DF pt )
{
	V4DF tpt, npt;

	vector4 ( tpt, pt[0], pt[1], pt[2], 1.0f );
	mult_mat_pt_4 ( npt, rot_mat, tpt );
	vector ( rpt, npt[0], npt[1], npt[2] );

}

void rotate_by_matrix ( V3DD rpt, Matrix4DD rot_mat, V3DD pt )
{
	V4DD tpt, npt;

	vector4 ( tpt, pt[0], pt[1], pt[2], 1.0 );
	mult_mat_pt_4 ( npt, rot_mat, tpt );
	vector ( rpt, npt[0], npt[1], npt[2] );

}

void get_point ( V3DF pt, V3DF p0, float t, V3DF v )
{
	pt[0] = p0[0] + t*v[0];
	pt[1] = p0[1] + t*v[1];
	pt[2] = p0[2] + t*v[2];
}

void get_point ( V3DD pt, V3DD p0, double t, V3DD v )
{
	pt[0] = p0[0] + t*v[0];
	pt[1] = p0[1] + t*v[1];
	pt[2] = p0[2] + t*v[2];
}

void get_point ( V3DF pt, V3DF p0, float t, V3DD v )
{
	pt[0] = p0[0] + t*(float)v[0];
	pt[1] = p0[1] + t*(float)v[1];
	pt[2] = p0[2] + t*(float)v[2];
}

void internal_subdivision_3D ( V3DI wpt, int v1, int v2, V3DI p1, V3DI p2 )
{
	wpt[0] = (v2*p1[0] + v1*p2[0])/(v1+v2);
	wpt[1] = (v2*p1[1] + v1*p2[1])/(v1+v2);
	wpt[2] = (v2*p1[2] + v1*p2[2])/(v1+v2);
}


void internal_subdivision_3D ( V3DF wpt, float v1, float v2, V3DF p1, V3DF p2 )
{
	wpt[0] = (v2*p1[0] + v1*p2[0])/(v1+v2);
	wpt[1] = (v2*p1[1] + v1*p2[1])/(v1+v2);
	wpt[2] = (v2*p1[2] + v1*p2[2])/(v1+v2);
}

void internal_subdivision_2D ( V2DI wpt, int v1, int v2, V2DI p1, V2DI p2 )
{
	wpt[0] = (v2*p1[0] + v1*p2[0])/(v1+v2);
	wpt[1] = (v2*p1[1] + v1*p2[1])/(v1+v2);
}

void internal_subdivision_2D ( V2DF wpt, float v1, float v2, V2DF p1, V2DF p2 )
{
	wpt[0] = (v2*p1[0] + v1*p2[0])/(v1+v2);
	wpt[1] = (v2*p1[1] + v1*p2[1])/(v1+v2);
}

void Matrix ( Matrix2DF m, float m11, float m12, float m21, float m22 )
{
	m[0][0] = m11;		m[0][1] = m12;
	m[1][0] = m21;		m[1][1] = m22;
}

void Matrix ( Matrix2DD m, double m11, double m12, double m21, double m22 )
{
	m[0][0] = m11;		m[0][1] = m12;
	m[1][0] = m21;		m[1][1] = m22;
}

void Matrix ( Matrix3DF m, float m11, float m12, float m13, float m21, float m22, float m23, float m31, float m32, float m33 )
{
	m[0][0] = m11;		m[0][1] = m12;		m[0][2] = m13;
	m[1][0] = m21;		m[1][1] = m22;		m[1][2] = m23;
	m[2][0] = m31;		m[2][1] = m32;		m[2][2] = m33;
}

void Matrix ( Matrix3DD m, double m11, double m12, double m13, double m21, double m22, double m23, double m31, double m32, double m33 )
{
	m[0][0] = m11;		m[0][1] = m12;		m[0][2] = m13;
	m[1][0] = m21;		m[1][1] = m22;		m[1][2] = m23;
	m[2][0] = m31;		m[2][1] = m32;		m[2][2] = m33;
}

void Matrix ( Matrix4DF m, float m11, float m12, float m13, float m14, float m21, float m22, float m23, float m24, float m31, float m32, float m33, float m34, float m41, float m42, float m43, float m44 )
{
	m[0][0] = m11;		m[0][1] = m12;		m[0][2] = m13;		m[0][3] = m14;
	m[1][0] = m21;		m[1][1] = m22;		m[1][2] = m23;		m[1][3] = m24;
	m[2][0] = m31;		m[2][1] = m32;		m[2][2] = m33;		m[2][3] = m34;
	m[3][0] = m41;		m[3][1] = m42;		m[3][2] = m43;		m[3][3] = m44;
}

void Matrix ( Matrix4DD m, double m11, double m12, double m13, double m14, double m21, double m22, double m23, double m24, double m31, double m32, double m33, double m34, double m41, double m42, double m43, double m44 )
{
	m[0][0] = m11;		m[0][1] = m12;		m[0][2] = m13;		m[0][3] = m14;
	m[1][0] = m21;		m[1][1] = m22;		m[1][2] = m23;		m[1][3] = m24;
	m[2][0] = m31;		m[2][1] = m32;		m[2][2] = m33;		m[2][3] = m34;
	m[3][0] = m41;		m[3][1] = m42;		m[3][2] = m43;		m[3][3] = m44;
}

float get_determinant2 ( float m11, float m12, float m21, float m22 )
{
	return m11*m22 - m12*m21;
}

double get_determinant2 ( double m11, double m12, double m21, double m22 )
{
	return m11*m22 - m12*m21;
}

float get_determinant3 ( float m11, float m12, float m13, float m21, float m22, float m23, float m31, float m32, float m33 )
{
	return m11*get_determinant2 ( m22, m23, m32, m33 ) - m12*get_determinant2 ( m21, m23, m31, m33 ) + m13*get_determinant2 ( m21, m22, m31, m32 );
}

double get_determinant3 ( double m11, double m12, double m13, double m21, double m22, double m23, double m31, double m32, double m33 )
{
	return m11*get_determinant2 ( m22, m23, m32, m33 ) - m12*get_determinant2 ( m21, m23, m31, m33 ) + m13*get_determinant2 ( m21, m22, m31, m32 );
}

float get_determinant4 ( float m11, float m12, float m13, float m14, float m21, float m22, float m23, float m24, float m31, float m32, float m33, float m34, float m41, float m42, float m43, float m44 )
{
	return m11*get_determinant3 ( m22, m23, m24, m32, m33, m34, m42, m43, m44 ) - m12*get_determinant3 ( m21, m23, m24, m31, m33, m34, m41, m43, m44 ) + m13*get_determinant3 ( m21, m22, m24, m31, m32, m34, m41, m42, m44 ) - m14*get_determinant3 ( m21, m22, m23, m31, m32, m33, m41, m42, m43 );
}

double get_determinant4 ( double m11, double m12, double m13, double m14, double m21, double m22, double m23, double m24, double m31, double m32, double m33, double m34, double m41, double m42, double m43, double m44 )
{
	return m11*get_determinant3 ( m22, m23, m24, m32, m33, m34, m42, m43, m44 ) - m12*get_determinant3 ( m21, m23, m24, m31, m33, m34, m41, m43, m44 ) + m13*get_determinant3 ( m21, m22, m24, m31, m32, m34, m41, m42, m44 ) - m14*get_determinant3 ( m21, m22, m23, m31, m32, m33, m41, m42, m43 );
}

void get_conj_m ( float **subm, int n, int x, int y, float **m )
{
	int i, j, ri, rj;

	for ( i = 0, ri = 0; i < n; i++ ) {
		if ( i == x )
			continue;
		for ( j = 0, rj = 0; j < n; j++ ) {
			if ( j == y )
				continue;
			subm[ri][rj] = m[i][j];
			rj++;
		}
		ri++;
	}
}

float get_determinant ( int n, float **m )
{
	int i;
	float det;
	float **subm;

	if ( n == 1 )
		return m[0][0];

	if ( n == 2 )
		return m[0][0]*m[1][1] - m[1][0]*m[0][1];

	subm = (float **) calloc ( n-1, sizeof(float *) );
	for ( i = 0; i < n-1; i++ )
		subm[i] = (float *) calloc ( n-1, sizeof(float) );

	for ( i = 0, det = 0.0; i < n; i++ ) {
		get_conj_m ( subm, n, 0, i, m );
		det += (float) pow(-1.0f, i)*m[0][i]*get_determinant ( n-1, subm );
	}

	for ( i = 0; i < n-1; i++ )
		free ( subm[i] );
	free ( subm );
	return det;
}

void get_inverse_4_4 ( Matrix4DF invm, Matrix4DF m )
{
	int i, j, k, l;
	int coi, coj;
	float cofactor, cof[3][3];
	float detM;

	detM = get_determinant4 ( m[0][0], m[0][1], m[0][2], m[0][3], m[1][0], m[1][1], m[1][2], m[1][3], m[2][0], m[2][1], m[2][2], m[2][3], m[3][0], m[3][1], m[3][2], m[3][3] );
	if ( detM == 0 ) {
		for ( i = 0; i < 4; i++ )
			for ( j = 0; j < 4; j++ )
				invm[i][j] = 0.0;
	}
	else {
		for ( i = 0; i < 4; i++ ) {
			for ( j = 0; j < 4; j++ ) {
				for ( k = 0, coi = 0; k < 4; k++ ) {
					for ( l = 0, coj = 0; l < 4; l++ ) {
						if ( k != i && l != j )
							cof[coi][coj++] = m[k][l];
					}
					if ( k != i )
						coi++;
				}
				cofactor = (float) pow(-1.0f, (float)i+(float)j)*get_determinant3 ( cof[0][0], cof[0][1], cof[0][2], cof[1][0], cof[1][1], cof[1][2], cof[2][0], cof[2][1], cof[2][2] );
				invm[j][i] = cofactor/detM;
			}
		}
	}
}

void get_inverse_4_4 ( Matrix4DD invm, Matrix4DD m )
{
	int i, j, k, l;
	int coi, coj;
	double cofactor, cof[3][3];
	double detM;

	detM = get_determinant4 ( m[0][0], m[0][1], m[0][2], m[0][3], m[1][0], m[1][1], m[1][2], m[1][3], m[2][0], m[2][1], m[2][2], m[2][3], m[3][0], m[3][1], m[3][2], m[3][3] );
	if ( detM == 0 ) {
		for ( i = 0; i < 4; i++ )
			for ( j = 0; j < 4; j++ )
				invm[i][j] = 0.0;
	}
	else {
		for ( i = 0; i < 4; i++ ) {
			for ( j = 0; j < 4; j++ ) {
				for ( k = 0, coi = 0; k < 4; k++ ) {
					for ( l = 0, coj = 0; l < 4; l++ ) {
						if ( k != i && l != j )
							cof[coi][coj++] = m[k][l];
					}
					if ( k != i )
						coi++;
				}
				cofactor = pow(-1.0, (double)i+(double)j)*get_determinant3 ( cof[0][0], cof[0][1], cof[0][2], cof[1][0], cof[1][1], cof[1][2], cof[2][0], cof[2][1], cof[2][2] );
				invm[j][i] = cofactor/detM;
			}
		}
	}
}

void get_inverse_3_3 ( Matrix3DF invm, Matrix3DF m )
{
	int i, j, k, l;
	int coi, coj;
	float cofactor, cof[2][2];
	float detM;

	detM = get_determinant3 ( m[0][0], m[0][1], m[0][2], m[1][0], m[1][1], m[1][2], m[2][0], m[2][1], m[2][2] );
	if ( detM == 0 ) {
		for ( i = 0; i < 3; i++ )
			for ( j = 0; j < 3; j++ )
				invm[i][j] = 0.0;
	}
	else {
		for ( i = 0; i < 3; i++ ) {
			for ( j = 0; j < 3; j++ ) {
				for ( k = 0, coi = 0; k < 3; k++ ) {
					for ( l = 0, coj = 0; l < 3; l++ ) {
						if ( k != i && l != j )
							cof[coi][coj++] = m[k][l];
					}
					if ( k != i )
						coi++;
				}
				cofactor = (float) pow(-1.0f, (float)i+(float)j)*get_determinant2 ( cof[0][0], cof[0][1], cof[1][0], cof[1][1] );
				invm[j][i] = cofactor/detM;
			}
		}
	}
}

void get_inverse_3_3 ( Matrix3DD invm, Matrix3DD m )
{
	int i, j, k, l;
	int coi, coj;
	double cofactor, cof[2][2];
	double detM;

	detM = get_determinant3 ( m[0][0], m[0][1], m[0][2], m[1][0], m[1][1], m[1][2], m[2][0], m[2][1], m[2][2] );
	if ( detM == 0 ) {
		for ( i = 0; i < 3; i++ )
			for ( j = 0; j < 3; j++ )
				invm[i][j] = 0.0;
	}
	else {
		for ( i = 0; i < 3; i++ ) {
			for ( j = 0; j < 3; j++ ) {
				for ( k = 0, coi = 0; k < 3; k++ ) {
					for ( l = 0, coj = 0; l < 3; l++ ) {
						if ( k != i && l != j )
							cof[coi][coj++] = m[k][l];
					}
					if ( k != i )
						coi++;
				}
				cofactor = pow(-1.0, (double)i+(double)j)*get_determinant2 ( cof[0][0], cof[0][1], cof[1][0], cof[1][1] );
				invm[j][i] = cofactor/detM;
			}
		}
	}
}

void get_inv_mat ( float **im, int n, float **m )
{
	int i, j;
	float det = get_determinant ( n, m );
	float **adj;
	float **tmp;

	adj = (float **) calloc ( n, sizeof(float *) );
	for ( i = 0; i < n; i++ )
		adj[i] = (float *) calloc ( n, sizeof(float) );

	tmp = (float **) calloc ( n - 1, sizeof(float *) );
	for ( i = 0; i < n - 1; i++ )
		tmp[i] = (float *) calloc ( n - 1, sizeof(float) );

	for ( i = 0; i < n; i++ ) {
		for ( j = 0; j < n; j++ ) {
			get_conj_m ( tmp, n, i, j, m );
			adj[i][j] = ((i+j)%2 == 0 ? 1.0f : -1.0f) * get_determinant ( n-1, tmp );
		}
	}

	for ( i = 0; i < n; i++ ) {
		for ( j = 0; j < n; j++ ) {
			im[i][j] = adj[j][i]/det;
		}
	}

	for ( i = 0; i < n; i++ )
		free ( adj[i] );
	free ( adj );

	for ( i = 0; i < n - 1; i++ )
		free ( tmp[i] );
	free ( tmp );
}

void mult_mat_pt_3 ( V3DF npt, Matrix3DF m, V3DF pt )
{
	int i, j;

	for ( i = 0; i < 3; i++ )
		for ( j = 0, npt[i] = 0.0f; j < 3; j++ )
			npt[i] += m[i][j]*pt[j];
}

void mult_mat_pt_3 ( V3DD npt, Matrix3DD m, V3DD pt )
{
	int i, j;

	for ( i = 0; i < 3; i++ )
		for ( j = 0, npt[i] = 0.0; j < 3; j++ )
			npt[i] += m[i][j]*pt[j];
}

void mult_mat_pt_4 ( V4DF npt, Matrix4DF m, V4DF pt )
{
	int i, j;

	for ( i = 0; i < 4; i++ )
		for ( j = 0, npt[i] = 0.0f; j < 4; j++ )
			npt[i] += m[i][j]*pt[j];
}

void mult_mat_pt_4 ( V4DD npt, Matrix4DD m, V4DD pt )
{
	int i, j;

	for ( i = 0; i < 4; i++ )
		for ( j = 0, npt[i] = 0.0; j < 4; j++ )
			npt[i] += m[i][j]*pt[j];
}

void mult_mat_pt ( V3DF npt, const Matrix4DF m, const V3DF pt )
{
	int i, j;
	V4DF pt4, npt4;

	vector ( pt4, pt[0], pt[1], pt[2], 1.0f );

	for ( i = 0; i < 4; i++ )
		for ( j = 0, npt4[i] = 0.0f; j < 4; j++ )
			npt4[i] += m[i][j]*pt4[j];

	if ( npt4[3] == 0.0f ) {
		printf("Wrong mult_mat_pt\n");
		exit ( 0 );
	}
	vector ( npt, npt4[0]/npt4[3], npt4[1]/npt4[3], npt4[2]/npt4[3] );
}

void mul_mat4 ( Matrix4DF m, Matrix4DF m1, Matrix4DF m2 )
{
	int i, j, k;

	for ( i = 0; i < 4; i++ ) {
		for ( j = 0; j < 4; j++ ) {
			for ( m[i][j] = 0.0f, k = 0; k < 4; k++ )
				m[i][j] += m1[i][k]*m2[k][j];
		}
	}
}

void mul_mat3 ( Matrix3DF m, Matrix3DF m1, Matrix3DF m2 )
{
	int i, j, k;

	for ( i = 0; i < 3; i++ ) {
		for ( j = 0; j < 3; j++ ) {
			for ( m[i][j] = 0.0f, k = 0; k < 3; k++ )
				m[i][j] += m1[i][k]*m2[k][j];
		}
	}
}

void mul_mat2 ( Matrix2DF m, Matrix2DF m1, Matrix2DF m2 )
{
	int i, j, k;

	for ( i = 0; i < 2; i++ ) {
		for ( j = 0; j < 2; j++ ) {
			for ( m[i][j] = 0.0f, k = 0; k < 2; k++ )
				m[i][j] += m1[i][k]*m2[k][j];
		}
	}
}

void multiply_transposed_matrix ( Matrix2DF QMQmat, Matrix2DF Qmat, Matrix2DF Mmat )
{
	Matrix2DF Qt, tmp;

	mmult ( tmp, Mmat, Qmat );
	mtrans ( Qt, Qmat );
	mmult ( QMQmat, Qt, tmp );
}

void multiply_transposed_matrix ( Matrix2DD QMQmat, Matrix2DD Qmat, Matrix2DD Mmat )
{
	Matrix2DD Qt, tmp;

	mmult ( tmp, Mmat, Qmat );
	mtrans ( Qt, Qmat );
	mmult ( QMQmat, Qt, tmp );
}

void multiply_transposed_matrix ( Matrix3DF QMQmat, Matrix3DF Qmat, Matrix3DF Mmat )
{
	Matrix3DF Qt, tmp;

	mmult ( tmp, Mmat, Qmat );
	mtrans ( Qt, Qmat );
	mmult ( QMQmat, Qt, tmp );
}

void multiply_transposed_matrix ( Matrix3DD QMQmat, Matrix3DD Qmat, Matrix3DD Mmat )
{
	Matrix3DD Qt, tmp;

	mmult ( tmp, Mmat, Qmat );
	mtrans ( Qt, Qmat );
	mmult ( QMQmat, Qt, tmp );
}

void multiply_transposed_matrix ( Matrix4DF QMQmat, Matrix4DF Qmat, Matrix4DF Mmat )
{
	Matrix4DF Qt, tmp;

	mmult ( tmp, Mmat, Qmat );
	mtrans ( Qt, Qmat );
	mmult ( QMQmat, Qt, tmp );
}

void multiply_transposed_matrix ( Matrix4DD QMQmat, Matrix4DD Qmat, Matrix4DD Mmat )
{
	Matrix4DD Qt, tmp;

	mmult ( tmp, Mmat, Qmat );
	mtrans ( Qt, Qmat );
	mmult ( QMQmat, Qt, tmp );
}

void copy_set ( int *des, int nsrc, int *src )
{
	int i;

//	des = (int *) calloc ( nsrc, sizeof(int) );
	for ( i = 0; i < nsrc; i++ )
		des[i] = src[i];
}

int in_set ( int a, int nset, int *set )
{
	int i;

	for ( i = 0; i < nset; i++ )
		if ( a == set[i] )
			return 1;

	return 0;
}

//	set1 and set2 are assumed to be sorted
int minus_set2i ( int *mset, int nset1, int nset2, int *set1, int *set2 )
{
	int i;
	int ntset;
	int *tset;

	tset = (int *) calloc ( max2i(nset1, nset2), sizeof(int) );

	for ( i = 0, ntset = 0; i < nset1; i++ ) {
		if ( !in_set ( set1[i], nset2, set2 ) )
			tset[ntset++] = set1[i];
	}
//	mset = (int *) calloc ( ntset, sizeof(int) );
	for ( i = 0; i < ntset; i++ )
		mset[i] = tset[i];

	return ntset;
}

int union_set2i ( int *tset, int num, int nset1, int *set1 )
{
	int i;

//	tset = (int *) calloc ( nset1 + 1, sizeof(int) );

	if ( in_set ( num, nset1, set1 ) )
		return nset1;

	i = 0;
	while ( set1[i] < num ) {
		tset[i] = set1[i];
		i++;
	}
	tset[i++] = num;
	while ( i <= nset1 ) {
		tset[i] = set1[i-1];
		i++;
	}

	return nset1+1;
}

int union_set2i ( int *uset, int nset1, int nset2, int *set1, int *set2 )
{
	int i, iuset, iset1, iset2;
	int nuset = nset1 + nset2;
	int *tset = (int *) calloc ( nuset, sizeof(int) );

	iset1 = iset2 = iuset = 0;
	do {
		if ( set1[iset1] < set2[iset2] ) {
			tset[iuset++] = set1[iset1++];
		}
		else if ( set1[iset1] > set2[iset2] ) {
			tset[iuset++] = set2[iset2++];
		}
		else {
			tset[iuset++] = set1[iset1++];
			iset2++;
		}
	} while ( iset1 < nset1 && iset2 < nset2 );
	if ( iset1 < nset1 ) {
		tset[iuset++] = set1[iset1++];
	}
	if ( iset2 < nset2 ) {
		tset[iuset++] = set2[iset2++];
	}
//	uset = (int *) calloc ( nuset, sizeof(int) );
	for ( i = 0; i < iuset; i++ )
		uset[i] = tset[i];

	return iuset;
}

int merge_set2i ( int n_set1, int n_set2, int *set1, int *set2 )
{
	int i;

	for ( i = 0; i < n_set2; i++ ) {
		if ( in_set ( set2[i], n_set1, set1 ) )
			continue;

		set1[n_set1++] = set2[i];
	}

	return n_set1;
}

int identical_set ( int n, int *set1, int *set2 )
{
	int i;

	for ( i = 0; i < n; i++ )
		if ( set1[i] != set2[i] )
			return 0;

	return 1;
}

int identical_set ( int n1, int n2, int *set1, int *set2 )
{
	if ( n1 != n2 )
		return 0;

	return identical_set ( n1, set1, set2 );
}

int identical_set ( int n, short *set1, short *set2 )
{
	int i;

	for ( i = 0; i < n; i++ )
		if ( set1[i] != set2[i] )
			return 0;

	return 1;
}

int identical_set ( int n1, int n2, short *set1, short *set2 )
{
	if ( n1 != n2 )
		return 0;

	return identical_set ( n1, set1, set2 );
}

//	return 1, if set1 > set2
int contain_seti ( int n1, int n2, int *set1, int *set2 )
{
	int i;

	for ( i = 0; i < n2; i++ ) {
		if ( !in_set ( set2[i], n1, set1 ) )
			return 0;
	}

	return 1;
}

int is_element_of_set ( int x, int n, int *set )
{
	for ( int i = 0; i < n; i++ )
		if ( x == set[i] )
			return 1;

	return 0;
}

void mult_quat ( V4DF quat, V4DF q1, V4DF q2 )
{
	quat[0] = q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3];
	quat[1] = q1[2]*q2[3] - q1[3]*q2[2] + q1[0]*q2[1] + q1[1]*q2[0];
	quat[2] = q1[3]*q2[1] - q1[1]*q2[3] + q1[0]*q2[2] + q1[2]*q2[0];
	quat[3] = q1[1]*q2[2] - q1[2]*q2[1] + q1[0]*q2[3] + q1[3]*q2[0];
}

void build_rot_quat ( V4DF quatq, V3DF pt, V3DF axis, float theta )
{
	int i;

	for ( i = 0, quatq[0] = (float)cos(theta); i < 3; i++ )
		quatq[i+1] = (float)sin(theta)*axis[i];
}

void rotate_point ( V3DF p2, V3DF p1, V3DF pt, V3DF axis, float theta )
{
	int i;
	V3DF np;
	V4DF quatp, quatq, quatqc;
	V4DF tmpq1, tmpq2;

	for ( i = 0, quatp[0] = 0.0; i < 3; i++ )
		quatp[i+1] = p1[i] - pt[i];

	build_rot_quat ( quatq, pt, axis, theta );
	for ( i = 0, quatqc[0] = quatq[0]; i < 3; i++ )
		quatqc[i+1] = -quatq[i+1];
	mult_quat ( tmpq1, quatq, quatp );
	mult_quat ( tmpq2, tmpq1, quatqc );

	for ( i = 0; i < 3; i++ )
		np[i] = tmpq2[i+1] + pt[i];

	vcopy ( p2, np );
}

void mscale ( Matrix2DF m, float s )
{
	for ( int i = 0; i < 2; i++ )
		for ( int j = 0; j < 2; j++ )
			m[i][j] *= s;
}

void mscale ( Matrix2DD m, double s )
{
	for ( int i = 0; i < 2; i++ )
		for ( int j = 0; j < 2; j++ )
			m[i][j] *= s;
}

void mscale ( Matrix3DF m, float s )
{
	for ( int i = 0; i < 3; i++ )
		for ( int j = 0; j < 3; j++ )
			m[i][j] *= s;
}

void mscale ( Matrix3DD m, double s )
{
	for ( int i = 0; i < 3; i++ )
		for ( int j = 0; j < 3; j++ )
			m[i][j] *= s;
}

void mscale ( Matrix4DF m, float s )
{
	for ( int i = 0; i < 4; i++ )
		for ( int j = 0; j < 4; j++ )
			m[i][j] *= s;
}

void mscale ( Matrix4DD m, double s )
{
	for ( int i = 0; i < 4; i++ )
		for ( int j = 0; j < 4; j++ )
			m[i][j] *= s;
}

void madd ( Matrix2DF m, Matrix2DF m1, Matrix2DF m2 )
{
	for ( int i = 0; i < 2; i++ )
		for ( int j = 0; j < 2; j++ )
			m[i][j] = m1[i][j] + m2[i][j];
}

void madd ( Matrix2DD m, Matrix2DD m1, Matrix2DD m2 )
{
	for ( int i = 0; i < 2; i++ )
		for ( int j = 0; j < 2; j++ )
			m[i][j] = m1[i][j] + m2[i][j];
}

void madd ( Matrix3DF m, Matrix3DF m1, Matrix3DF m2 )
{
	for ( int i = 0; i < 3; i++ )
		for ( int j = 0; j < 3; j++ )
			m[i][j] = m1[i][j] + m2[i][j];
}

void madd ( Matrix3DD m, Matrix3DD m1, Matrix3DD m2 )
{
	for ( int i = 0; i < 3; i++ )
		for ( int j = 0; j < 3; j++ )
			m[i][j] = m1[i][j] + m2[i][j];
}

void madd ( Matrix4DF m, Matrix4DF m1, Matrix4DF m2 )
{
	for ( int i = 0; i < 4; i++ )
		for ( int j = 0; j < 4; j++ )
			m[i][j] = m1[i][j] + m2[i][j];
}


void madd ( Matrix4DD m, Matrix4DD m1, Matrix4DD m2 )
{
	for ( int i = 0; i < 4; i++ )
		for ( int j = 0; j < 4; j++ )
			m[i][j] = m1[i][j] + m2[i][j];
}

void msub ( Matrix2DF m, Matrix2DF m1, Matrix2DF m2 )
{
	for ( int i = 0; i < 2; i++ )
		for ( int j = 0; j < 2; j++ )
			m[i][j] = m1[i][j] - m2[i][j];
}

void msub ( Matrix2DD m, Matrix2DD m1, Matrix2DD m2 )
{
	for ( int i = 0; i < 2; i++ )
		for ( int j = 0; j < 2; j++ )
			m[i][j] = m1[i][j] - m2[i][j];
}

void msub ( Matrix3DF m, Matrix3DF m1, Matrix3DF m2 )
{
	for ( int i = 0; i < 3; i++ )
		for ( int j = 0; j < 3; j++ )
			m[i][j] = m1[i][j] - m2[i][j];
}

void msub ( Matrix3DD m, Matrix3DD m1, Matrix3DD m2 )
{
	for ( int i = 0; i < 3; i++ )
		for ( int j = 0; j < 3; j++ )
			m[i][j] = m1[i][j] - m2[i][j];
}

void msub ( Matrix4DF m, Matrix4DF m1, Matrix4DF m2 )
{
	for ( int i = 0; i < 4; i++ )
		for ( int j = 0; j < 4; j++ )
			m[i][j] = m1[i][j] - m2[i][j];
}

void msub ( Matrix4DD m, Matrix4DD m1, Matrix4DD m2 )
{
	for ( int i = 0; i < 4; i++ )
		for ( int j = 0; j < 4; j++ )
			m[i][j] = m1[i][j] - m2[i][j];
}

void mmult ( Matrix2DF m, Matrix2DF m1, Matrix2DF m2 )
{
	int i, j, k;

	for ( i = 0; i < 2; i++ )
		for ( j = 0; j < 2; j++ ) 
			for ( m[i][j] = 0.0, k = 0; k < 2; k++ )
				m[i][j] += m1[i][k]*m2[k][j];
}

void mmult ( Matrix2DD m, Matrix2DD m1, Matrix2DD m2 )
{
	int i, j, k;

	for ( i = 0; i < 2; i++ )
		for ( j = 0; j < 2; j++ ) 
			for ( m[i][j] = 0.0, k = 0; k < 2; k++ )
				m[i][j] += m1[i][k]*m2[k][j];
}

void mmult ( Matrix3DF m, Matrix3DF m1, Matrix3DF m2 )
{
	int i, j, k;

	for ( i = 0; i < 3; i++ )
		for ( j = 0; j < 3; j++ ) 
			for ( m[i][j] = 0.0, k = 0; k < 3; k++ )
				m[i][j] += m1[i][k]*m2[k][j];
}

void mmult ( Matrix3DD m, Matrix3DD m1, Matrix3DD m2 )
{
	int i, j, k;

	for ( i = 0; i < 3; i++ )
		for ( j = 0; j < 3; j++ ) 
			for ( m[i][j] = 0.0, k = 0; k < 3; k++ )
				m[i][j] += m1[i][k]*m2[k][j];
}

void mmult ( Matrix4DF m, Matrix4DF m1, Matrix4DF m2 )
{
	int i, j, k;

	for ( i = 0; i < 4; i++ )
		for ( j = 0; j < 4; j++ ) 
			for ( m[i][j] = 0.0, k = 0; k < 4; k++ )
				m[i][j] += m1[i][k]*m2[k][j];
}

void mmult ( Matrix4DD m, Matrix4DD m1, Matrix4DD m2 )
{
	int i, j, k;

	for ( i = 0; i < 4; i++ )
		for ( j = 0; j < 4; j++ ) 
			for ( m[i][j] = 0.0, k = 0; k < 4; k++ )
				m[i][j] += m1[i][k]*m2[k][j];
}

void mtrans ( Matrix2DF mt, Matrix2DF m )
{
	for ( int i = 0; i < 2; i++ )
		for ( int j = 0; j < 2; j++ )
			mt[i][j] = m[j][i];
}

void mtrans ( Matrix2DD mt, Matrix2DD m )
{
	for ( int i = 0; i < 2; i++ )
		for ( int j = 0; j < 2; j++ )
			mt[i][j] = m[j][i];
}

void mtrans ( Matrix3DF mt, Matrix3DF m )
{
	for ( int i = 0; i < 3; i++ )
		for ( int j = 0; j < 3; j++ )
			mt[i][j] = m[j][i];
}

void mtrans ( Matrix3DD mt, Matrix3DD m )
{
	for ( int i = 0; i < 3; i++ )
		for ( int j = 0; j < 3; j++ )
			mt[i][j] = m[j][i];
}

void mtrans ( Matrix4DF mt, Matrix4DF m )
{
	for ( int i = 0; i < 4; i++ )
		for ( int j = 0; j < 4; j++ )
			mt[i][j] = m[j][i];
}

void mtrans ( Matrix4DD mt, Matrix4DD m )
{
	for ( int i = 0; i < 4; i++ )
		for ( int j = 0; j < 4; j++ )
			mt[i][j] = m[j][i];
}

void load_zero ( Matrix2DF mat )
{
	for ( int i = 0; i < 2; i++ )
		for ( int j = 0; j < 2; j++ )
			mat[i][j] = 0.0f;
}

void load_zero ( Matrix2DD mat )
{
	for ( int i = 0; i < 2; i++ )
		for ( int j = 0; j < 2; j++ )
			mat[i][j] = 0.0;
}

void load_zero ( Matrix3DF mat )
{
	for ( int i = 0; i < 3; i++ )
		for ( int j = 0; j < 3; j++ )
			mat[i][j] = 0.0f;
}

void load_zero ( Matrix3DD mat )
{
	for ( int i = 0; i < 3; i++ )
		for ( int j = 0; j < 3; j++ )
			mat[i][j] = 0.0;
}

void load_zero ( Matrix4DF mat )
{
	for ( int i = 0; i < 4; i++ )
		for ( int j = 0; j < 4; j++ )
			mat[i][j] = 0.0f;
}

void load_zero ( Matrix4DD mat )
{
	for ( int i = 0; i < 4; i++ )
		for ( int j = 0; j < 4; j++ )
			mat[i][j] = 0.0;
}

void load_identity ( Matrix2DF mat )
{
	int i, j;

	for ( i = 0; i < 2; i++ )
		for ( j = 0; j < 2; j++ ) 
			mat[i][j] = 0.0f;

	for ( i = 0; i < 2; i++ ) 
		mat[i][i] = 1.0f;
}

void load_identity ( Matrix2DD mat )
{
	int i, j;

	for ( i = 0; i < 2; i++ )
		for ( j = 0; j < 2; j++ ) 
			mat[i][j] = 0.0;

	for ( i = 0; i < 2; i++ ) 
		mat[i][i] = 1.0;
}

void load_identity ( Matrix3DF mat )
{
	int i, j;

	for ( i = 0; i < 3; i++ )
		for ( j = 0; j < 3; j++ ) 
			mat[i][j] = 0.0f;

	for ( i = 0; i < 3; i++ ) 
		mat[i][i] = 1.0f;
}

void load_identity ( Matrix3DD mat )
{
	int i, j;

	for ( i = 0; i < 3; i++ )
		for ( j = 0; j < 3; j++ ) 
			mat[i][j] = 0.0;

	for ( i = 0; i < 3; i++ ) 
		mat[i][i] = 1.0;
}

void load_identity ( Matrix4DF mat )
{
	int i, j;

	for ( i = 0; i < 4; i++ )
		for ( j = 0; j < 4; j++ ) 
			mat[i][j] = 0.0f;

	for ( i = 0; i < 4; i++ ) 
		mat[i][i] = 1.0f;
}

void load_identity ( Matrix4DD mat )
{
	int i, j;

	for ( i = 0; i < 4; i++ )
		for ( j = 0; j < 4; j++ ) 
			mat[i][j] = 0.0;

	for ( i = 0; i < 4; i++ ) 
		mat[i][i] = 1.0;
}

void get_matrix ( Matrix3DF mat, V3DF v )
{
	for ( int i = 0; i < 3; i++ )
		for ( int j = 0; j < 3; j++ )
			mat[i][j] = v[i]*v[j];
}

void get_matrix ( Matrix3DD mat, V3DD v )
{
	for ( int i = 0; i < 3; i++ )
		for ( int j = 0; j < 3; j++ )
			mat[i][j] = v[i]*v[j];
}

void get_matrix ( Matrix4DF mat, V4DF v )
{
	for ( int i = 0; i < 4; i++ )
		for ( int j = 0; j < 4; j++ )
			mat[i][j] = v[i]*v[j];
}

void get_matrix ( Matrix4DD mat, V4DD v )
{
	for ( int i = 0; i < 4; i++ )
		for ( int j = 0; j < 4; j++ )
			mat[i][j] = v[i]*v[j];
}

void translatemat ( Matrix4DF mat, float p0, float p1, float p2 )
{
	int i, j;

	for ( i = 0; i < 4; i++ ) 
		for ( j = 0; j < 4; j++ ) 
			mat[i][j] = 0.0f;

	for ( i = 0; i < 4; i++ ) 
		mat[i][i] = 1.0f;

	mat[0][3] = p0;
	mat[1][3] = p1;
	mat[2][3] = p2;
}

void translatemat ( Matrix4DD mat, double p0, double p1, double p2 )
{
	int i, j;

	for ( i = 0; i < 4; i++ ) 
		for ( j = 0; j < 4; j++ ) 
			mat[i][j] = 0.0;

	for ( i = 0; i < 4; i++ ) 
		mat[i][i] = 1.0;

	mat[0][3] = p0;
	mat[1][3] = p1;
	mat[2][3] = p2;
}

void rotatemat ( Matrix4DF mat, float ang, V3DF axis )
{
	int i, j;
    float s = (float) -sin ( ang );
    float c = (float) cos ( ang );

    for ( i = 0; i < 4; i++ )
        for ( j = 0; j < 4; j++ )
            mat[i][j] = 0.0f;

    for ( i = 0; i < 4; i++ )
        mat[i][i] = 1.0f;

    mat[0][0] = (1 - c)*axis[0]*axis[0] + c;
    mat[0][1] = (1 - c)*axis[0]*axis[1] + s*axis[2];
    mat[0][2] = (1 - c)*axis[0]*axis[2] - s*axis[1];
    mat[1][0] = (1 - c)*axis[0]*axis[1] - s*axis[2];
    mat[1][1] = (1 - c)*axis[1]*axis[1] + c;
    mat[1][2] = (1 - c)*axis[1]*axis[2] + s*axis[0];
    mat[2][0] = (1 - c)*axis[0]*axis[2] + s*axis[1];
    mat[2][1] = (1 - c)*axis[1]*axis[2] - s*axis[0];
    mat[2][2] = (1 - c)*axis[2]*axis[2] + c;
}

void rotatemat ( Matrix4DD mat, double ang, V3DD axis )
{
	int i, j;
    double s = -sin ( ang );
    double c = cos ( ang );

    for ( i = 0; i < 4; i++ )
        for ( j = 0; j < 4; j++ )
            mat[i][j] = 0.0;

    for ( i = 0; i < 4; i++ )
        mat[i][i] = 1.0;

    mat[0][0] = (1 - c)*axis[0]*axis[0] + c;
    mat[0][1] = (1 - c)*axis[0]*axis[1] + s*axis[2];
    mat[0][2] = (1 - c)*axis[0]*axis[2] - s*axis[1];
    mat[1][0] = (1 - c)*axis[0]*axis[1] - s*axis[2];
    mat[1][1] = (1 - c)*axis[1]*axis[1] + c;
    mat[1][2] = (1 - c)*axis[1]*axis[2] + s*axis[0];
    mat[2][0] = (1 - c)*axis[0]*axis[2] + s*axis[1];
    mat[2][1] = (1 - c)*axis[1]*axis[2] - s*axis[0];
    mat[2][2] = (1 - c)*axis[2]*axis[2] + c;
}

/* mat1 = mat1 * mat2	*/
void mul_mat ( Matrix4DF mat1, Matrix4DF mat2 )
{
	int i, j, k;
	Matrix4DF mat3;

	for ( i = 0; i < 4; i++ )
		for ( j = 0; j < 4; j++ ) 
			for ( mat3[i][j] = 0.0, k = 0; k < 4; k++ )
				mat3[i][j] += mat1[i][k]*mat2[k][j];

	for ( i = 0; i < 4; i++ )
		for ( j = 0; j < 4; j++ ) 
			mat1[i][j] = mat3[i][j];
}

void mul_mat ( Matrix4DD mat1, Matrix4DD mat2 )
{
	int i, j, k;
	Matrix4DD mat3;

	for ( i = 0; i < 4; i++ )
		for ( j = 0; j < 4; j++ ) 
			for ( mat3[i][j] = 0.0, k = 0; k < 4; k++ )
				mat3[i][j] += mat1[i][k]*mat2[k][j];

	for ( i = 0; i < 4; i++ )
		for ( j = 0; j < 4; j++ ) 
			mat1[i][j] = mat3[i][j];
}

/*	make a rotation matrix that rotates at p around axis a by phi	*/
void get_rot_mat ( Matrix4DF mat, V3DF a, V3DF p, float phi )
{
	Matrix4DF temp;

	load_identity ( mat );
	translatemat ( temp, p[0], p[1], p[2] );
	mul_mat ( mat, temp );
	rotatemat ( temp, phi, a );
	mul_mat ( mat, temp );
	translatemat ( temp, -p[0], -p[1], -p[2] );
	mul_mat ( mat, temp );
}

void get_rot_mat ( Matrix4DD mat, V3DD a, V3DD p, double phi )
{
	Matrix4DD temp;

	load_identity ( mat );
	translatemat ( temp, p[0], p[1], p[2] );
	mul_mat ( mat, temp );
	rotatemat ( temp, phi, a );
	mul_mat ( mat, temp );
	translatemat ( temp, -p[0], -p[1], -p[2] );
	mul_mat ( mat, temp );
}

void get_rot_mat ( Matrix4DF mat, float src1, float src2, float src3, float dst1, float dst2, float dst3, float p1, float p2, float p3 )
{
	V3DF src, dst;
	V3DF axis, origin;
	float angle;

	vector ( origin, p1, p2, p3 );
	nvector ( src, src1, src2, src3 );
	nvector ( dst, dst1, dst2, dst3 );
	if ( is_same_vector ( src, dst ) ) {
		angle = 0.0;
		vector ( axis, 1.0, 0.0, 0.0 );
	}
	else {
		angle = (float) acos( vdot( src, dst ) );
		nvcross ( axis, src, dst );
	}
	get_rot_mat ( mat, axis, origin, angle );	
}

void get_rot_mat ( Matrix4DF mat, V3DF src, V3DF dst )
{
	V3DF axis, origin;
	float angle;

	vzero ( origin );
	if ( is_same_vector ( src, dst ) ) {
		angle = 0.0;
		vector ( axis, 1.0, 0.0, 0.0 );
	}
	else {
		angle = (float) acos( vdot( src, dst ) );
		nvcross ( axis, src, dst );
	}
	get_rot_mat ( mat, axis, origin, angle );	
}

void matrix3to4 ( Matrix4DF mat4, Matrix3DF mat3 )
{
	int i, j;

    for ( i = 0; i < 4; i++ )
        for ( j = 0; j < 4; j++ )
            mat4[i][j] = 0.0;

    for ( i = 0; i < 4; i++ )
        mat4[i][i] = 1.0;

	for ( i = 0; i < 3; i++ )
		for ( j = 0; j < 3; j++ )
			mat4[i][j] = mat3[i][j];
}

void get_rot_mat ( Matrix4DF mat, Matrix3DF mat3, V3DF p )
{
	Matrix4DF temp;

	load_identity ( mat );
	translatemat ( temp, p[0], p[1], p[2] );
	mul_mat ( mat, temp );
	matrix3to4 ( temp, mat3 );
	mul_mat ( mat, temp );
	translatemat ( temp, -p[0], -p[1], -p[2] );
	mul_mat ( mat, temp );
}

void get_normals ( V3DF *TrNmls, int nTrs, int nTrPts, V3DI *Trs, V3DF *TrPts )
{
	int i;
	V3DF v1, v2, nml;

	for ( i = 0; i < nTrPts; i++ )
		vzero ( TrNmls[i] );

	for ( i = 0; i < nTrs; i++ ) {
		//	Compute the normals of each faces
		nvector ( v1, TrPts[Trs[i][1]], TrPts[Trs[i][0]] );
		nvector ( v2, TrPts[Trs[i][2]], TrPts[Trs[i][0]] );
		nvcross ( nml, v1, v2 );
		//	Add the normals of the faces to each nml
		vadd ( TrNmls[Trs[i][0]], TrNmls[Trs[i][0]], nml );
		vadd ( TrNmls[Trs[i][1]], TrNmls[Trs[i][1]], nml );
		vadd ( TrNmls[Trs[i][2]], TrNmls[Trs[i][2]], nml );
	}

//	Normalize nml
	for ( i = 0; i < nTrPts; i++ )
		vnorm ( TrNmls[i] );
}

int get_a_projected_point_on_a_plane ( V3DF ppt, V3DF pt, V3DF p1, V3DF nml )
{
	float t;
	V3DF v;

	vector ( v, p1, pt );
	t = vdot ( v, nml )/vdot(nml,nml);
	if ( t >= 0.0f ) {
		get_point ( ppt, pt, t, nml );
		return -1;
//		printf("Wrong projected point\n");
//		exit ( 0 );
	}
	get_point ( ppt, pt, t, nml );
	return 1;
}

int get_a_projected_point_on_a_plane ( V3DD ppt, V3DD pt, V3DD p1, V3DD nml )
{
	double t;
	V3DD v;

	vector ( v, p1, pt );
	t = vdot ( v, nml )/vdot(nml,nml);
	if ( t >= 0.0 ) {
		get_point ( ppt, pt, t, nml );
		return -1;
//		printf("Wrong projected point\n");
//		exit ( 0 );
	}
	get_point ( ppt, pt, t, nml );
	return 1;
}

void get_point_from_barycentric ( V3DF pt, V3DF bary, V3DF p1, V3DF p2, V3DF p3 )
{
	vector ( pt, bary[0]*p1[0]+bary[1]*p2[0]+bary[2]*p3[0],  bary[0]*p1[1]+bary[1]*p2[1]+bary[2]*p3[1],  bary[0]*p1[2]+bary[1]*p2[2]+bary[2]*p3[2] );
}

void get_barycentric ( V3DF bary, V3DF pt, V3DF p1, V3DF p2, V3DF p3 )
{
	float l, s, t;
	V3DF v, ppt;
	V3DF v1, v2, npt;
	nvector ( v, pt, p1 );
//	l = ((p2[0]-p1[0])*(p3[1]-p2[1])-(p2[1]-p1[1])*(p3[0]-p2[0]))/(v[0]*(p3[1]-p2[1])-v[1]*(p3[0]-p2[0]));
//	printf("%f, %f, %f\n",  ((p2[0]-p1[0])*(p3[1]-p2[1])-(p2[1]-p1[1])*(p3[0]-p2[0]))/(v[0]*(p3[1]-p2[1])-v[1]*(p3[0]-p2[0])), 
//							((p2[1]-p1[1])*(p3[2]-p2[2])-(p2[2]-p1[2])*(p3[1]-p2[1]))/(v[1]*(p3[2]-p2[2])-v[2]*(p3[1]-p2[1])), 
//							((p2[2]-p1[2])*(p3[0]-p2[0])-(p2[0]-p1[0])*(p3[2]-p2[2]))/(v[2]*(p3[0]-p2[0])-v[0]*(p3[2]-p2[2])) );

	l = get_determinant2 ( p2[0]-p1[0], p2[0]-p3[0], p2[1]-p1[1], p2[1]-p3[1] )/get_determinant2 ( v[0], p2[0]-p3[0], v[1], p2[1]-p3[1] );
	get_point ( ppt, p1, l, v );
/*
	float l1, k1, l2, k2, l3, k3;
	l1 = get_determinant2 ( p2[0]-p1[0], p2[0]-p3[0], p2[1]-p1[1], p2[1]-p3[1] )/get_determinant2 ( v[0], p2[0]-p3[0], v[1], p2[1]-p3[1] );
	k1 = get_determinant2 ( v[0], p2[0]-p1[0], v[1], p2[1]-p1[1] )/get_determinant2 ( v[0], p2[0]-p3[0], v[1], p2[1]-p3[1] );
	printf("%f, %f > %f, %f, %f\n", l1, k1, l1*v[0]+k1*(p2[0]-p3[0]), l1*v[1]+k1*(p2[1]-p3[1]), l1*v[2]+k1*(p2[2]-p3[2]) );

	l2 = get_determinant2 ( p2[0]-p1[0], p2[0]-p3[0], p2[2]-p1[2], p2[2]-p3[2] )/get_determinant2 ( v[0], p2[0]-p3[0], v[2], p2[2]-p3[2] );
	k2 = get_determinant2 ( v[0], p2[0]-p1[0], v[2], p2[2]-p1[2] )/get_determinant2 ( v[0], p2[0]-p3[0], v[2], p2[2]-p3[2] );
	printf("%f, %f > %f, %f, %f\n", l1, k1, l1*v[0]+k1*(p2[0]-p3[0]), l1*v[1]+k1*(p2[1]-p3[1]), l1*v[2]+k1*(p2[2]-p3[2]) );

	l3 = get_determinant2 ( p2[1]-p1[1], p2[1]-p3[1], p2[2]-p1[2], p2[2]-p3[2] )/get_determinant2 ( v[1], p2[1]-p3[1], v[2], p2[2]-p3[2] );
	k3 = get_determinant2 ( v[1], p2[1]-p1[1], v[2], p2[2]-p1[2] )/get_determinant2 ( v[1], p2[1]-p3[1], v[2], p2[2]-p3[2] );
	printf("%f, %f > %f, %f, %f\n", l1, k1, l1*v[0]+k1*(p2[0]-p3[0]), l1*v[1]+k1*(p2[1]-p3[1]), l1*v[2]+k1*(p2[2]-p3[2]) );

	get_point ( ppt, p1, l, v );
	nvector ( v1, p2, ppt );
	nvector ( v2, p3, ppt );
	printf("%11.10f, %11.10f, %11.10f\n", v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]);

	get_point ( ppt, p1, l1, v );
	nvector ( v1, p2, ppt );
	nvector ( v2, p3, ppt );
	printf("%11.10f, %11.10f, %11.10f\n", v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]);

	get_point ( ppt, p1, l2, v );
	nvector ( v1, p2, ppt );
	nvector ( v2, p3, ppt );
	printf("%11.10f, %11.10f, %11.10f\n", v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]);

	get_point ( ppt, p1, l3, v );
	nvector ( v1, p2, ppt );
	nvector ( v2, p3, ppt );
	printf("%11.10f, %11.10f, %11.10f\n", v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]);
*/

	t = pdist(pt, p1)/l;	
	internal_subdivision_3D ( npt, t, 1.0f-t, p1, ppt );
//	printf("[%f]%f, %f, %f -> ", t, npt[0] - pt[0], npt[1] - pt[1], npt[2] - pt[2]);

	nvector ( v1, ppt, p2 );
	nvector ( v2, p3, p2 );
	if ( vdot(v1,v2) > 0.0f )
		s = pdist(ppt, p2)/pdist(p2,p3);
	else
		s = -pdist(ppt, p2)/pdist(p2,p3);
	internal_subdivision_3D ( npt, s, 1.0f-s, p2, p3 );
//	printf("[%f]%f, %f, %f\n", s, npt[0] - ppt[0], npt[1] - ppt[1], npt[2] - ppt[2]);
	internal_subdivision_3D ( ppt, t, 1.0f-t, p1, npt );
//	printf("%f, %f, %f > ", ppt[0], ppt[1], ppt[2]);
	vector ( bary, 1.0f-t, t*(1.0f-s), s*t );
//	printf("%f, %f, %f\n", bary[0]*p1[0]+bary[1]*p2[0]+bary[2]*p3[0], bary[0]*p1[1]+bary[1]*p2[1]+bary[2]*p3[1], bary[0]*p1[2]+bary[1]*p2[2]+bary[2]*p3[2]);
}

int inside_triangle ( V3DF bary, V3DF pt, V3DF p1, V3DF p2, V3DF p3 )
{
	get_barycentric ( bary, pt, p1, p2, p3 );
//	printf("%f,%f,%f\n", bary[0], bary[1], bary[2]);
//	printf("[%f,%f,%f][%f, %f, %f]:[%f, %f, %f]\n", bary[0], bary[1], bary[2], pt[0], pt[1], pt[2], bary[0]*p1[0]+bary[1]*p2[0]+bary[2]*p3[0], bary[0]*p1[1]+bary[1]*p2[1]+bary[2]*p3[1], bary[0]*p1[2]+bary[1]*p2[2]+bary[2]*p3[2]);
	return ( bary[0] >= 0.0f && bary[1] >= 0.0f && bary[2] >= 0.0f );
}

int is_a_projected_point_on_a_directed_triangle ( V3DF bary, V3DF pt, V3DF p1, V3DF p2, V3DF p3 )
{
	V3DF vec, v1, v2, v3,  nml;
	V3DF cent, ppt;

	vzero ( bary );
//	compute normal of the triangle
	nvector ( v1, p2, p1 );
	nvector ( v2, p3, p1 );
	nvector ( v3, p3, p2 );
	nvcross ( nml, v1, v2 );

//	on the back halfspace
	vector ( cent, (p1[0]+p2[0]+p3[0])/3.0f, (p1[1]+p2[1]+p3[1])/3.0f, (p1[2]+p2[2]+p3[2])/3.0f );
	nvector ( vec, cent, pt );
	if ( vdot ( vec, nml ) >= 0.0f ) {
//		printf("Back:[%f,%f,%f][%f,%f,%f]%f\n", vec[0], vec[1], vec[2], nml[0], nml[1], nml[2], vdot(vec, nml));
		return 0;
	}
//	printf("front\n");
//	return 1;
	//	get a projected point on the plane
	int imsi = get_a_projected_point_on_a_plane ( ppt, pt, p1, nml );

	//	if the point is inside the triangle
//	int imsi = inside_triangle ( bary, ppt, p1, p2, p3 );
//	return 1;

	if ( inside_triangle ( bary, ppt, p1, p2, p3 ) )
		return 1;
	else
		return 0;

}

void transpose ( Matrix4DF mat )
{
	int i, j;
	Matrix4DF tmpmat;

	for ( i = 0; i < 4; i++ )
		for ( j = 0; j < 4; j++ )
			tmpmat[i][j] = mat[i][j];

	for ( i = 0; i < 4; i++ )
		for ( j = 0; j < 4; j++ )
			mat[i][j] = tmpmat[j][i];
}

void vector_rotate ( V3DF tvec, V3DF svec, V3DF axis, V3DF cent, float ang )
{
	V3DF vpt;
	V4DF tpt4, rpt4;
	Matrix4DF mat4;

	vadd ( vpt, svec, cent );
	get_rot_mat ( mat4, axis, cent, ang );
	vector4 ( tpt4, vpt[0], vpt[1], vpt[2], 1.0f );
	mult_mat_pt_4 ( rpt4, mat4, tpt4 );
	vector ( tvec, rpt4[0] - cent[0], rpt4[1] - cent[1], rpt4[2] - cent[2] );
}

void vector_rotate ( V3DD tvec, V3DD svec, V3DD axis, V3DD cent, double ang )
{
	V3DF vpt;
	V4DF tpt4, rpt4;
	Matrix4DF mat4;

	V3DF tvecf, svecf, axisf, centf;
	float angf;

	vector ( tvecf, (float) tvec[0], (float) tvec[1], (float) tvec[2] );
	vector ( svecf, (float) svec[0], (float) svec[1], (float) svec[2] );
	vector ( axisf, (float) axis[0], (float) axis[1], (float) axis[2] );
	vector ( centf, (float) cent[0], (float) cent[1], (float) cent[2] );
	angf = (float) ang;

	vadd ( vpt, svecf, centf );
	get_rot_mat ( mat4, axisf, centf, angf );
	vector4 ( tpt4, vpt[0], vpt[1], vpt[2], 1.0f );
	mult_mat_pt_4 ( rpt4, mat4, tpt4 );
	vector ( tvecf, rpt4[0] - centf[0], rpt4[1] - centf[1], rpt4[2] - centf[2] );

	vector ( tvec, (double) tvecf[0], (double) tvecf[1], (double) tvecf[2] );
}

void vector_rotate ( V3DF tvec, V3DF svec, Matrix4DF mat )
{
	V4DF tv, sv;

	vector4 ( sv, svec[0], svec[1], svec[2], 1.0f );

	mult_mat_pt_4 ( tv, mat, sv );

	vector ( tvec, tv[0]/tv[3], tv[1]/tv[3], tv[2]/tv[3] );
}

void vector_rotate ( V3DF tvec, V3DF svec, V3DF src, V3DF dst )
{
	Matrix4DF rmat;

	build_rot_mat ( rmat, src, dst );
	vector_rotate ( tvec, svec, rmat );
}

void mzero3 ( Matrix3DF m )
{
	m[0][0] = m[0][1] = m[0][2] = 0.0f;
	m[1][0] = m[1][1] = m[1][2] = 0.0f;
	m[2][0] = m[2][1] = m[2][2] = 0.0f;
}

void mzero4 ( Matrix4DF m )
{
	m[0][0] = m[0][1] = m[0][2] = m[0][3] = 0.0f;
	m[1][0] = m[1][1] = m[1][2] = m[1][3] = 0.0f;
	m[2][0] = m[2][1] = m[2][2] = m[2][3] = 0.0f;
	m[3][0] = m[3][1] = m[3][2] = m[3][3] = 0.0f;
}

void compute_intersection_xy ( V3DF pos, V3DF p1, V3DF p2, V3DF q1, V3DF q2 )
{
	float a, b, c, d, v, w, s, t;

	a = p2[0] - p1[0];
	b = q1[0] - q2[0];
	c = p2[1] - p1[1];
	d = q1[1] - q2[1];
	v = q1[0] - p1[0];
	w = q1[1] - p1[1];

	char error_message[256] = "Wrong intersection computation of xy lines\n";

	error_exit ( (a*d-b*c == 0.0f), error_message);

	t = (d*v-b*w)/(a*d-b*c);
	s = (a*w-c*v)/(a*d-b*c);

	V3DF c1, c2;
	vector ( c1, p1[0] + t*(p2[0]-p1[0]), p1[1] + t*(p2[1]-p1[1]), p1[2] + t*(p2[2]-p1[2]) );
	vector ( c2, q1[0] + s*(q2[0]-q1[0]), q1[1] + s*(q2[1]-q1[1]), q1[2] + s*(q2[2]-q1[2]) );
//	if ( !is_same_vector ( c1, c2 ) ) {
//		printf("[[%15.14f, %15.14f, %15.14f][%15.14f, %15.14f, %15.14f]]\n", c1[0], c1[1], c1[2], c2[0], c2[1], c2[2]);
//	}
//	error_exit ( !is_same_vector ( c1, c2 ), "2nd Wrong intersection computation of xy line\n");
	vcopy ( pos, c1 );
}
/*
void compute_intersection_xy ( V3DD pos, V3DD p1, V3DD p2, V3DD q1, V3DD q2 )
{
	double a, b, c, d, v, w, s, t;

	a = p2[0] - p1[0];
	b = q1[0] - q2[0];
	c = p2[1] - p1[1];
	d = q1[1] - q2[1];
	v = q1[0] - p1[0];
	w = q1[1] - p1[1];

	if ( a*d-b*c == 0 )
		printf("%f, %f, %f, %f\n", a, b, c, d);
	error_exit ( (a*d-b*c == 0.0), "Wrong intersection computation of xy lines\n");

//	t = (d*v-b*w)/(a*d-b*c);
//	s = (a*w-c*v)/(a*d-b*c);
	t = ((q1[1] - q2[1])*(q1[0] - p1[0])-(q1[0] - q2[0])*(q1[1] - p1[1]))/((p2[0] - p1[0])*(q1[1] - q2[1])-(q1[0] - q2[0])*(p2[1] - p1[1]));
	s = ((p2[0] - p1[0])*(q1[1] - p1[1])-(p2[1] - p1[1])*(q1[0] - p1[0]))/((p2[0] - p1[0])*(q1[1] - q2[1])-(q1[0] - q2[0])*(p2[1] - p1[1]));
	V3DD c1, c2;
	vector ( c1, p1[0] + t*(p2[0]-p1[0]), p1[1] + t*(p2[1]-p1[1]), p1[2] + t*(p2[2]-p1[2]) );
	vector ( c2, q1[0] + s*(q2[0]-q1[0]), q1[1] + s*(q2[1]-q1[1]), q1[2] + s*(q2[2]-q1[2]) );
	if ( !is_same_vector ( c1, c2 ) ) {
		printf("[%25.24f, %25.24f, %25.24f]\n[%25.24f, %25.24f, %25.24f]\n", c1[0], c1[1], c1[2], c2[0], c2[1], c2[2]);
	}
	error_exit ( !is_same_vector ( c1, c2 ), "2nd Wrong intersection computation of xy line\n");
	vcopy ( pos, c1 );
}
*/
void compute_intersection_xy ( V3DD pos, V3DD p1, V3DD p2, V3DD q1, V3DD q2 )
{
	char message[256] = "Wrong intersection computation of xy lines\n";
	double a, b, c, d, v, w, s, t;

	a = p2[0] - p1[0];
	b = q1[0] - q2[0];
	c = p2[1] - p1[1];
	d = q1[1] - q2[1];
	v = q1[0] - p1[0];
	w = q1[1] - p1[1];

	if ( a*d-b*c == 0 )
		printf("%f, %f, %f, %f\n", a, b, c, d);
	error_exit ( (a*d-b*c == 0.0), message);

//	t = (d*v-b*w)/(a*d-b*c);
//	s = (a*w-c*v)/(a*d-b*c);
	t = ((q1[1] - q2[1])*(q1[0] - p1[0])-(q1[0] - q2[0])*(q1[1] - p1[1]))/((p2[0] - p1[0])*(q1[1] - q2[1])-(q1[0] - q2[0])*(p2[1] - p1[1]));
	s = ((p2[0] - p1[0])*(q1[1] - p1[1])-(p2[1] - p1[1])*(q1[0] - p1[0]))/((p2[0] - p1[0])*(q1[1] - q2[1])-(q1[0] - q2[0])*(p2[1] - p1[1]));

	if ( s < 0.0 || fabs(s) <= 0.00001 )
		s = 0.001;
	if ( s > 1.0 || fabs(s-1.0) <= 0.00001 )
		s = 0.999;
	if ( q1[0] == q2[0] ) {
		vector ( pos, q1[0], q1[1] + s*(q2[1]-q1[1]), q1[2] );
	}
	else if ( q1[1] == q2[1] ) {
		vector ( pos, q1[0] + s*(q2[0]-q1[0]), q1[1], q1[2] );
	}
	else {
			printf("Wrong edge for computing intersection\n");
	}
}

void error_exit ( int value, char *str )
{
	if ( value ) {
		printf("%s", str);
		exit ( 0 );
	}
}

void error_report ( int value, char *str )
{
	if ( value ) 
		printf("%s", str);
}

//	return 1, if the line segment intersects the circle only at one point
//	return 0, if the line segment intersects the circle at zero or two points
int intersection_circle_line_segment ( V3DF ip, V3DF p1, V3DF p2, V3DF p, float radius, int dir )
{
//	if ( (pdist(p1, p) - radius)*(pdist(p2,p) - radius) > 0.0 )
//		return 0;

	float a, b, c;
	a = pdist2 (p1, p2);
	b = (p2[0]-p1[0])*(p1[0]-p[0])+(p2[1]-p1[1])*(p1[1]-p[1])+(p2[2]-p1[2])*(p1[2]-p[2]);
	c = pdist2(p1, p) - radius*radius;

	if (b*b - a*c < 0)
		return 0;

//	float t = (-b - sqrt(b*b - a*c))/a;
	float t = (dir > 0 ) ? (-b + (float) sqrt(b*b - a*c))/a : (-b - (float) sqrt(b*b - a*c))/a;
	V3DF vec;
	vector ( vec, p2, p1 );
	get_point ( ip, p1, t, vec);
	return 1;
}

void build_change_of_basis ( Matrix3DD TSmat, Matrix3DD Smat, Matrix3DD Tmat )
{
	Matrix3DD SImat;

	get_inverse_3_3 ( SImat, Smat );
	mmult ( TSmat, Tmat, SImat );
}

float area_triangle ( V3DF p1, V3DF p2, V3DF p3 )
{
	V3DF v1, v2, v;

	vector ( v1, p1, p2 );
	vector ( v2, p3, p2 );
	vcross ( v, v1, v2 );

	return (float) sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

double area_triangle ( V3DD p1, V3DD p2, V3DD p3 )
{
	V3DD v1, v2, v;

	vector ( v1, p1, p2 );
	vector ( v2, p3, p2 );
	vcross ( v, v1, v2 );

	return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

//	return cosine of the angle between (p3 - p2) & (p1 - p2)
float get_cosine ( V3DF p1, V3DF p2, V3DF p3 )
{
	V3DF v1, v2;

	nvector ( v1, p3, p2 );
	nvector ( v2, p1, p2 );

	return vdot ( v1, v2 );
}

//	returns 1, if above
//			0, if on
//		   -1, if below
int is_a_point_above_a_plane ( V3DF pt, V3DF p, V3DF n )
{
	float angle;
	V3DF v;

	nvector ( v, pt, p );
	angle = vdot ( v, n );
	if ( angle > 0.0 )
		return 1;
	else if ( angle == 0.0f )
		return 0;
	else
		return -1;
}

//	returns 0, if no intersection
//	returns 1, if intersection
int get_an_intersection_bw_edge_N_plane ( V3DF ip, V3DF p1, V3DF p2, V3DF p, V3DF n )
{
	float t;
	V3DF v1, v2;

	vector ( v1, p, p1 );
	vector ( v2, p2, p1 );

	t = vdot(n, v1)/vdot(n, v2);

	get_point ( ip, p1, t, v2 );

	return ( 0.0f <= t && t <= 1.0f );
}

int get_an_intersection_bw_edge_N_plane_pure ( V3DF ip, V3DF p1, V3DF p2, V3DF p, V3DF n )
{
	float t;
	V3DF v1, v2;

	vector ( v1, p, p1 );
	vector ( v2, p2, p1 );

	t = vdot(n, v1)/vdot(n, v2);

	get_point ( ip, p1, t, v2 );

	return ( 0.0f < t && t < 1.0f );
}

//	line: pt & vec
//	line segment: p1 & p2
//	ip = pt + t * vec = p1 + s * (p2 - p1)
//	t * vec - s * (p2 - p1) = p1 - pt
//  (t * vec - s * (p2 - p1))^2 = (p1 - pt)^2
//	(vec DOT vec) t^2 - 2 vec DOT (p2 - p1) t s + (p2 - p1) DOT (p2 - p1) s^2 - (p1 - pt)^2 = 0
//		a <= (vec DOT vec)
//		b <= -2 vec DOT (p2 - p1)
//		c <= (p2 - p1) DOT (p2 - p1)
//		d <= -(p1 - pt) DOT (p1 - pt)
//	a t^2 + b s t + c s^2 + d = 0
//	a (t^2 + b/a st + c/a s^2) + d = 0
//	a (t^2 + b/a st + (b/2a)^2 s^2 + (c/a - (b/2a)^2 s^2) + d = 0
//	a ( t + (b/2a) s)^2 + a(c/a - (b/2a)^2) s^2 + d = 0
//
//	s^2 = -d * (a(c/a - (b/2a)^2))^(-1),
//
//	invalid, if -d * (a(c/a - (b/2a)^2))^(-1) < 0 or 
//				a == 0 or 
//				c/a - (b/2a)^2 == 0
//
//	s = +- sqrt (-d * (a(c/a - (b/2a)^2))^(-1))
//	t = -(b/2a) s
//
int get_an_intersection_line_linesegment ( V3DF ip, V3DF pt, V3DF vec, V3DF p1, V3DF p2 )
{
	float a, b, c, d;
	float t, s, s2;
	V3DF vec1, vec2;

	vector ( vec1, p2, p1 );
	vector ( vec2, p1, pt );

	if ( is_zero_vector ( vec ) )
		return -1;

	a = vdot ( vec, vec );
	b = -2.0f * vdot ( vec, vec1 );
	c = vdot ( vec1, vec1 );
	d = - vdot ( vec2, vec2 );

	s2 = -d / ( a * (c/a - (b / (2*a))*(b / (2*a))) );
	if ( s2 < 0.0f )
		return -1;

	s = sqrt ( s2 );
	t = - ( b/ (2*a) ) * s;
	if ( 0.0f <= s && s <= 1.0f ) {
		get_point ( ip, pt, t, vec );
		return 1;
	}
	
	return -1;

}

//	line: pt + t vec
//	plane: p1 & p2 are points on the plane
//		   vec1 is a vector on the plane
//
//	step 1. find a normal vector of the plane: nm = (p2 - p1) X vec1
//	step 2. find an intersection point of the line & the plane
//	step 3. project the intersection point onto the line segment (p1 & p2)
int get_an_intersection_line_linesegment_nonplanar ( V3DF ip, V3DF p0, V3DF vec, V3DF p1, V3DF p2, V3DF vec1 )
{
	V3DF vec2, nm;

//	step 1.
	nvector ( vec2, p2, p1 );
	nvcross ( nm, vec1, vec2 );

//	step 2. 
//		pt on the line = p0 + t vec
//		(pt - p1) DOT nm = 0
//		(t vec + p0 - p1 ) DOT nm = 0
//		t vec DOT nm = (p1 - p0) DOT nm
//		t = (p1 - p0) DOT nm / vec DOT nm
	if ( vdot ( vec, nm ) == 0.0f )
		return -1;

	V3DF pt;
	float t;
	vector ( vec2, p1, p0 );
	t = vdot ( vec2, nm ) / vdot ( vec, nm );
	get_point ( pt, p0, t, vec );

//	step 3.
//		ppt = p1 + t (p2 - p1)
//		(ppt - pt) DOT (p2 - p1) = 0
//		t (p2 - p1) DOT (p2 - p1) + (p1 - pt) DOT (p2 - p1) = 0
//		t = (pt - p1) DOT (p2 - p1) / (p2 - p1) DOT (p2 - p1)
	V3DF vec3;
	vector ( vec2, p2, p1 );
	vector ( vec3, pt, p1 );
	t = vdot ( vec3, vec2 ) / vdot ( vec2, vec2 );
	get_point ( ip, p1, t, vec2 );
	if ( 0.0f <= t && t <= 1.0f )
		return 1;

	return -1;
}

int between ( float l, float v, float r )
{
	return ( l <= v && v <= r );
}

int a_cube_contains_a_point ( V3DF p1, V3DF p7, V3DF q )
{
	return ( between ( p1[0], q[0], p7[0] ) && between ( p1[1], q[1], p7[1] ) && between ( p1[2], q[2], p7[2] ) );
}

int is_cube_p_contains_cube_q ( V3DF p1, V3DF p2, V3DF p3, V3DF p4, V3DF p5, V3DF p6, V3DF p7, V3DF p8, 
						   V3DF q1, V3DF q2, V3DF q3, V3DF q4, V3DF q5, V3DF q6, V3DF q7, V3DF q8 )
{
	return ( a_cube_contains_a_point ( p1, p7, q1 ) && 
			 a_cube_contains_a_point ( p1, p7, q2 ) && 
			 a_cube_contains_a_point ( p1, p7, q3 ) && 
			 a_cube_contains_a_point ( p1, p7, q4 ) && 
			 a_cube_contains_a_point ( p1, p7, q5 ) && 
			 a_cube_contains_a_point ( p1, p7, q6 ) && 
			 a_cube_contains_a_point ( p1, p7, q7 ) && 
			 a_cube_contains_a_point ( p1, p7, q8 ) );
}

//	p1: Xmin, Ymin, Zmin
//	p2: Xmin, Ymin, Zmax
//	p3: Xmax, Ymin, Zmax
//	p4: Xmax, Ymin, Zmin
//	p5: Xmin, Ymax, Zmin
//	p6: Xmin, Ymax, Zmax
//	p7: Xmax, Ymax, Zmax
//	p8: Xmax, Ymax, Zmin
//	simple intersection test
int cube_p_intersects_cube_q ( V3DF p1, V3DF p2, V3DF p3, V3DF p4, V3DF p5, V3DF p6, V3DF p7, V3DF p8, 
								V3DF q1, V3DF q2, V3DF q3, V3DF q4, V3DF q5, V3DF q6, V3DF q7, V3DF q8 )
{
	return ( a_cube_contains_a_point ( p1, p7, q1 ) ||
			 a_cube_contains_a_point ( p1, p7, q2 ) || 
			 a_cube_contains_a_point ( p1, p7, q3 ) || 
			 a_cube_contains_a_point ( p1, p7, q4 ) || 
			 a_cube_contains_a_point ( p1, p7, q5 ) || 
			 a_cube_contains_a_point ( p1, p7, q6 ) || 
			 a_cube_contains_a_point ( p1, p7, q7 ) || 
			 a_cube_contains_a_point ( p1, p7, q8 ) );
}

int an_edge_intersects_a_triangle ( V3DF p1, V3DF p2, V3DF t1, V3DF t2, V3DF t3 )
{
	V3DF nml;
	V3DF ip, bary;

	get_normal ( nml, t1, t2, t3 );
	if ( get_an_intersection_bw_edge_N_plane ( ip, p1, p2, t1, nml ) )
		return inside_triangle ( bary, ip, t1, t2, t3 );

	return 0;

}

int an_edge_intersects_a_triangle ( V3DF ip, V3DF p1, V3DF p2, V3DF t1, V3DF t2, V3DF t3 )
{
	V3DF nml;
	V3DF bary;

	get_normal ( nml, t1, t2, t3 );
	if ( get_an_intersection_bw_edge_N_plane ( ip, p1, p2, t1, nml ) )
		return inside_triangle ( bary, ip, t1, t2, t3 );

	return 0;

}

int X_ray_intersects_a_triangle ( V3DF ip, V3DF p1, V3DF p2, V3DF t1, V3DF t2, V3DF t3 )
{
	V3DF nml;
	V3DF bary;
	V3DF Max, Min;

//	1. Get bounding box of triangle
	for ( int i = 0; i < 3; i++ ) {
		Max[i] = max3f ( t1[i], t2[i], t3[i] );
		Min[i] = min3f ( t1[i], t2[i], t3[i] );
	}
	if ( between ( Min[1], p1[1], Max[1] ) && between ( Min[2], p1[2], Max[2] ) ) {
		get_normal ( nml, t1, t2, t3 );
		if ( get_an_intersection_bw_edge_N_plane ( ip, p1, p2, t1, nml ) )
			return inside_triangle ( bary, ip, t1, t2, t3 );

		return 0;
	}
	else
		return 0;

}

int Y_ray_intersects_a_triangle ( V3DF ip, V3DF p1, V3DF p2, V3DF t1, V3DF t2, V3DF t3 )
{
	V3DF nml;
	V3DF bary;
	V3DF Max, Min;

//	1. Get bounding box of triangle
	for ( int i = 0; i < 3; i++ ) {
		Max[i] = max3f ( t1[i], t2[i], t3[i] );
		Min[i] = min3f ( t1[i], t2[i], t3[i] );
	}
	if ( between ( Min[0], p1[0], Max[0] ) && between ( Min[2], p1[2], Max[2] ) ) {
		get_normal ( nml, t1, t2, t3 );
		if ( get_an_intersection_bw_edge_N_plane ( ip, p1, p2, t1, nml ) )
			return inside_triangle ( bary, ip, t1, t2, t3 );

		return 0;
	}
	else
		return 0;

}

int Z_ray_intersects_a_triangle ( V3DF ip, V3DF p1, V3DF p2, V3DF t1, V3DF t2, V3DF t3 )
{
	V3DF nml;
	V3DF bary;
	V3DF Max, Min;

//	1. Get bounding box of triangle
	for ( int i = 0; i < 3; i++ ) {
		Max[i] = max3f ( t1[i], t2[i], t3[i] );
		Min[i] = min3f ( t1[i], t2[i], t3[i] );
	}
	if ( between ( Min[0], p1[0], Max[0] ) && between ( Min[1], p1[1], Max[1] ) ) {
		get_normal ( nml, t1, t2, t3 );
		if ( get_an_intersection_bw_edge_N_plane ( ip, p1, p2, t1, nml ) )
			return inside_triangle ( bary, ip, t1, t2, t3 );

		return 0;
	}
	else
		return 0;

}

int does_edge_intersect_plane ( V3DF p1, V3DF p2, V3DF n, V3DF p )
{
	V3DF v1, v2;

	nvector ( v1, p1, p );
	nvector ( v2, p2, p );

	return ( vdot(v1, n) * vdot(v2, n) < 0.0f );
}

//	return 0, if different set
//	return 1, if same elements & same order
//	return -1,if same elements & reversed order
int is_same_triangle (  V3DI t1, V3DI t2 )
{
	int i;
	int cnt1, cnt2;

	for ( i = 0, cnt1 = 0; i < 3; i++ )
		cnt1 += in_set ( t1[i], 3, t2 );

	for ( i = 0, cnt2 = 0; i < 3; i++ )
		cnt2 += in_set ( t2[i], 3, t1 );
	
//	same set
	if ( cnt1 == 3 && cnt2 == 3 ) {
		int idx;
		for ( i = 0; i < 3; i++ ) {
			if ( t2[i] == t1[0] ) {
				idx = i;
				break;
			}
		}
		if ( t1[1] == t2[(idx+1)%3] )
			return 1;
		else
			return -1;
	}
	else
		return 0;
}


int does_triangle_include_edge ( V3DI f, V2DI e )
{
	if ( f[0] == e[0] && f[1] == e[1] || f[0] == e[1] && f[1] == e[0] )
		return 1;

	if ( f[1] == e[0] && f[2] == e[1] || f[1] == e[1] && f[2] == e[0] )
		return 1;

	if ( f[2] == e[0] && f[0] == e[1] || f[2] == e[1] && f[0] == e[0] )
		return 1;

	return 0;
}

int does_triangle_include_edge ( V3DI f, int v1, int v2 )
{
	if ( f[0] == v2 && f[1] == v1 || f[0] == v1 && f[1] == v2 )
		return 1;

	if ( f[1] == v2 && f[2] == v1 || f[1] == v1 && f[2] == v2 )
		return 1;

	if ( f[2] == v2 && f[0] == v1 || f[2] == v1 && f[0] == v2 )
		return 1;

	return 0;
}

int Union_set ( int *uset, int n_set1, int n_set2, int *set1, int *set2 )
{
	int i;
	int n_uset;

	n_uset = n_set1;
	for ( i = 0; i < n_set1; i++ )
		uset[i] = set1[i];

	for ( i = 0; i < n_set2; i++ ) {
		if ( !in_a_set ( set2[i], n_set1, set1 ) )
			uset[n_uset++] = set2[i];
	}

	return n_uset;
}

int Union_set ( int *uset, int n_set1, int n_set2, int *set2 )
{
	int i;
	int n_uset;

	n_uset = n_set1;

	for ( i = 0; i < n_set2; i++ ) {
		if ( !in_a_set ( set2[i], n_set1, uset ) )
			uset[n_uset++] = set2[i];
	}

	return n_uset;
}

void interpolate_point ( V3DF p, V3DF p1, V3DF p2, int a, int b )
{
	int i;

	if ( a+b == 0 )
		printf("Cannot compute interpolation for a & b: %d, %d\n", a, b);

	for ( i = 0; i < 3; i++ )
		p[i] = (b*p1[i] + a*p2[i])/(a+b);
}

void interpolate_point ( V3DF p, V3DF p1, V3DF p2, float a, float b )
{
	int i;

	if ( a+b == 0.0f )
		printf("Cannot compute interpolation for a & b: %f, %f\n", a, b);

	for ( i = 0; i < 3; i++ )
		p[i] = (b*p1[i] + a*p2[i])/(a+b);
}

//	vt = irat[0] * va + irat[1] * vt
void get_interpolation_ratio ( V2DF irat, V3DF va, V3DF vb, V3DF vt )
{
	irat[0] = pdist ( vb, vt ) / pdist ( va, vb );
	irat[1] = pdist ( vt, va ) / pdist ( va, vb );
}
//	returns -1, if no solution such as 0 x^2 + 0 x + c = 0
//	return   0, if imaginary solution such as b^2 - 4 a c < 0
//	return   1, if real solutions
//	sol[0] < sol[1]
int quadratic_polynomial ( V2DF sol, float a, float b, float c )
{
	if ( a == 0.0f ) {
		if ( b == 0.0f )
			return -1;
		else
			sol[0] = -1.0f * c / b;
			sol[1] = -1.0f * c / b;

			return 1;
	}

	if ( b * b - 4.0f * a * c < 0 ) 
		return 0;
	else if ( b * b - 4.0f * a * c == 0 ) {
		sol[0] = (-1.0f * b)/(2.0f * a);
		sol[1] = (-1.0f * b)/(2.0f * a);

		return 1;
	}

	sol[0] = (-1.0f * b + sqrt ( b * b - 4.0f * a * c ))/(2.0f * a);
	sol[1] = (-1.0f * b - sqrt ( b * b - 4.0f * a * c ))/(2.0f * a);

	if ( sol[1] < sol[0] )
		swap ( &(sol[0]), &(sol[1]) );
//	printf("%f, %f\n", a*sol[0]*sol[0] + b*sol[0] + c, a*sol[1]*sol[1] + b*sol[1] + c );
	return 1;
}

int quadratic_polynomial ( V2DD sol, double a, double b, double c )
{
	if ( a == 0.0f ) {
		if ( b == 0.0 )
			return -1;
		else
			sol[0] = -1.0 * c / b;
			sol[1] = -1.0 * c / b;

			return 1;
	}

	if ( b * b - 4.0 * a * c < 0 ) 
		return 0;
	else if ( b * b - 4.0 * a * c == 0 ) {
		sol[0] = (-1.0 * b)/(2.0f * a);
		sol[1] = (-1.0 * b)/(2.0f * a);

		return 1;
	}

	sol[0] = (-1.0 * b + sqrt ( b * b - 4.0 * a * c ))/(2.0 * a);
	sol[1] = (-1.0 * b - sqrt ( b * b - 4.0 * a * c ))/(2.0 * a);

	if ( sol[1] < sol[0] )
		swap ( &(sol[0]), &(sol[1]) );
//	printf("%f, %f\n", a*sol[0]*sol[0] + b*sol[0] + c, a*sol[1]*sol[1] + b*sol[1] + c );
	return 1;
}

void hermite2D ( V2DF pt, float t, V2DF p1, V2DF p2, V2DF v1, V2DF v2 )
{
/*
	Matrix4DF Hmat;
	Matrix ( Hmat,	1.0f,	0.0f,	-3.0f,		 2.0f,
					0.0f,	0.0f,	3.0f,		-2.0f,
					0.0f,	1.0f,	-2.0f,		 1.0f,
					0.0f,	0.0f,	-1.0f,		 1.0f );

	V4DF Tvec;
	vector4 ( Tvec, 1.0f, t, t*t, t*t*t );

	V4DF HT;
	mult_mat_pt_4 ( HT, Hmat, Tvec );

	V4DF g1, g2;
	vector4 ( g1, p1[0], p2[0], v1[0], v2[0] );
	vector4 ( g2, p1[1], p2[1], v1[1], v2[1] );

	pt[0] = vdot4 ( g1, HT );
	pt[1] = vdot4 ( g2, HT );
*/
	V4DF HT;
	vector4 ( HT, 1.0f - 3.0f * t * t + 2.0f * t * t * t,
				  3.0f * t * t - 2.0f * t * t * t,
				  t - 2.0f * t * t + t * t * t,
				  -1.0f * t * t + t * t * t );

	V4DF g1, g2;
	vector4 ( g1, p1[0], p2[0], v1[0], v2[0] );
	vector4 ( g2, p1[1], p2[1], v1[1], v2[1] );

	pt[0] = vdot4 ( g1, HT );
	pt[1] = vdot4 ( g2, HT );
}

void hermite2D ( V2DD pt, double t, V2DD p1, V2DD p2, V2DD v1, V2DD v2 )
{
	V4DD HT;
	vector4 ( HT, 1.0 - 3.0 * t * t + 2.0 * t * t * t,
				  3.0 * t * t - 2.0 * t * t * t,
				  t - 2.0 * t * t + t * t * t,
				  -1.0 * t * t + t * t * t );

	V4DD g1, g2;
	vector4 ( g1, p1[0], p2[0], v1[0], v2[0] );
	vector4 ( g2, p1[1], p2[1], v1[1], v2[1] );

	pt[0] = vdot4 ( g1, HT );
	pt[1] = vdot4 ( g2, HT );
}

void hermite3D ( V3DF pt, float t, V3DF p1, V3DF p2, V3DF v1, V3DF v2 )
{
/*
	Matrix4DF Hmat;
	Matrix ( Hmat,	1.0f,	0.0f,	-3.0f,		 2.0f,
					0.0f,	0.0f,	3.0f,		-2.0f,
					0.0f,	1.0f,	-2.0f,		 1.0f,
					0.0f,	0.0f,	-1.0f,		 1.0f );

	V4DF Tvec;
	vector4 ( Tvec, 1.0f, t, t*t, t*t*t );

	int i, j;
	V4DF HT;
	mult_mat_pt_4 ( HT, Hmat, Tvec );

	V4DF g1, g2, g3;
	vector4 ( g1, p1[0], p2[0], v1[0], v2[0] );
	vector4 ( g2, p1[1], p2[1], v1[1], v2[1] );
	vector4 ( g3, p1[2], p2[2], v1[2], v2[2] );

	pt[0] = vdot ( g1, HT );
	pt[1] = vdot ( g2, HT );
	pt[2] = vdot ( g3, HT );
*/
	V4DF HT;
	vector4 ( HT, 1.0f - 3.0f * t * t + 2.0f * t * t * t,
				  3.0f * t * t - 2.0f * t * t * t,
				  t - 2.0f * t * t + t * t * t,
				  -1.0f * t * t + t * t * t );

	V4DF g1, g2, g3;
	vector4 ( g1, p1[0], p2[0], v1[0], v2[0] );
	vector4 ( g2, p1[1], p2[1], v1[1], v2[1] );
	vector4 ( g3, p1[2], p2[2], v1[2], v2[2] );

	pt[0] = vdot4 ( g1, HT );
	pt[1] = vdot4 ( g2, HT );
	pt[2] = vdot4 ( g3, HT );
}

void hermite3D ( V3DD pt, double t, V3DD p1, V3DD p2, V3DD v1, V3DD v2 )
{
	V4DD HT;
	vector4 ( HT, 1.0 - 3.0 * t * t + 2.0 * t * t * t,
				  3.0 * t * t - 2.0 * t * t * t,
				  t - 2.0 * t * t + t * t * t,
				  -1.0 * t * t + t * t * t );

	V4DD g1, g2, g3;
	vector4 ( g1, p1[0], p2[0], v1[0], v2[0] );
	vector4 ( g2, p1[1], p2[1], v1[1], v2[1] );
	vector4 ( g3, p1[2], p2[2], v1[2], v2[2] );

	pt[0] = vdot4 ( g1, HT );
	pt[1] = vdot4 ( g2, HT );
	pt[2] = vdot4 ( g3, HT );
}

//	k = |x' y'' - x'' y'|/(x'^2 + y'^2)^(3/2)
float hermite2DCurvature ( float t, V2DF p1, V2DF p2, V2DF v1, V2DF v2 )
{
	V2DF pt1, pt2;
	V4DF HT1, HT2;
	vector4 ( HT1, -6.0f * t + 6.0f * t * t,
				    6.0f * t - 6.0f * t * t,
				   1.0f - 4.0f * t + 3.0f * t * t,
				   -2.0f * t + 3.0f * t * t );
	vector4 ( HT2, -6.0f + 12.0f * t,
				    6.0f - 12.0f * t,
				   -4.0f + 6.0f * t,
				   -2.0f + 6.0f * t );

	V4DF g1, g2;
	vector4 ( g1, p1[0], p2[0], v1[0], v2[0] );
	vector4 ( g2, p1[1], p2[1], v1[1], v2[1] );

	pt1[0] = vdot4 ( g1, HT1 );		//	x'
	pt1[1] = vdot4 ( g2, HT1 );		//	y'

	pt2[0] = vdot4 ( g1, HT2 );		//	x''
	pt2[1] = vdot4 ( g2, HT2 );		//	y''

	float norm = sqrt(pt1[0]*pt1[0] + pt1[1]*pt1[1]);
	return fabs ( pt1[0]*pt2[1] - pt2[0]*pt1[1])/ (norm*norm*norm);
}

double hermite2DCurvature ( double t, V2DD p1, V2DD p2, V2DD v1, V2DD v2 )
{
	V2DD pt1, pt2;
	V4DD HT1, HT2;
	vector4 ( HT1, -6.0 * t + 6.0 * t * t,
				    6.0 * t - 6.0 * t * t,
				   1.0 - 4.0 * t + 3.0 * t * t,
				   -2.0 * t + 3.0 * t * t );
	vector4 ( HT2, -6.0 + 12.0 * t,
				    6.0 - 12.0 * t,
				   -4.0 + 6.0 * t,
				   -2.0 + 6.0 * t );

	V4DD g1, g2;
	vector4 ( g1, p1[0], p2[0], v1[0], v2[0] );
	vector4 ( g2, p1[1], p2[1], v1[1], v2[1] );

	pt1[0] = vdot4 ( g1, HT1 );		//	x'
	pt1[1] = vdot4 ( g2, HT1 );		//	y'

	pt2[0] = vdot4 ( g1, HT2 );		//	x''
	pt2[1] = vdot4 ( g2, HT2 );		//	y''

	double norm = sqrt(pt1[0]*pt1[0] + pt1[1]*pt1[1]);
	return fabs ( pt1[0]*pt2[1] - pt2[0]*pt1[1])/ (norm*norm*norm);
}

int Remove_set ( int n_set, int *set, int elt )
{
	int i;
	int cnt;

	for ( i = 0, cnt = 0; i < n_set; i++ ) {
		if ( set[i] == elt ) {
			set[i] = -1;
			cnt++;
		}
	}
	int ni;
	int *nset = (int *) calloc ( n_set, sizeof(int) );
	for ( i = 0, ni = 0; i < n_set; i++ ) {
		if ( set[i] != -1 ) {
			nset[ni++] = set[i];
		}
	}

	for ( i = 0; i < ni; i++ )
		set[i] = nset[i];

	return ni;
}

//	Find a rotation matrix that rotates v1 to v2
void build_rot_mat ( Matrix4DF rmat, V3DF v1, V3DF v2 )
{
	V3DF axis;
	float angle;
	V3DF org;

	vzero ( org );

	nvcross ( axis, v1, v2 );
	angle = angle_r ( v1, v2 );

	get_rot_mat ( rmat, axis, org, angle );

}

//	returns 0, if no intersection
//	returns 1, if one vertex of the triangle is on the plane (p1)
//	returns 2, if two vertices of the triangle are on the plane (p1 & p2)
//	returns 3, if the triangle is on the plane 
//	returns 4, if vertex t1 and edge (t2, t3) are on the plane
//	returns 5, if vertex t2 and edge (t1, t3) are on the plane
//	returns 6, if vertex t3 and edge (t1, t2) are on the plane
//	returns 7, if edge (t1, t2) and edge (t1, t3) are on the plane
//	returns 8, if edge (t2, t1) and edge (t2, t3) are on the plane
//	returns 9, if edge (t3, t1) and edge (t3, t2) are on the plane
int intersection_triangle_plane ( V3DF p1, V3DF p2, V3DF t1, V3DF t2, V3DF t3, V3DF pt, V3DF nm )
{
	V3DF v1, v2, v3;

	vector ( v1, t1, pt );
	vector ( v2, t2, pt );
	vector ( v3, t3, pt );

	if ( (vdot (v1, nm) == 0.0f) && (vdot (v2, nm) == 0.0f) && (vdot (v3, nm) == 0.0f) )
		return 3;

	if ( (vdot (v1, nm) > 0.0f && vdot (v2, nm) > 0.0f && vdot (v3, nm) > 0.0f) ||
		 (vdot (v1, nm) < 0.0f && vdot (v2, nm) < 0.0f && vdot (v3, nm) < 0.0f) )	//	No intersection
		 return 0;

	if ( vdot (v1, nm) == 0.0f ) {								//	t1 is on the plane
		if ( vdot (v2, nm) == 0.0f ) {
			vcopy ( p1, t1 );
			vcopy ( p2, t2 );
			return 2;
		}

		if ( vdot (v3, nm) == 0.0f ) {
			vcopy ( p1, t1 );
			vcopy ( p2, t3 );
			return 2;
		}

		if ( vdot (v2, nm) * vdot (v3, nm) > 0.0f ) {
			vcopy ( p1, t1 );
			return 1;
		}

		vcopy ( p1, t1 );
		if ( get_an_intersection_bw_edge_N_plane_pure ( p2, t2, t3, pt, nm ) ) {
			return 4;
		}

		return -1;
	}

	if ( vdot (v2, nm) == 0.0f ) {								//	t2 is on the plane
		if ( vdot (v3, nm) == 0.0f ) {
			vcopy ( p1, t2 );
			vcopy ( p2, t3 );
			return 2;
		}

		if ( vdot (v1, nm) == 0.0f ) {
			vcopy ( p1, t1 );
			vcopy ( p2, t2 );
			return 2;
		}

		if ( vdot (v1, nm) * vdot (v3, nm) > 0.0f ) {
			vcopy ( p1, t2 );
			return 1;
		}

		vcopy ( p1, t2 );
		if ( get_an_intersection_bw_edge_N_plane_pure ( p2, t1, t3, pt, nm ) ) {
			return 5;
		}

		return -1;
	}

	if ( vdot ( v3, nm) == 0.0f ) { 							//	t3 is on the plane
		if ( vdot (v1, nm) == 0.0f ) {
			vcopy ( p1, t1 );
			vcopy ( p2, t3 );
			return 2;
		}

		if ( vdot (v2, nm) == 0.0f ) {
			vcopy ( p1, t2 );
			vcopy ( p2, t3 );
			return 2;
		}

		if ( vdot (v1, nm) * vdot (v2, nm) > 0.0f ) {
			vcopy ( p1, t3 );
			return 1;
		}

		vcopy ( p1, t3 );
		if ( get_an_intersection_bw_edge_N_plane_pure ( p2, t1, t2, pt, nm ) ) {
			return 6;
		}

		return -1;
	}

	if ( (vdot (v1, nm) > 0.0f && vdot(v2, nm) < 0.0f && vdot(v3, nm) < 0.0f) || (vdot (v1, nm) < 0.0f && vdot(v2, nm) > 0.0f && vdot(v3, nm) > 0.0f) ) {
		if ( get_an_intersection_bw_edge_N_plane_pure ( p1, t1, t2, pt, nm ) && 
			 get_an_intersection_bw_edge_N_plane_pure ( p2, t1, t3, pt, nm ) ) 
			 return 7;
	}

	if ( (vdot (v2, nm) > 0.0f && vdot(v3, nm) < 0.0f && vdot(v1, nm) < 0.0f) || (vdot (v2, nm) < 0.0f && vdot(v3, nm) > 0.0f && vdot(v1, nm) > 0.0f) ) {
		if ( get_an_intersection_bw_edge_N_plane_pure ( p1, t1, t2, pt, nm ) && 
			 get_an_intersection_bw_edge_N_plane_pure ( p2, t2, t3, pt, nm ) ) 
			 return 8;
	}

	if ( (vdot (v3, nm) > 0.0f && vdot(v1, nm) < 0.0f && vdot(v2, nm) < 0.0f) || (vdot (v3, nm) < 0.0f && vdot(v1, nm) > 0.0f && vdot(v2, nm) > 0.0f) ) {
		if ( get_an_intersection_bw_edge_N_plane_pure ( p1, t1, t3, pt, nm ) && 
			 get_an_intersection_bw_edge_N_plane_pure ( p2, t2, t3, pt, nm ) ) 
			 return 9;
	}
	return -1;
}

float get_avg ( int n, float *set )
{
	int i;
	float avg;

	if ( n == 0 )
		return 0.0f;

	for ( i = 0, avg = 0.0f; i < n; i++ )
		avg += set[i];

	return avg / (float) n;
}

float get_std ( int n, float *set )
{
	int i;
	float avg, std;

	if ( n == 0 )
		return 0.0f;

	avg = get_avg ( n, set );
	for ( i = 0, std = 0.0f; i < n; i++ )
		std += (set[i] - avg)*(set[i] - avg);

	return sqrt(std / (float) n);
}

int get_extreme ( V2DF xtr, V2DF axtr, V2DF p0, V2DF p1, V2DF p2, V2DF p3 )
{
	V2DF v1, v2;
	V2DF sol;		//	solution of the differentiated equation
	V2DF rsol1;		//	first solution of the D.E. rsol1[0] -> x, rsol1[1] -> y
	V2DF rsol2;		//	second solution of the D.E. rsol2[0] -> x, rsol2[1] -> y
	float avgv;

	vector2 ( v1, p2, p0 );
	vscale2 ( v1, 0.5f );
	vector2 ( v2, p3, p1 );
	vscale2 ( v2, 0.5f );

//	Build a hermite curve using p1, p2, v1, v2
	float a, b, c;

	a = 6.0f * p1[1] - 6.0f * p2[1] + 3.0f * v1[1] + 3.0f * v2[1];
	b = -6.0f * p1[1] + 6.0f * p2[1] - 4.0f * v1[1] - 2.0f * v2[1];
	c = v1[1];

//	get solution
	vavg2 ( rsol1, p1, p2 );
	vavg2 ( rsol2, p1, p2 );
	avgv = (p1[1] + p2[1])/2.0f;
	vector2 ( xtr, avgv, avgv );
	vector2 ( axtr, -1.0f, -1.0f );

	//	sol[0] > sol[1]
	//	a > 0 ==> sol[1] for emax & sol[0] for emin
	//	a < 0 ==> sol[1] for emin & sol[0] for emax
	int flag = 0;
	if ( quadratic_polynomial ( sol, a, b, c ) > 0 ) {
		if ( (0.0f <= sol[0]) && (sol[0] <= 1.0f) ) {
			hermite2D ( rsol1, sol[0], p1, p2, v1, v2 );
			if ( a > 0.0f ) {			//	determination of local maximum
				xtr[1] = rsol1[1];
				axtr[1] = rsol1[0];
				flag = 1;
			}
			else if ( a < 0.0f ) {		//	determination of local minimum
				xtr[0] = rsol1[1];
				axtr[0] = rsol1[0];
				flag = 2;
			}
		}
		if ( (0.0f <= sol[1]) && (sol[1] <= 1.0f) ) {
			hermite2D ( rsol2, sol[1], p1, p2, v1, v2 );
			if ( a < 0.0f ) {			//	determination of extreme maximum
				xtr[1] = rsol2[1];
				axtr[1] = rsol2[0];
				flag = 1;
			}
			else if ( a > 0.0f ) {		//	determination of local minimum
				xtr[0] = rsol2[1];
				axtr[0] = rsol2[0];
				flag = 2;
			}
		}
	}

	return flag;
}

void get_max_grads ( V2DF max_grad, int n_grad, V2DF *grads )
{
	int i;
	int a, b, c, d;
	V2DF xtr, axtr;
	int cnt = n_grad;
	int flag;

	for ( i = 0; i < n_grad; i++ ) {
		a = (i - 1 + n_grad)%n_grad;
		b = i;
		c = (i + 1)%n_grad;
		d = (i + 2)%n_grad;

		flag = get_extreme ( xtr, axtr, grads[a], grads[b], grads[c], grads[d] );
		if ( flag == 1 ) {
			grads[cnt][0] = axtr[1];
			grads[cnt][1] = xtr[1];
			cnt++;
		}
	}

	int maxi;
	float maxv;
	for ( i = 0, maxi = -1, maxv = 0.0f; i < cnt; i++ ) {
		if ( grads[i][1] > maxv ) {
			maxi = i;
			maxv = grads[i][1];
		}
	}

	vcopy2 ( max_grad, grads[maxi] );
}

float pdist2 ( V2DI pt1, V2DI pt2 )
{
	int i;
	float tmp;

	for ( i = 0, tmp = 0.0; i < 2; i++ )
		tmp += (pt1[i]-pt2[i])*(pt1[i]-pt2[i]);

	return (float) sqrt ( tmp );
}

////////////////////CGVM_MATH ���� �߰�////////////////////

void build_zmat ( Matrix4DF zmat, float sx, float sy, float sz )
{
	zmat[0][0] = sx;	zmat[0][1] = 0.0f;	zmat[0][2] = 0.0f;	zmat[0][3] = 0.0f;
	zmat[1][0] = 0.0f;	zmat[1][1] = sy;	zmat[1][2] = 0.0f;	zmat[1][3] = 0.0f;
	zmat[2][0] = 0.0f;	zmat[2][1] = 0.0f;	zmat[2][2] = sz;	zmat[2][3] = 0.0f;
	zmat[3][0] = 0.0f;	zmat[3][1] = 0.0f;	zmat[3][2] = 0.0f;	zmat[3][3] = 1.0f;
}

//	vector length square
float vlsquare ( V3DF v )
{
	return vdot ( v, v );
}

//	vector length
float vlength ( V3DF v )
{
	return sqrt ( vlsquare ( v ) );
}

//	vector scalar multiplication
void vscalar ( V3DF v, float sca )
{
	v[0] *= sca;
	v[1] *= sca;
	v[2] *= sca;
}
void vscalar ( V3DI v, int sca )
{
	v[0] *= sca;
	v[1] *= sca;
	v[2] *= sca;
}

//	vector add and assignment
//	dest = dest + src
void vadd ( V3DF dest, V3DF src ) 
{
	dest[0] += src[0];
	dest[1] += src[1];
	dest[2] += src[2];
}

void vzero( cv::Vec3f &v )
{
	v[0] = 0.0f;
	v[1] = 0.0f;
	v[2] = 0.0f;
}

void vcopy( V3DF dst, cv::Vec3f &src ) 
{
	dst[0] = src[0];
	dst[1] = src[1];
	dst[2] = src[2];
}


void vcopy( cv::Vec3f &dst, V3DF src )
{
	dst[0] = src[0];
	dst[1] = src[1];
	dst[2] = src[2];
}
