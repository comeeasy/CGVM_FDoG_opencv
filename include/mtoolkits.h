#pragma once


//------------------
#ifndef _MTOOLKITS_H_
#define _MTOOLKITS_H_
//------------------

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <opencv2/opencv.hpp>


#ifndef PI
#define PI 3.141592653589793238462643383279 
#endif

#ifndef PIF
#define PIF 3.141592653589793238462643383279f
#endif

#ifndef PIF2
#define PIF2 2.0f*PIF
#endif

#ifndef EXPF
#define EXPF 2.718281828f
#endif

#ifndef INFINITY
#define INFINITY 65535	//2^16-1(signed)
#endif

#define RED_WEIGHT 0.299
#define GREEN_WEIGHT 0.587
#define BLUE_WEIGHT 0.114


typedef int V2DI[2];
typedef int V3DI[3];
typedef int V4DI[4];
typedef int V5DI[5];
typedef float V2DF[2];
typedef float V3DF[3];
typedef float V4DF[4];
typedef double V2DD[2];
typedef double V3DD[3];
typedef double V4DD[4];
typedef float Matrix2DF[2][2];
typedef float Matrix3DF[3][3];
typedef float Matrix4DF[4][4];
typedef double Matrix2DD[2][2];
typedef double Matrix3DD[3][3];
typedef double Matrix4DD[4][4];


typedef struct _TRGL {
	V3DF pts[3];
} TRGL;

float avg2f ( float a, float b );
double avg2d ( double a, double b );
float avg3f ( float a, float b, float c );
int in_a_set ( int v, int nSet, int *Set );
int sign ( float v );
void swap2i ( V2DI v1, V2DI v2 );
void swap2f ( V2DF v1, V2DF v2 );
void swap3f ( V3DF v1, V3DF v2 );
void swap ( int *l1, int *l2 );
void swap ( float *l1, float *l2 );
void swap ( double *l1, double *l2 );

void vscale2 ( V2DF v, float t );
void vscale2 ( V2DD v, double t );
void vscale ( V3DF v, float t );
void vscale ( V3DD v, double t );
void vscale ( V3DD tar, double t, V3DD src );
void vscale ( V3DF tar, float t, V3DF src );

void vzero( cv::Vec3f &v );
void vcross ( cv::Vec3f &vc, V3DF v1, V3DF v2 );
void vcopy( V3DF dst, cv::Vec3f &src );
void vcopy( cv::Vec3f &dst, V3DF src );

void vcopy ( V3DF des, V3DF src );
void vcopy ( V3DD des, V3DD src );
void vcopy2 ( V2DF des, V2DF src );
void vcopy2 ( V2DD des, V2DD src );
void vcopy4 ( V4DF des, V4DF src );
void vcopy3i ( V3DI des, V3DI src );
void vcopy2i ( V2DI des, V2DI src );
void vcopy4i ( V4DI des, V4DI src );
void vsum ( V3DF v, V3DF v1 );
void vadd ( V3DF v, V3DF v1, V3DF v2 );
void vadd ( V3DD v, V3DD v1, V3DD v2 );
void vsub ( V3DF v, V3DF v1, V3DF v2 );
void vsub ( V3DD v, V3DD v1, V3DD v2 );
void vminus ( V3DF v, V3DF v1 );
void vzero ( V3DI v );
void vzero4i ( V4DI v );
void vzero ( V3DF v );
void vzero ( V3DD v );
void vnegate ( V3DF v );
void vnegate ( V3DD v );
void vnegative ( V3DF v1, V3DF v2 );
void vnegative ( V3DD v1, V3DD v2 );
void vnorm ( V3DF v );
void vnorm ( V3DD v );
void vector2 ( V2DI v, int f1, int f2 );
void vector2 ( V2DF v, float f1, float f2 );
void vector2 ( V2DD v, double f1, double f2 );
void vector ( V3DF v, float f1, float f2, float f3 );
void vector ( V4DF v, float f1, float f2, float f3, float f4 );
void vector ( V3DD v, double f1, double f2, double f3 );
void vector2 ( V2DF v, V2DF p1, V2DF p2 );
void vector2 ( V2DD v, V2DD p1, V2DD p2 );
void vector ( V3DF v, V3DF p1, V3DF p2 );
void vector ( V3DD v, V3DD p1, V3DD p2 );
void vector4 ( V4DF v, float f1, float f2, float f3, float f4 );
void vector4 ( V4DD v, double f1, double f2, double f3, double f4 );
void nvector ( V3DF v, float f1, float f2, float f3 );
void nvector ( V3DD v, double f1, double f2, double f3 );
void nvector ( V3DF v, V3DI p1, V3DI p2 );
void nvector ( V3DF v, V3DF p1, V3DF p2 );
void nvector ( V3DD v, V3DD p1, V3DD p2 );
void nvector ( V3DD v, V3DF p1, V3DF p2 );
void vcross ( V3DF v, V3DF p1, V3DF p2 );
void vcross ( V3DD v, V3DD p1, V3DD p2 );
void nvcross ( V3DF vc, V3DF v1, V3DF v2 );
void nvcross ( V3DD vc, V3DD v1, V3DD v2 );
void nvcross ( V3DF vc, V3DD v1, V3DF v2 );
void nvcross ( V3DF vc, V3DF v1, V3DD v2 );
float vdot4 ( V4DF v1, V4DF v2 );
double vdot4 ( V4DD v1, V4DD v2 );
float vdot ( V3DF v1, V3DF v2 );
double vdot ( V3DD v1, V3DD v2 );
float vdot_scaled ( V3DF v1, V3DF v2 );
void vectori ( V2DI v, int p1, int p2 );
void vectori ( V3DI v, int p1, int p2, int p3 );
void vectori ( V5DI v, int p1, int p2, int p3, int p4, int p5 );
void vectori ( V3DI v, V3DI v1, V3DI v2 );
void nvectori ( V3DF v, V3DI v1, V3DI v2 );
void vscaled_sum3 ( V3DF w, float wt, V3DF v );
void vscaled_sum4 ( V4DF w, float wt, V4DF v );
float vleng ( V3DF v );
float vleng2 (V3DF v );
float vlsquare ( V3DF v );
float vlength ( V3DF v );
void vscalar ( V3DF v, float sca );
void vscalar ( V3DI v, int sca );
void vadd ( V3DF dest, V3DF src );

void vavg3 ( V3DF mpt, V3DF v1, V3DF v2, V3DF v3 );
void vavg2 ( V2DF mpt, V2DF v1, V2DF v2 );
void vavg2 ( V2DD mpt, V2DD v1, V2DD v2 );
void vavg ( V3DF mpt, V3DF v1, V3DF v2 );
void vavg ( V3DD mpt, V3DD v1, V3DD v2 );
void nvavg ( V3DF vt, V3DF v1, V3DF v2 );
int match_V4DI ( V4DI v1, V4DI v2 );
void transpose ( Matrix4DF mat );

void get_projected_v ( V3DF pvt, V3DF vt, V3DF pt, V3DF nm );
float get_projected ( V3DF ppt, V3DF pt, V3DF v, V3DF pt0 );
double get_projected ( V3DD ppt, V3DD pt, V3DD v, V3DD pt0 );
double get_projected_2 ( V3DD pt, V3DD pt0, V3DD pt1 );
float angle_d ( V3DF p1, V3DF p2, V3DF p3 );
float angle_d_2PIF ( V3DF p1, V3DF p2, V3DF p3, V3DF nm );
float angle_d_2PIF ( V3DF v1, V3DF v2, V3DF nm );
float angle_d ( V3DF v1, V3DF v2 );
double angle_d ( V3DD v1, V3DD v2 );
float angle_r ( V3DF p1, V3DF p2, V3DF p3 );
double angle_r ( V3DD p1, V3DD p2, V3DD p3 );
float angle_r ( V3DF v1, V3DF v2 );
float angle_r_2PIF ( V3DF v1, V3DF v2, V3DF nm );
double angle_r_2PI ( V3DD v1, V3DD v2, V3DD nm );
double angle_r ( V3DD v1, V3DD v2 );
void get_plane ( V3DF p0, V3DF p1, V3DF p2, V3DF nml );
//void get_normal ( V3DF v1, V3DF v2, V3DF v3, V4DF nm );
void get_normal ( V3DF nm, V3DF p1, V3DF p2, V3DF p3 );
void get_normal ( V3DD nm, V3DD p1, V3DD p2, V3DD p3 );
void get_point ( V3DF pt, V3DF p0, float t, V3DF v );
void get_point ( V3DD pt, V3DD p0, double t, V3DD v );
void get_point ( V3DF pt, V3DF p0, float t, V3DD v );

int is_same_vector ( V3DF v1, V3DF v2 );
int is_same_vector ( V3DD v1, V3DD v2 );
int is_same_vector ( V3DI v1, V3DI v2 );
int is_same_vector4i ( V4DI v1, V4DI v2 );
int is_zero_vector ( V3DF v );
int is_zero_vector ( cv::Vec3f &v );
int is_zero_vector ( V3DD v );
int is_near_zero ( float val, float thr );
int is_near_zero ( double val, double thr );
int is_similar_vector ( V3DF p1, V3DF p2, float thr );
int is_similar_vector ( V3DD p1, V3DD p2, double thr );
int is_similar_reverse_vector ( V3DF p1, V3DF p2, float thr );
int In_Out_Test( V3DF pt, V3DF p1, V3DF p2, V3DF p3, V3DF pnorm );
int in_set ( int a, int nset, int *set );

int is_same_triangle ( V3DI t1, V3DI t2 );
int does_triangle_include_edge ( V3DI f, V2DI e );
int does_triangle_include_edge ( V3DI f, int v1, int v2 );

double max_d ( int n, V2DD *list );
double min_d ( int n, V2DD *list );
double max_d ( int n, V3DD *list );
double min_d ( int n, V3DD *list );

int min2i ( int a, int b );
int min4i ( int a, int b, int c, int d );
int max2i ( int a, int b );
float min2f ( float a, float b );
float max2f ( float a, float b );
float max3f ( float a, float b, float c );
float min3f ( float a, float b, float c );
float max4f ( float a, float b, float c, float d );
float min4f ( float a, float b, float c, float d );
float max5f ( float a, float b, float c, float d, float e );
float min5f ( float a, float b, float c, float d, float e );
double min2d ( double a, double b );
double max2d ( double a, double b );
double max4d ( double a, double b, double c, double d );
int max4fi ( float v[4] );
int min4fip ( float v[4] );
int project_round ( float v );

float pdist ( V3DI p1, V3DI p2 );
float pdist ( V3DF pt1, V3DF pt2 );
double pdist ( V3DD pt1, V3DD pt2 );
float pdist2 ( V3DF pt1, V3DF pt2 );
double pdist2 ( V3DD pt1, V3DD pt2 );
float edist ( V3DF pt1, V3DF pt2, V3DF elp );
float ldist( V3DF pt, V3DF lp, V3DF rp);  
float lsdist ( V3DF pt, V3DF lp, V3DF rp);		//	dist from line segment 
float ldist2 ( V3DF pt, V3DF lp, V3DF rp);		//	dist from line
double ldist2 ( V3DD pt, V3DD lp, V3DD rp); 
float ldist3 ( V3DF pt, V3DF p0, V3DF vec );
float tdist(V3DF pt, V3DF p1, V3DF p2, V3DF p3);
float area ( V3DF p1, V3DF p2, V3DF p3 );

void rotate_vtx_x ( V3DF vtx, float ang );
void rotate_vtx_y ( V3DF vtx, float ang );
void rotate_vtx_z ( V3DF vtx, float ang );
void rotate_vtx ( V3DF vtx, V3DF pt, V3DF ang, V3DF org );
void rotate_pt ( V3DF pts, V3DF ang, V3DF org );

void internal_subdivision_3D ( V3DI wpt, int v1, int v2, V3DI p1, V3DI p2 );
void internal_subdivision_2D ( V2DI wpt, int v1, int v2, V2DI p1, V2DI p2 );
void internal_subdivision_3D ( V3DF wpt, float v1, float v2, V3DF p1, V3DF p2 );
void internal_subdivision_2D ( V2DF wpt, float v1, float v2, V2DF p1, V2DF p2 );

void Matrix ( Matrix2DF M, float m11, float m12, float m21, float m22 );
void Matrix ( Matrix2DD m, double m11, double m12, double m21, double m22 );
void Matrix ( Matrix3DF m, float m11, float m12, float m13, float m21, float m22, float m23, float m31, float m32, float m33 );
void Matrix ( Matrix3DD m, double m11, double m12, double m13, double m21, double m22, double m23, double m31, double m32, double m33 );
void Matrix ( Matrix4DF m, float m11, float m12, float m13, float m14, float m21, float m22, float m23, float m24, float m31, float m32, float m33, float m34, float m41, float m42, float m43, float m44 );
void Matrix ( Matrix4DD m, double m11, double m12, double m13, double m14, double m21, double m22, double m23, double m24, double m31, double m32, double m33, double m34, double m41, double m42, double m43, double m44 );
float get_determinant2 ( float m11, float m12, float m21, float m22 );
float get_determinant3 ( float m11, float m12, float m13, float m21, float m22, float m23, float m31, float m32, float m33 );
float get_determinant4 ( float m11, float m12, float m13, float m14, float m21, float m22, float m23, float m24, float m31, float m32, float m33, float m34, float m41, float m42, float m43, float m44 );
float get_determinant ( int n, float **m );
void get_inv_mat ( float **im, int n, float **m );
void get_inverse_4_4 ( Matrix4DD invm, Matrix4DD m );
void get_inverse_4_4 ( Matrix4DF invm, Matrix4DF m );
void get_inverse_3_3 ( Matrix3DD invm, Matrix3DD m );
void get_inverse_3_3 ( Matrix3DF invm, Matrix3DF m );
void mult_mat_pt ( V3DF npt, const Matrix4DF m, const V3DF pt );
void mult_mat_pt_3 ( V3DF npt, Matrix3DF m, V3DF pt );
void mult_mat_pt_3 ( V3DD npt, Matrix3DD m, V3DD pt );
void mult_mat_pt_4 ( V4DF npt, Matrix4DF m, V4DF pt );
void mult_mat_pt_4 ( V4DD npt, Matrix4DD m, V4DD pt );
void mul_mat2 ( Matrix2DF m, Matrix2DF m1, Matrix2DF m2 );
void mul_mat3 ( Matrix3DF m, Matrix3DF m1, Matrix3DF m2 );
void mul_mat4 ( Matrix4DF m, Matrix4DF m1, Matrix4DF m2 );
void rotate_by_matrix ( V3DF rpt, Matrix4DF rot_mat, V3DF pt );
void rotate_by_matrix ( V3DD rpt, Matrix4DD rot_mat, V3DD pt );

void copy_set ( int *des, int nsrc, int *src );
int minus_set2i ( int *mset, int nset1, int nset2, int *set1, int *set2 );
int union_set2i ( int *tset, int num, int nset1, int *set1 );
int union_set2i ( int *uset, int nset1, int nset2, int *set1, int *set2 );
int merge_set2i ( int n_set1, int n_set2, int *set1, int *set2 );
int identical_set ( int n, int *set1, int *set2 );
int identical_set ( int n1, int n2, int *set1, int *set2 );
int identical_set ( int n, short *set1, short *set2 );
int identical_set ( int n1, int n2, short *set1, short *set2 );
int is_element_of_set ( int x, int n, int *set );
int contain_seti ( int n1, int n2, int *set1, int *set2 );

void mult_quat ( V4DF quat, V4DF q1, V4DF q2 );
void build_rot_quat ( V4DF quatq, V3DF pt, V3DF axis, float theta );
void rotate_point ( V3DF p2, V3DF p1, V3DF pt, V3DF axis, float theta );

void mscale ( Matrix2DF m, float s );
void mscale ( Matrix2DD m, double s );
void mscale ( Matrix3DF m, float s );
void mscale ( Matrix3DD m, double s );
void mscale ( Matrix4DF m, float s );
void mscale ( Matrix4DD m, double s );
void madd ( Matrix2DF m, Matrix2DF m1, Matrix2DF m2 );
void madd ( Matrix2DD m, Matrix2DD m1, Matrix2DD m2 );
void madd ( Matrix3DF m, Matrix3DF m1, Matrix3DF m2 );
void madd ( Matrix3DD m, Matrix3DD m1, Matrix3DD m2 );
void madd ( Matrix4DF m, Matrix4DF m1, Matrix4DF m2 );
void madd ( Matrix4DD m, Matrix4DD m1, Matrix4DD m2 );
void msub ( Matrix2DF m, Matrix2DF m1, Matrix2DF m2 );
void msub ( Matrix2DD m, Matrix2DD m1, Matrix2DD m2 );
void msub ( Matrix3DF m, Matrix3DF m1, Matrix3DF m2 );
void msub ( Matrix3DD m, Matrix3DD m1, Matrix3DD m2 );
void msub ( Matrix4DF m, Matrix4DF m1, Matrix4DF m2 );
void msub ( Matrix4DD m, Matrix4DD m1, Matrix4DD m2 );
void mmult ( Matrix2DF m, Matrix2DF m1, Matrix2DF m2 );
void mmult ( Matrix2DD m, Matrix2DD m1, Matrix2DD m2 );
void mmult ( Matrix3DF m, Matrix3DF m1, Matrix3DF m2 );
void mmult ( Matrix3DD m, Matrix3DD m1, Matrix3DD m2 );
void mmult ( Matrix4DF m, Matrix4DF m1, Matrix4DF m2 );
void mmult ( Matrix4DD m, Matrix4DD m1, Matrix4DD m2 );
void mtrans ( Matrix2DF mt, Matrix2DF m );
void mtrans ( Matrix2DD mt, Matrix2DD m );
void mtrans ( Matrix3DF mt, Matrix3DF m );
void mtrans ( Matrix3DD mt, Matrix3DD m );
void mtrans ( Matrix4DF mt, Matrix4DF m );
void mtrans ( Matrix4DD mt, Matrix4DD m );
void load_zero ( Matrix2DF mat );
void load_zero ( Matrix2DD mat );
void load_zero ( Matrix3DF mat );
void load_zero ( Matrix3DD mat );
void load_zero ( Matrix4DF mat );
void load_zero ( Matrix4DD mat );
void load_identity ( Matrix2DF mat );
void load_identity ( Matrix2DD mat );
void load_identity ( Matrix3DF mat );
void load_identity ( Matrix3DD mat );
void load_identity ( Matrix4DF mat );
void load_identity ( Matrix4DD mat );
void get_matrix ( Matrix3DF mat, V3DF v );
void get_matrix ( Matrix3DD mat, V3DD v );
void get_matrix ( Matrix4DF mat, V4DF v );
void get_matrix ( Matrix4DD mat, V4DD v );
void multiply_transposed_matrix ( Matrix2DF QMQmat, Matrix2DF Qmat, Matrix2DF Mmat );
void multiply_transposed_matrix ( Matrix2DD QMQmat, Matrix2DD Qmat, Matrix2DD Mmat );
void multiply_transposed_matrix ( Matrix3DF QMQmat, Matrix3DF Qmat, Matrix3DF Mmat );
void multiply_transposed_matrix ( Matrix3DD QMQmat, Matrix3DD Qmat, Matrix3DD Mmat );
void multiply_transposed_matrix ( Matrix4DF QMQmat, Matrix4DF Qmat, Matrix4DF Mmat );
void multiply_transposed_matrix ( Matrix4DD QMQmat, Matrix4DD Qmat, Matrix4DD Mmat );
void translatemat ( Matrix4DF mat, float p0, float p1, float p2 );
void translatemat ( Matrix4DD mat, double p0, double p1, double p2 );
void rotatemat ( Matrix4DF mat, float ang, V3DF axis );
void rotatemat ( Matrix4DD mat, double ang, V3DD axis );
void mul_mat ( Matrix4DF mat1, Matrix4DF mat2 );
void mul_mat ( Matrix4DD mat1, Matrix4DD mat2 );
void get_rot_mat ( Matrix4DF mat, V3DF a, V3DF p, float phi );
void get_rot_mat ( Matrix4DD mat, V3DD a, V3DD p, double phi );
void get_rot_mat ( Matrix4DF mat, Matrix3DF mat3, V3DF p );
void get_rot_mat ( Matrix4DF mat, float src1, float src2, float src3, float dst1, float dst2, float dst3, float p1, float p2, float p3 );
void get_rot_mat ( Matrix4DF mat, V3DF src, V3DF dst );
void vector_rotate ( V3DF tvec, V3DF svec, V3DF axis, V3DF cent, float ang );
void vector_rotate ( V3DD tvec, V3DD svec, V3DD axis, V3DD cent, double ang );
void vector_rotate ( V3DF tvec, V3DF svec, Matrix4DF mat );
void vector_rotate ( V3DF tvec, V3DF svec, V3DF src, V3DF dst );
void mzero3 ( Matrix3DF m );
void mzero4 ( Matrix4DF m );

void V3DF_interpolate(V3DF** src, V3DF dst, float px, float py, int w, int h);
void V3DF_interpolate(cv::Mat &src, V3DF dst, float px, float py, int w, int h);
float float_interpolate(cv::Mat &src, float px, float py, int w, int h);

void get_point_from_barycentric ( V3DF pt, V3DF bary, V3DF p1, V3DF p2, V3DF p3 );
void get_normals ( V3DF *TrNmls, int nTrs, int nTrPts, V3DI *Trs, V3DF *TrPts );
int is_a_projected_point_on_a_directed_triangle ( V3DF bary, V3DF pt, V3DF p1, V3DF p2, V3DF p3 );
int get_a_projected_point_on_a_plane ( V3DF ppt, V3DF pt, V3DF p1, V3DF nml );
int get_a_projected_point_on_a_plane ( V3DD ppt, V3DD pt, V3DD p1, V3DD nml );
void get_barycentric ( V3DF bary, V3DF pt, V3DF p1, V3DF p2, V3DF p3 );

void compute_intersection_xy ( V3DF pos, V3DF p1, V3DF p2, V3DF q1, V3DF q2 );
void compute_intersection_xy ( V3DD pos, V3DD p1, V3DD p2, V3DD q1, V3DD q2 );

void error_exit ( int value, char *str );
void error_report ( int value, char *str );

int intersection_circle_line_segment ( V3DF ip, V3DF p1, V3DF p2, V3DF p, float radius, int dir );

void build_change_of_basis ( Matrix3DD TSmat, Matrix3DD Smat, Matrix3DD Tmat );


float area_triangle ( V3DF p1, V3DF p2, V3DF p3 );
double area_triangle ( V3DD p1, V3DD p2, V3DD p3 );
float get_cosine ( V3DF p1, V3DF p2, V3DF p3 );

int is_a_point_above_a_plane ( V3DF pt, V3DF p, V3DF n );
int get_an_intersection_bw_edge_N_plane ( V3DF ip, V3DF p1, V3DF p2, V3DF p, V3DF n );
int get_an_intersection_bw_edge_N_plane_pure ( V3DF ip, V3DF p1, V3DF p2, V3DF p, V3DF n );
int get_an_intersection_line_linesegment ( V3DF ip, V3DF pt, V3DF vec, V3DF p1, V3DF p2 );
int get_an_intersection_line_linesegment_nonplanar ( V3DF ip, V3DF p0, V3DF vec, V3DF p1, V3DF p2, V3DF vec1 );

int between ( float l, float v, float r );
int cube_p_intersects_cube_q ( V3DF p1, V3DF p2, V3DF p3, V3DF p4, V3DF p5, V3DF p6, V3DF p7, V3DF p8, 
							   V3DF q1, V3DF q2, V3DF q3, V3DF q4, V3DF q5, V3DF q6, V3DF q7, V3DF q8 );

int inside_triangle ( V3DF bary, V3DF pt, V3DF p1, V3DF p2, V3DF p3 );
int an_edge_intersects_a_triangle ( V3DF p1, V3DF p2, V3DF t1, V3DF t2, V3DF t3 );
int an_edge_intersects_a_triangle ( V3DF ip, V3DF p1, V3DF p2, V3DF t1, V3DF t2, V3DF t3 );
int X_ray_intersects_a_triangle ( V3DF ip, V3DF p1, V3DF p2, V3DF t1, V3DF t2, V3DF t3 );
int Y_ray_intersects_a_triangle ( V3DF ip, V3DF p1, V3DF p2, V3DF t1, V3DF t2, V3DF t3 );
int Z_ray_intersects_a_triangle ( V3DF ip, V3DF p1, V3DF p2, V3DF t1, V3DF t2, V3DF t3 );

int does_edge_intersect_plane ( V3DF p1, V3DF p2, V3DF n, V3DF p );

int Union_set ( int *uset, int n_set1, int n_set2, int *set1, int *set2 );
int Union_set ( int *uset, int n_set1, int n_set2, int *set2 );
int Remove_set ( int n_set, int *set, int elt );

void interpolate_point ( V3DF p, V3DF p1, V3DF p2, int a, int b );
void interpolate_point ( V3DF p, V3DF p1, V3DF p2, float a, float b );

void get_interpolation_ratio ( V2DF irat, V3DF va, V3DF vb, V3DF vt );

int quadratic_polynomial ( V2DF sol, float a, float b, float c );
int quadratic_polynomial ( V2DD sol, double a, double b, double c );

void hermite2D ( V2DF pt, float t, V2DF p1, V2DF p2, V2DF v1, V2DF v2 );
void hermite2D ( V2DD pt, double t, V2DD p1, V2DD p2, V2DD v1, V2DD v2 );

void hermite3D ( V3DF pt, float t, V3DF p1, V3DF p2, V3DF v1, V3DF v2 );
void hermite3D ( V3DD pt, double t, V3DD p1, V3DD p2, V3DD v1, V3DD v2 );

float hermite2DCurvature ( float t, V2DF p1, V2DF p2, V2DF v1, V2DF v2 );
double hermite2DCurvature ( double t, V2DD p1, V2DD p2, V2DD v1, V2DD v2 );


void build_rot_mat ( Matrix4DF rmat, V3DF v1, V3DF v2 );

int intersection_triangle_plane ( V3DF p1, V3DF p2, V3DF t1, V3DF t2, V3DF t3, V3DF pt, V3DF nm );

float get_avg ( int n, float *set );
float get_std ( int n, float *set );

void get_max_grads ( V2DF max_grad, int n_grad, V2DF *grads );

float pdist2 ( V2DI pt1, V2DI pt2 );

////////////////////CGVM_MATH ���� �߰�////////////////////
void build_zmat ( Matrix4DF zmat, float sx, float sy, float sz );

//------------------
#endif
//------------------