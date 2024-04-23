#include "image_processing.h"


void apply_FBL_filter(
    V3DF** src_cim, 
	FlowPath** src_fpath,
    V3DF** cimFBL,
    int w, int h,
    float sigma_e, float gamma_e, float sigma_g, float gamma_g,	int iteration,
    int threshold_T
)
{	// Flow-Based Bilateral Filter
	printf("FBLFilter..");
	int it, i,j, r;
	float px,py;
	V3DF clr,dclr;
	int MaskSize=threshold_T;
	float G_sigma_e, G_sigma_g;
	float H_gamma_e, H_gamma_g;
	V3DF Ce, Cg;
	float normal_term;

	V3DF **InImg = new V3DF *[w];
	V3DF **OImg = new V3DF *[w];
	for (i=0; i<w; i++){
		InImg[i] = new V3DF[h];
		OImg[i] = new V3DF[h];
		for (j=0; j<h; j++){
			vcopy(InImg[i][j], src_cim[i][j]);
			vcopy(OImg[i][j], src_cim[i][j]);
		}
	}
	conv_RGBtoLAB(InImg, w, h);
	conv_RGBtoLAB(OImg, w, h);

	float s, div_g,div_g2;
	float cdist, div_s, div_s2;	// color distance
	for (it=0; it<iteration; it++){
		printf("%d..",it+1);
		// 1. Ce : linear bilateral filter along the edge(or ETF)
		div_g	= sqrt(2*PI)*sigma_e;	div_g2	= 2*sigma_e*sigma_e;
		div_s	= sqrt(2*PI)*gamma_e;	div_s2	= 2*gamma_e*gamma_e;
		for (i=0; i<w; i++){
			for (j=0; j<h; j++){
				normal_term = 0.0f;	vzero(Ce);
				for (r=0; r<src_fpath[i][j].sn; r++){
					px = src_fpath[i][j].Alpha[r][0];	py = src_fpath[i][j].Alpha[r][1];
					//if ( px<0||px>=MAX_X || py<0||py>=MAX_Y)	continue;
					// 1) color
					GetValAtPoint_V3DF(InImg, w, h, px,py, clr);

					// 2) spatial
					s =	sqrt((i-px)*(i-px) + (j-py)*(j-py));	// (i,j)~(px,py)
					G_sigma_e = ( exp(-(s*s)/div_g2) / div_g );

					// 3) color
					vsub(dclr, InImg[i][j], clr);	// clr(i,j) ~ clr(px,py)
					cdist = vlength(dclr);
					H_gamma_e = ( exp(-(cdist*cdist)/div_s2) / div_s );

					// 4) 최종
					vscalar(clr, G_sigma_e * H_gamma_e);
					vadd(Ce, clr);
					normal_term += (G_sigma_e * H_gamma_e);
				}
				if ( normal_term==0.0f )	continue;

				vscalar(Ce, 1.0f/normal_term);
				vcopy(OImg[i][j], Ce);
			}
		}

		init_buf(InImg, OImg, w, h);

		// 2. Cg : linear bilateral filter along the gradient direction
		div_g	= sqrt(2*PI)*sigma_g;	div_g2	= 2*sigma_g*sigma_g;
		div_s	= sqrt(2*PI)*gamma_g;	div_s2	= 2*gamma_g*gamma_g;
		for (i=0; i<w; i++){
			for (j=0; j<h; j++){
				normal_term = 0.0f;	vzero(Cg);
				for (r=0; r<src_fpath[i][j].tn; r++){
					px = src_fpath[i][j].Beta[r][0];	py = src_fpath[i][j].Beta[r][1];
					//if ( px<0||px>=MAX_X || py<0||py>=MAX_Y)	continue;
					// 1) color
					GetValAtPoint_V3DF(InImg, w,h, px,py, clr);

					// 2) spatial
					s =	sqrt((i-px)*(i-px) + (j-py)*(j-py));	// (i,j)~(px,py)
					G_sigma_g = ( exp(-(s*s)/div_g2) / div_g );

					// 3) color
					vsub(dclr, InImg[i][j], clr);	// clr(i,j) ~ clr(px,py)
					cdist = vlength(dclr);
					H_gamma_g = ( exp(-(cdist*cdist)/div_s2) / div_s );

					// 4) 최종
					vscalar(clr, G_sigma_g * H_gamma_g);
					vadd(Cg, clr);
					normal_term += (G_sigma_g *H_gamma_g);
				}
				if ( normal_term==0.0f )	continue;

				vscalar(Cg, 1.0f/normal_term);
				vcopy(OImg[i][j], Cg);
			}
		}

		init_buf(InImg, OImg, w, h);
	}

	conv_LABtoRGB(OImg, w, h);
	for (i=0; i<w; i++)
		for (j=0; j<h; j++)
			vcopy(cimFBL[i][j], OImg[i][j]);

	for (i=0; i<w; i++){
		delete [] InImg[i];
		delete [] OImg[i];
	}
	delete [] InImg;
	delete [] OImg;

	printf("done\n");
}
void apply_FBL_filter(
    cv::Mat &src_cim, 
	FlowPath** src_fpath,
    cv::Mat &dst_cimFBL,
    int w, int h,
    float sigma_e, float gamma_e, float sigma_g, float gamma_g,	int iteration,
    int threshold_T
)
{	// Flow-Based Bilateral Filter
	printf("FBLFilter..");
	int it, i,j, r;
	float px,py;
	V3DF clr,dclr;

	int MaskSize=threshold_T;
	float G_sigma_e, G_sigma_g;
	float H_gamma_e, H_gamma_g;
	V3DF Ce, Cg;
	float normal_term;

	V3DF **InImg = new V3DF *[w];
	V3DF **OImg = new V3DF *[w];
	for (i=0; i<w; i++){
		InImg[i] = new V3DF[h];
		OImg[i] = new V3DF[h];
		for (j=0; j<h; j++){
			vcopy(InImg[i][j], src_cim.at<cv::Vec3f>(j, i));
			vcopy(OImg[i][j], src_cim.at<cv::Vec3f>(j, i));
		}
	}
	conv_RGBtoLAB(InImg, w, h);
	conv_RGBtoLAB(OImg, w, h);

	float s, div_g,div_g2;
	float cdist, div_s, div_s2;	// color distance
	for (it=0; it<iteration; it++){
		printf("%d..",it+1);
		// 1. Ce : linear bilateral filter along the edge(or ETF)
		div_g	= sqrt(2*PI)*sigma_e;	div_g2	= 2*sigma_e*sigma_e;
		div_s	= sqrt(2*PI)*gamma_e;	div_s2	= 2*gamma_e*gamma_e;
		for (i=0; i<w; i++){
			for (j=0; j<h; j++){
				normal_term = 0.0f;	vzero(Ce);
				for (r=0; r<src_fpath[i][j].sn; r++){
					px = src_fpath[i][j].Alpha[r][0];	py = src_fpath[i][j].Alpha[r][1];
					//if ( px<0||px>=MAX_X || py<0||py>=MAX_Y)	continue;
					// 1) color
					GetValAtPoint_V3DF(InImg, w, h, px,py, clr);

					// 2) spatial
					s =	sqrt((i-px)*(i-px) + (j-py)*(j-py));	// (i,j)~(px,py)
					G_sigma_e = ( exp(-(s*s)/div_g2) / div_g );

					// 3) color
					vsub(dclr, InImg[i][j], clr);	// clr(i,j) ~ clr(px,py)
					cdist = vlength(dclr);
					H_gamma_e = ( exp(-(cdist*cdist)/div_s2) / div_s );

					// 4) 최종
					vscalar(clr, G_sigma_e * H_gamma_e);
					vadd(Ce, clr);
					normal_term += (G_sigma_e * H_gamma_e);
				}
				if ( normal_term==0.0f )	continue;

				vscalar(Ce, 1.0f/normal_term);
				vcopy(OImg[i][j], Ce);

				printf("Ce: [%.3f, %.3f, %.3f]\n", Ce[0], Ce[1], Ce[2]);
			}
		}

		init_buf(InImg, OImg, w, h);

		// 2. Cg : linear bilateral filter along the gradient direction
		div_g	= sqrt(2*PI)*sigma_g;	div_g2	= 2*sigma_g*sigma_g;
		div_s	= sqrt(2*PI)*gamma_g;	div_s2	= 2*gamma_g*gamma_g;
		for (i=0; i<w; i++){
			for (j=0; j<h; j++){
				normal_term = 0.0f;	vzero(Cg);
				for (r=0; r<src_fpath[i][j].tn; r++){
					px = src_fpath[i][j].Beta[r][0];	py = src_fpath[i][j].Beta[r][1];
					//if ( px<0||px>=MAX_X || py<0||py>=MAX_Y)	continue;
					// 1) color
					GetValAtPoint_V3DF(InImg, w,h, px,py, clr);

					// 2) spatial
					s =	sqrt((i-px)*(i-px) + (j-py)*(j-py));	// (i,j)~(px,py)
					G_sigma_g = ( exp(-(s*s)/div_g2) / div_g );

					// 3) color
					vsub(dclr, InImg[i][j], clr);	// clr(i,j) ~ clr(px,py)
					cdist = vlength(dclr);
					H_gamma_g = ( exp(-(cdist*cdist)/div_s2) / div_s );

					// 4) 최종
					vscalar(clr, G_sigma_g * H_gamma_g);
					vadd(Cg, clr);
					normal_term += (G_sigma_g *H_gamma_g);
				}
				if ( normal_term==0.0f )	continue;

				vscalar(Cg, 1.0f/normal_term);
				vcopy(OImg[i][j], Cg);
			}
		}

		init_buf(InImg, OImg, w, h);
	}

	conv_LABtoRGB(OImg, w, h);
	for (i=0; i<w; i++)
		for (j=0; j<h; j++)
			vcopy(dst_cimFBL.at<cv::Vec3f>(j, i), OImg[i][j]);

	for (i=0; i<w; i++){
		delete [] InImg[i];
		delete [] OImg[i];
	}
	delete [] InImg;
	delete [] OImg;

	printf("done\n");
}

void get_gradient(
	float** src_gray_im, 
	V3DF** dst_grad, 
	int w, int h, 
	float grad_thr
) {
	float **imG = new float *[w];	//Gaussian����� gray image
	for (int i=0; i<w; i++){
		imG[i] = new float[h];
	}
	
	GaussianFilter(src_gray_im, imG, w, h);

	int cnt=0;
	int row,col;
	int colOffset,rowOffset;
	int colTotal,rowTotal;
	float GxMask[3][3];				// Sobel mask in the x direction
	float GyMask[3][3];				// Sobel mask in the y direction
	float Gx,Gy;	// Sum of Sobel mask products values in the x,y direction

	/* Declare Sobel masks */
	GxMask[0][2] = -1.0f; GxMask[1][2] = 0.0f; GxMask[2][2] = 1.0f;
	GxMask[0][1] = -2.0f; GxMask[1][1] = 0.0f; GxMask[2][1] = 2.0f;
	GxMask[0][0] = -1.0f; GxMask[1][0] = 0.0f; GxMask[2][0] = 1.0f;

	GyMask[0][2] =  1.0f; GyMask[1][2] =  2.0f; GyMask[2][2] =  1.0f;
	GyMask[0][1] =  0.0f; GyMask[1][1] =  0.0f; GyMask[2][1] =  0.0f;
	GyMask[0][0] = -1.0f; GyMask[1][0] = -2.0f; GyMask[2][0] = -1.0f;

	float vmax = -INFINITY;
	/* Determine edge directions and gradient strengths */
	for (col=0; col<w; col++){	//col
		for (row=0; row<h; row++){	//row
			Gx = 0;
			Gy = 0;
			/* Calculate the sum of the Sobel mask times the nine surrounding pixels in the x and y direction */
			for (colOffset=-1; colOffset<=1; colOffset++) {
				for (rowOffset=-1; rowOffset<=1; rowOffset++) {
					colTotal = col+colOffset;
					rowTotal = row+rowOffset;
					//���ó��(��Ī)
					if (rowTotal<0) rowTotal *= (-1.0f);
					else if (rowTotal>=h)	rowTotal = h-(rowTotal-(h-1));
					if (colTotal<0) colTotal *= (-1.0f);
					else if (colTotal>=w)	colTotal = w-(colTotal-(w-1));

					Gx = Gx + (imG[colTotal][rowTotal] * GxMask[colOffset+1][rowOffset+1]);//??
					Gy = Gy + (imG[colTotal][rowTotal] * GyMask[colOffset+1][rowOffset+1]);//??
				}
			}
			dst_grad[col][row][0] = Gx;
			dst_grad[col][row][1] = Gy;
			dst_grad[col][row][2] = sqrt(Gx*Gx + Gy*Gy);

			if ( dst_grad[col][row][2] < grad_thr ){	//���� gradient ����	//[���]
				vzero(dst_grad[col][row]);
			} 
		}
	}

	for (int i=0; i<w; i++)	delete [] imG[i];	delete	[] imG;
}

void get_gradient(
	cv::Mat &src_gray_img,
	cv::Mat &dst_grad,
	int w, int h, 
	float grad_thr
) {
	float **imG 				= new float *[w];	//Gaussian����� gray image
	float **src_gray_img_arr	= new float *[w];
	for (int i=0; i<w; i++){
		imG[i] 				= new float[h];
		src_gray_img_arr[i] = new float[h];
	}

	for(int i=0; i<w; ++i) {
		for(int j=0; j<h; ++j) {
			src_gray_img_arr[i][j] = src_gray_img.at<float>(j, i);
		}
	}
	
	GaussianFilter(src_gray_img_arr, imG, w, h);

	int cnt=0;
	int row,col;
	int colOffset,rowOffset;
	int colTotal,rowTotal;
	float GxMask[3][3];				// Sobel mask in the x direction
	float GyMask[3][3];				// Sobel mask in the y direction
	float Gx,Gy;	// Sum of Sobel mask products values in the x,y direction

	/* Declare Sobel masks */
	GxMask[0][2] = -1.0f; GxMask[1][2] = 0.0f; GxMask[2][2] = 1.0f;
	GxMask[0][1] = -2.0f; GxMask[1][1] = 0.0f; GxMask[2][1] = 2.0f;
	GxMask[0][0] = -1.0f; GxMask[1][0] = 0.0f; GxMask[2][0] = 1.0f;

	GyMask[0][2] =  1.0f; GyMask[1][2] =  2.0f; GyMask[2][2] =  1.0f;
	GyMask[0][1] =  0.0f; GyMask[1][1] =  0.0f; GyMask[2][1] =  0.0f;
	GyMask[0][0] = -1.0f; GyMask[1][0] = -2.0f; GyMask[2][0] = -1.0f;

	float vmax = -INFINITY;
	/* Determine edge directions and gradient strengths */
	for (col=0; col<w; col++){	//col
		for (row=0; row<h; row++){	//row
			Gx = 0;
			Gy = 0;
			/* Calculate the sum of the Sobel mask times the nine surrounding pixels in the x and y direction */
			for (colOffset=-1; colOffset<=1; colOffset++) {
				for (rowOffset=-1; rowOffset<=1; rowOffset++) {
					colTotal = col+colOffset;
					rowTotal = row+rowOffset;
					//���ó��(��Ī)
					if (rowTotal<0) rowTotal *= (-1.0f);
					else if (rowTotal>=h)	rowTotal = h-(rowTotal-(h-1));
					if (colTotal<0) colTotal *= (-1.0f);
					else if (colTotal>=w)	colTotal = w-(colTotal-(w-1));

					Gx = Gx + (imG[colTotal][rowTotal] * GxMask[colOffset+1][rowOffset+1]);//??
					Gy = Gy + (imG[colTotal][rowTotal] * GyMask[colOffset+1][rowOffset+1]);//??
				}
			}
			cv::Vec3f &pixel = dst_grad.at<cv::Vec3f>(row, col);
			pixel[0] = Gx;
			pixel[1] = Gy;
			pixel[2] = sqrt(Gx*Gx + Gy*Gy);

			if ( pixel[2] < grad_thr ){	//���� gradient ����	//[���]
				pixel[0] = 0;
				pixel[1] = 0;
				pixel[2] = 0;
			} 
		}
	}

	for (int i=0; i<w; i++)	{
		delete [] imG[i];
		delete [] src_gray_img_arr[i];
	}	
	delete	[] imG;
	delete	[] src_gray_img_arr;
}

void get_tangent(
	V3DF** src_grad, 
	V3DF** dst_tangent, 
	int w, int h
) {
	int i,j;
	V3DF gradient;
	V3DF axis_z = {0.0f, 0.0f, 1.0f};

	for (i=0; i<w; i++){
		for (j=0; j<h; j++){
			if ( src_grad[i][j][2]==0.0f ){
				vzero( dst_tangent[i][j] );
				continue;
			}
			vector(gradient, src_grad[i][j][0], src_grad[i][j][1], 0.0f);
			vcross(dst_tangent[i][j], gradient, axis_z);
		}
	}
}

void get_tangent(
	cv::Mat &src_grad, 
	cv::Mat &dst_tangent, 
	int w, int h
) {
	V3DF gradient;
	V3DF axis_z = {0.0f, 0.0f, 1.0f};

	for (int i=0; i<w; i++){
		for (int j=0; j<h; j++){
			cv::Vec3f &src_grad_pixel = src_grad.at<cv::Vec3f>(j, i);
			cv::Vec3f &dst_tangent_pixel = dst_tangent.at<cv::Vec3f>(j, i);
			if ( src_grad_pixel[2]==0.0f ){
				vzero( dst_tangent_pixel );
				continue;
			}
			vector(gradient, src_grad_pixel[0], src_grad_pixel[1], 0.0f);
			vcross(dst_tangent_pixel, gradient, axis_z);
		}
	}
}

void get_ETF(
	V3DF** src_grad, V3DF** &dst_etf,
	int nbhd, int w, int h
) {
	printf("GetETF.. ");
	float start = clock();
	int i,j,a,b;
	V3DI xvec,yvec;
	V3DF tcur,tsum;
	float W;
	float r = (float)nbhd;

	V3DF **oldetf = new V3DF *[w];
	for (i=0; i<w; i++){
		oldetf[i] = new V3DF[h];
		for (j=0; j<h; j++)
			vcopy(oldetf[i][j], dst_etf[i][j]);	//old etf �ʱ�ȭ
	}

	for (i=0; i<w; i++){
		for (j=0; j<h; j++){
			if ( src_grad[i][j][2] == 0.0f ){
				vzero(dst_etf[i][j]);
				continue;
			}

			vectori(xvec, i,j,0);
			vzero(tsum);
			for (a=-nbhd; a<=nbhd; a++){
				for (b=-nbhd; b<=nbhd; b++){
					vectori(yvec, i+a, j+b, 0);

					// ���ó��
					if (yvec[0]<0) yvec[0] *= (-1.0f);
					else if (yvec[0]>=w)	yvec[0] = w-(yvec[0]-(w-1));
					if (yvec[1]<0) yvec[1] *= (-1.0f);
					else if (yvec[1]>=h)	yvec[1] = h-(yvec[1]-(h-1));

					if (src_grad[yvec[0]][yvec[1]][2] == 0.0f)	continue;

					//1. w_s(x,y) = 0, if |xvec-yvec| >= r
					//if ( pdist2(xvec,yvec) >= r)
					if ( pdist2(xvec,yvec) > r)
						continue;

					//2. w_m = 0.5 (1 + tanh[(g(y) - g(x))])
					W = 0.5f * (1.0f + (src_grad[yvec[0]][yvec[1]][2] - src_grad[xvec[0]][xvec[1]][2]));

					//3. w_d = | etf(x) DOT etf(y) |
					W *= fabs( vdot( oldetf[xvec[0]][xvec[1]], oldetf[yvec[0]][yvec[1]] ) );

					//4. phi
					W *= (vdot ( oldetf[xvec[0]][xvec[1]], oldetf[yvec[0]][yvec[1]] ) > 0.0f) ? 1.0f : -1.0f;

					//5. t_cur(y)
					vector(tcur, W*oldetf[yvec[0]][yvec[1]][0], W*oldetf[yvec[0]][yvec[1]][1], W*oldetf[yvec[0]][yvec[1]][2] );

					vsum(tsum, tcur);
				}
			}
			//normalizing tsum
			vnorm(tsum);
			vcopy(dst_etf[i][j], tsum);
		}
	}
	for (i=0; i<w; i++)		delete [] oldetf[i];	delete [] oldetf;

	float finish = clock();	
	printf("Took %.2f seconds.\n", ((finish-start)/CLOCKS_PER_SEC ));
}
void get_ETF(
	cv::Mat &src_grad, cv::Mat &dst_etf,
	int nbhd, int w, int h
) {
	printf("GetETF.. ");
	float start = clock();
	int i,j,a,b;
	V3DI xvec,yvec;
	V3DF tcur,tsum;
	float W;
	float r = (float)nbhd;

	V3DF **oldetf = new V3DF *[w];
	for (i=0; i<w; i++){
		oldetf[i] = new V3DF[h];
		for (j=0; j<h; j++) {
			cv::Vec3f &dst_etf_pixel = dst_etf.at<cv::Vec3f>(j, i);
			vcopy(oldetf[i][j], dst_etf_pixel);	
		}
	}

	for (i=0; i<w; i++){
		for (j=0; j<h; j++){
			cv::Vec3f &src_grad_pixel = src_grad.at<cv::Vec3f>(j, i);
			cv::Vec3f &dst_etf_pixel = dst_etf.at<cv::Vec3f>(j, i);
			if ( src_grad_pixel[2] == 0.0f ){
				vzero(dst_etf_pixel);
				continue;
			}

			vectori(xvec, i,j,0);
			vzero(tsum);
			for (a=-nbhd; a<=nbhd; a++){
				for (b=-nbhd; b<=nbhd; b++){
					vectori(yvec, i+a, j+b, 0);
					
					// ���ó��
					if (yvec[0]<0) yvec[0] *= (-1.0f);
					else if (yvec[0]>=w)	yvec[0] = w-(yvec[0]-(w-1));
					if (yvec[1]<0) yvec[1] *= (-1.0f);
					else if (yvec[1]>=h)	yvec[1] = h-(yvec[1]-(h-1));

					cv::Vec3f &src_grad_pixel_yvec = src_grad.at<cv::Vec3f>(yvec[1], yvec[0]);
					cv::Vec3f &src_grad_pixel_xvec = src_grad.at<cv::Vec3f>(xvec[1], xvec[0]);

					// if (src_grad[yvec[0]][yvec[1]][2] == 0.0f)	continue;
					if (src_grad_pixel_yvec[2] == 0.0f) continue;

					//1. w_s(x,y) = 0, if |xvec-yvec| >= r
					//if ( pdist2(xvec,yvec) >= r)
					if ( pdist2(xvec,yvec) > r)
						continue;

					//2. w_m = 0.5 (1 + tanh[(g(y) - g(x))])
					W = 0.5f * (1.0f + (src_grad_pixel_yvec[2] - src_grad_pixel_xvec[2]));

					//3. w_d = | etf(x) DOT etf(y) |
					W *= fabs( vdot( oldetf[xvec[0]][xvec[1]], oldetf[yvec[0]][yvec[1]] ) );

					//4. phi
					W *= (vdot ( oldetf[xvec[0]][xvec[1]], oldetf[yvec[0]][yvec[1]] ) > 0.0f) ? 1.0f : -1.0f;

					//5. t_cur(y)
					vector(tcur, W*oldetf[yvec[0]][yvec[1]][0], W*oldetf[yvec[0]][yvec[1]][1], W*oldetf[yvec[0]][yvec[1]][2] );

					vsum(tsum, tcur);
				}
			}
			//normalizing tsum
			vnorm(tsum);
			vcopy(dst_etf_pixel, tsum);
		}
	}
	for (i=0; i<w; i++)		delete [] oldetf[i];	delete [] oldetf;

	float finish = clock();	
	printf("Took %.2f seconds.\n", ((finish-start)/CLOCKS_PER_SEC ));
}

void get_flow_path(
	V3DF** src_etf, 
	V3DF** src_grad, 
	FlowPath** dst_fpath, 
	int w, int h, int threshold_S, int threshold_T
) {
	printf("SetFlowPath.. ");
	float start = clock();

	int i,j,k;

	if ( dst_fpath==NULL ){
		dst_fpath = new FlowPath *[w];
		for (i=0; i<w; i++){
			dst_fpath[i] = new FlowPath[h];
		}
	}
	else{
		for (i=0; i<w; i++)
			for (j=0; j<h; j++)
				dst_fpath[i][j].Init();
	}

	V3DF zerov = {0.0f, 0.0f, 0.0f};
	for (i=0; i<w; i++){
		for (j=0; j<h; j++){

			if ( is_similar_vector(src_etf[i][j], zerov, 0.000001f) ){
				dst_fpath[i][j].sn = dst_fpath[i][j].tn = 0;
				continue;
			}

			dst_fpath[i][j].sn = threshold_S;
			if ( dst_fpath[i][j].Alpha )	delete [] dst_fpath[i][j].Alpha;
			dst_fpath[i][j].Alpha = new V3DF[threshold_S];
			cl_set_flow_at_point_S(src_etf, &(dst_fpath[i][j]), i, j, w, h, threshold_S);

			if ( dst_fpath[i][j].Beta )		delete [] dst_fpath[i][j].Beta;
			dst_fpath[i][j].Beta = new V3DF[threshold_T];
			dst_fpath[i][j].tn = threshold_T;
			cl_set_flow_at_point_T(src_grad, &(dst_fpath[i][j]), i, j, w, h, threshold_T);
		}
	}

	printf("fpath[0][0].sn = %3d, fpath[0][0].sn = %3d\n", dst_fpath[0][0].sn, dst_fpath[0][0].tn);

	float finish = clock();	//�����ð�����
	printf("Took %.2f seconds.\n", ((finish-start)/CLOCKS_PER_SEC ));
}
void get_flow_path(
	cv::Mat &src_etf, 
	cv::Mat &src_grad, 
	FlowPath** dst_fpath, 
	int w, int h, int threshold_S, int threshold_T
)
{
	printf("SetFlowPath.. ");
	float start = clock();

	int i,j,k;

	if ( dst_fpath==NULL ){
		dst_fpath = new FlowPath *[w];
		for (i=0; i<w; i++){
			dst_fpath[i] = new FlowPath[h];
		}
	}
	else{
		for (i=0; i<w; i++)
			for (j=0; j<h; j++)
				dst_fpath[i][j].Init();
	}

	V3DF zerov = {0.0f, 0.0f, 0.0f};
	V3DF tmp = {0.0f, 0.0f, 0.0f};
	for (i=0; i<w; i++){
		for (j=0; j<h; j++){
			cv::Vec3f &src_etf_pixel = src_etf.at<cv::Vec3f>(j, i);
			vcopy(tmp, src_etf_pixel);

			if ( is_similar_vector(tmp, zerov, 0.000001f) ){
				dst_fpath[i][j].sn = dst_fpath[i][j].tn = 0;
				continue;
			}

			dst_fpath[i][j].sn = threshold_S;
			if ( dst_fpath[i][j].Alpha )	delete [] dst_fpath[i][j].Alpha;
			dst_fpath[i][j].Alpha = new V3DF[threshold_S];
			cl_set_flow_at_point_S(src_etf, &(dst_fpath[i][j]), i, j, w, h, threshold_S);

			if ( dst_fpath[i][j].Beta )		delete [] dst_fpath[i][j].Beta;
			dst_fpath[i][j].Beta = new V3DF[threshold_T];
			dst_fpath[i][j].tn = threshold_T;
			cl_set_flow_at_point_T(src_grad, &(dst_fpath[i][j]), i, j, w, h, threshold_T);
		}
	}

	float finish = clock();	
	printf("Took %.2f seconds.\n", ((finish-start)/CLOCKS_PER_SEC ));
}

void get_coherent_line(
	float** src_gray_im, V3DF** src_etf, FlowPath** src_fpath, float** dst_imCL,
 	int w, int h, float threshold_T, float CL_tanh_he_thr,
	float sigmaC, float sigmaM, float P, int iterations
) {
	printf("SetCoherentLine..");
	int i,j,k;

	float** oldim;
	oldim = new float*[w];
	for (int i=0; i<w; i++){
		oldim[i] = new float[h];
	}
	for (int i=0; i<w; i++){
		for (int j=0; j<h; j++){
			oldim[i][j] = src_gray_im[i][j];
		}
	}

	// For Hg
	int MaskSize = threshold_T;
	float *Mc = new float[MaskSize];
	float *Ms = new float[MaskSize];
	float *DogMask = new float[MaskSize];
	float sigmaS=sigmaC*1.6f;
	// Dog mask;	f(t) = G(sigma_c)(t) - P*G(sigma_s)(t)
	BuildGaussianMask(MaskSize,sigmaC,Mc);	//G(sigma_c)(t)
	BuildGaussianMask(MaskSize,sigmaS,Ms);	//G(sigma_s)(t)
	ScalarBuf(Ms,P,MaskSize);	//P*G(sigma_s)(t)
	DifferenceBuf(DogMask,Mc,Ms,MaskSize);	//f(t)

	float px,py,intensity;
	float Hg,He;
	float **clhg = new float*[w];
	float **clhe = new float*[w];

	for (i=0; i<w; i++){
		clhg[i] = new float[h];
		clhe[i] = new float[h];
	}

	for(int iter=0; iter<iterations; ++iter) {
		// For He
		for (i=0; i<w; i++){
			for (j=0; j<h; j++){
				clhg[i][j] = Hg = 0.0f;

				if ( is_zero_vector(src_etf[i][j]) ){
					continue;
				}

				// 1. Hg
				if ( src_fpath[i][j].tn < threshold_T ){
					for (k=0; k<src_fpath[i][j].tn; k++){
						px = src_fpath[i][j].Beta[k][0];	py = src_fpath[i][j].Beta[k][1];
						if ( px<0||px>=w || py<0||py>=h)	continue;
						GetValAtPoint(oldim, w, h, px,py, &intensity);
						Hg += (intensity*DogMask[k]);
					}
				}
				else{
					for (k=0; k<threshold_T; k++){
						px = src_fpath[i][j].Beta[k][0];	py = src_fpath[i][j].Beta[k][1];
						if ( px<0||px>=w || py<0||py>=h)	continue;
						GetValAtPoint(oldim, w, h, px,py, &intensity);
						Hg += (intensity*DogMask[k]);
					}
				}


				clhg[i][j] = Hg;
			}
		}

		float s, div,div2;
		div = sqrt(2*PI)*sigmaM;
		div2 = 2*sigmaM*sigmaM;
		for (i=0; i<w; i++){
			for (j=0; j<h; j++){
				clhe[i][j] = 0.0f;

				// S passing (Obama image 51.01 sec -> 12.89 sec)
				if ( is_zero_vector(src_etf[i][j]) ){
					dst_imCL[i][j] = 1.0f;
					continue;
				}

				// 2. He
				//|		He(x) = (-S~+S) G(s) Hg(cx(s)) ds
				//|	��� T�� ���� ���xDog filter
				He = 0.0f;


				for (k=0; k<src_fpath[i][j].sn; k++){
					px = src_fpath[i][j].Alpha[k][0];	py = src_fpath[i][j].Alpha[k][1];
					if ( px<0||px>=w || py<0||py>=h)	continue;
					GetValAtPoint(clhg, w, h, px,py, &Hg);

					//He += (GaussianMask[k]*Hg);
					s =	sqrt((i-px)*(i-px) + (j-py)*(j-py));	// (i,j)~(px,py)
					He += ( ( exp(-(s*s)/div2) / div )*Hg );
				}

				// 3. set line(0 or 1)
				if ( He < 0 && 1.0f+tanh(He) < CL_tanh_he_thr )	// 0.997f or 1.0f
					//if ( 1.0f+tanh(He)<srcCL_tanh_he_thr )	// 0.997f or 1.0f
					dst_imCL[i][j] = 0.0f;
				else
					dst_imCL[i][j] = 1.0f;

				// 4. set old im -> line����(�ݺ�����)
				if (dst_imCL[i][j]==0.0f)
					oldim[i][j] = 0.0f;

				clhe[i][j] = He;
			}
		}
	}

	// Normalize(clhg, MAX_X, MAX_Y);
	// Normalize(clhe, MAX_X, MAX_Y);

	// char name_1[20] = "1_hg.bmp";
	// char name_2[20] = "2_he.bmp";

	// FLOATtoBMP(name_1, MAX_X, MAX_Y, clhg);
	// FLOATtoBMP(name_2, MAX_X, MAX_Y, clhe);

	delete [] Mc;		delete [] Ms;
	delete [] DogMask;
	for (i=0; i<w; i++){
		delete [] clhg[i];
		delete [] clhe[i];
	}
	delete [] clhg;		delete [] clhe;

	printf("done\n");
}
void get_coherent_line(
	cv::Mat &src_gray_im, cv::Mat &src_etf, FlowPath** src_fpath, cv::Mat &dst_imCL,
 	int w, int h, float threshold_T, float CL_tanh_he_thr,
	float sigmaC, float sigmaM, float P, int iterations
)
{
	printf("SetCoherentLine..");
	int i,j,k;

	float** oldim;
	oldim = new float*[w];
	for (int i=0; i<w; i++){
		oldim[i] = new float[h];
	}
	for (int i=0; i<w; i++){
		for (int j=0; j<h; j++){
			oldim[i][j] = src_gray_im.at<float>(j, i);
		}
	}

	// For Hg
	int MaskSize = threshold_T;
	float *Mc = new float[MaskSize];
	float *Ms = new float[MaskSize];
	float *DogMask = new float[MaskSize];
	float sigmaS=sigmaC*1.6f;
	// Dog mask;	f(t) = G(sigma_c)(t) - P*G(sigma_s)(t)
	BuildGaussianMask(MaskSize,sigmaC,Mc);	//G(sigma_c)(t)
	BuildGaussianMask(MaskSize,sigmaS,Ms);	//G(sigma_s)(t)
	ScalarBuf(Ms,P,MaskSize);	//P*G(sigma_s)(t)
	DifferenceBuf(DogMask,Mc,Ms,MaskSize);	//f(t)

	float px,py,intensity;
	float Hg,He;
	float **clhg = new float*[w];
	float **clhe = new float*[w];

	for (i=0; i<w; i++){
		clhg[i] = new float[h];
		clhe[i] = new float[h];
	}

	for(int iter=0; iter<iterations; ++iter) {
		// For He
		for (i=0; i<w; i++){
			for (j=0; j<h; j++){
				clhg[i][j] = Hg = 0.0f;

				if ( is_zero_vector(src_etf.at<cv::Vec3f>(j, i)) ){
					continue;
				}

				// 1. Hg
				if ( src_fpath[i][j].tn < threshold_T ){
					for (k=0; k<src_fpath[i][j].tn; k++){
						px = src_fpath[i][j].Beta[k][0];	py = src_fpath[i][j].Beta[k][1];
						if ( px<0||px>=w || py<0||py>=h)	continue;
						GetValAtPoint(oldim, w, h, px,py, &intensity);
						Hg += (intensity*DogMask[k]);
					}
				}
				else{
					for (k=0; k<threshold_T; k++){
						px = src_fpath[i][j].Beta[k][0];	py = src_fpath[i][j].Beta[k][1];
						if ( px<0||px>=w || py<0||py>=h)	continue;
						GetValAtPoint(oldim, w, h, px,py, &intensity);
						Hg += (intensity*DogMask[k]);
					}
				}


				clhg[i][j] = Hg;
			}
		}

		float s, div,div2;
		div = sqrt(2*PI)*sigmaM;
		div2 = 2*sigmaM*sigmaM;
		for (i=0; i<w; i++){
			for (j=0; j<h; j++){
				clhe[i][j] = 0.0f;

				// S passing (Obama image 51.01 sec -> 12.89 sec)
				if ( is_zero_vector(src_etf.at<cv::Vec3f>(j, i)) ){
					dst_imCL.at<float>(j, i) = 1.0f;
					continue;
				}

				// 2. He
				//|		He(x) = (-S~+S) G(s) Hg(cx(s)) ds
				//|	��� T�� ���� ���xDog filter
				He = 0.0f;


				for (k=0; k<src_fpath[i][j].sn; k++){
					px = src_fpath[i][j].Alpha[k][0];	py = src_fpath[i][j].Alpha[k][1];
					if ( px<0||px>=w || py<0||py>=h)	continue;
					GetValAtPoint(clhg, w, h, px,py, &Hg);

					//He += (GaussianMask[k]*Hg);
					s =	sqrt((i-px)*(i-px) + (j-py)*(j-py));	// (i,j)~(px,py)
					He += ( ( exp(-(s*s)/div2) / div )*Hg );
				}

				// 3. set line(0 or 1)
				if ( He < 0 && 1.0f+tanh(He) < CL_tanh_he_thr )	// 0.997f or 1.0f
					//if ( 1.0f+tanh(He)<srcCL_tanh_he_thr )	// 0.997f or 1.0f
					dst_imCL.at<float>(j, i) = 0.0f;
				else
					dst_imCL.at<float>(j, i) = 1.0f;

				// 4. set old im -> line����(�ݺ�����)
				if (dst_imCL.at<float>(j, i)==0.0f)
					oldim[i][j] = 0.0f;

				clhe[i][j] = He;
			}
		}
	}

	delete [] Mc;		delete [] Ms;
	delete [] DogMask;
	for (i=0; i<w; i++){
		delete [] clhg[i];
		delete [] clhe[i];
	}
	delete [] clhg;		delete [] clhe;

	printf("done\n");
}


void cl_set_flow_at_point_S(
	V3DF** src_etf, FlowPath* dst_fpath,
	int px, int py, int w, int h, int threshold_S /*, int* nPoints, V3DF *Alpha*/
) {
	float x,y, x_,y_;	//ù��, ����
	V3DF etf_,  etf__;
	vcopy(etf__, src_etf[px][py]);

	int nN, nP, i;
	int init_n = threshold_S;
	int my_n = 2;
	V3DF **alpha = new V3DF*[my_n];
	for (i=0; i<my_n; i++)
		alpha[i] = new V3DF[init_n/2+1];

	nP=nN=0;
	x=px;	y=py;

	float angle;
	int threshold_Angle = 15.0f;
	int threshold_Angle2 = 135.0f;

	float sqrt2 = sqrt(2.0f);	//sqrt(2)��ŭ ����.
	//-S
	int limitation = 0;
	
	while ( nN < init_n/2 ){
		if (limitation++ > 1000) break;

		if ( x<0||x>=w||y<0||y>=h )	break;	//image ��

		//1. �ش� ������ etf
		V3DF_interpolate(src_etf, etf_, x, y, w, h);
		//2. normalize

		vnorm(etf_);
		if ( is_zero_vector(etf_) )	break;

		angle = angle_r(etf_,etf__)*180.0f/PI;
		if ( angle > threshold_Angle2 ){
			vnegate(etf_);
		}
		if ( angle > threshold_Angle ) {
			//������ etf�� �̷�� ���� �Ӱ�ġ �̻��̸�
			break;
		}	

		x_ = x-etf_[0];		y_ = y-etf_[1];

		if ( x==x_ && y==y_ )	break;	// ���ѷ��� ����
		if ( ((int) x != (int) x_) || ((int) y != (int) y_) )	// �ʹ� �����ϰ� ���� �ʱ� ����
			vector( alpha[NEGATIVE][nN++], x_, y_, 0.0f );
		x=x_;	y=y_;

		vcopy(etf__,etf_);
	}

	//+S
	x=px;	y=py;
	vcopy(etf__, src_etf[px][py]);
	limitation = 0;
	while ( nP<init_n/2 ){
		if (limitation++ > 1000) break;

		if ( x<0||x>=w||y<0||y>=h )	break;	//image ��

		//1. �ش� ������ etf
		V3DF_interpolate(src_etf, etf_, x, y, w, h);
		//2. normalize
		vnorm(etf_);
		if ( is_zero_vector(etf_) )	break;

		angle = angle_r(etf_,etf__)*180.0f/PI;
		if ( angle > threshold_Angle2 ){
			vnegate(etf_);
		}
		if ( angle > threshold_Angle )	//������ etf�� �̷�� ���� �Ӱ�ġ �̻��̸�
			break;

		x_ = x+etf_[0];		y_ = y+etf_[1];
		if ( x==x_ && y==y_ )	break;	// ���ѷ��� ����
		if ( ((int) x != (int) x_) || ((int) y != (int) y_) )	// �ʹ� �����ϰ� ���� �ʱ� ����
			vector( alpha[POSITIVE][nP++], x_, y_, 0.0f );
		x=x_;	y=y_;

		vcopy(etf__,etf_);
	}

	dst_fpath->sn = nN+nP+1;
	if ( dst_fpath->sn < 10 ){
		dst_fpath->sn = 0;
		for (i=0; i<2; i++)		delete [] alpha[i];	delete [] alpha;
		return;
	}

	//���� alpha(2D) -> Alpha(1D)(-S ~ +S)
	for (i=0; i<init_n; i++)	vzero( dst_fpath->Alpha[i] );	//0. �ʱ�ȭ
	for (i=0; i<nN; i++)	
		vcopy( dst_fpath->Alpha[i], alpha[NEGATIVE][nN-1-i] );	//1. negative ��
	vector(dst_fpath->Alpha[i], px,py,0.0f);	//2. �߽ɰ�
	int n=nN+1;	//���ݱ��� ���ǵ� Alpha ����
	for (i=0; i<nP; i++)	
		vcopy( dst_fpath->Alpha[n+i], alpha[POSITIVE][i] );	//3. positive ��

	for (i=0; i<my_n; i++)		delete [] alpha[i];	delete [] alpha;
}
void cl_set_flow_at_point_S(
	cv::Mat &src_etf, FlowPath* dst_fpath,
	int px, int py, int w, int h, int threshold_S /*, int* nPoints, V3DF *Alpha*/
)
{
	float x,y, x_,y_;	//ù��, ����
	V3DF etf_,  etf__;

	cv::Vec3f &src_etf_pixel = src_etf.at<cv::Vec3f>(py, px);
	vcopy(etf__, src_etf_pixel);

	int nN, nP, i;
	int init_n = threshold_S;
	int my_n = 2;
	V3DF **alpha = new V3DF*[my_n];
	for (i=0; i<my_n; i++)
		alpha[i] = new V3DF[init_n/2+1];

	nP=nN=0;
	x=px;	y=py;

	float angle;
	int threshold_Angle = 15.0f;
	int threshold_Angle2 = 135.0f;

	float sqrt2 = sqrt(2.0f);	//sqrt(2)��ŭ ����.
	//-S
	int limitation = 0;
	
	while ( nN < init_n/2 ){
		if (limitation++ > 1000) break;

		if ( x<0||x>=w||y<0||y>=h )	break;	//image ��

		//1. �ش� ������ etf
		V3DF_interpolate(src_etf, etf_, x, y, w, h);
		//2. normalize

		vnorm(etf_);
		if ( is_zero_vector(etf_) )	break;

		angle = angle_r(etf_,etf__)*180.0f/PI;
		if ( angle > threshold_Angle2 ){
			vnegate(etf_);
		}
		if ( angle > threshold_Angle ) {
			//������ etf�� �̷�� ���� �Ӱ�ġ �̻��̸�
			break;
		}	

		x_ = x-etf_[0];		y_ = y-etf_[1];

		if ( x==x_ && y==y_ )	break;	// ���ѷ��� ����
		if ( ((int) x != (int) x_) || ((int) y != (int) y_) )	// �ʹ� �����ϰ� ���� �ʱ� ����
			vector( alpha[NEGATIVE][nN++], x_, y_, 0.0f );
		x=x_;	y=y_;

		vcopy(etf__,etf_);
	}

	//+S
	x=px;	y=py;
	// vcopy(etf__, src_etf[px][py]);
	vcopy(etf__, src_etf.at<cv::Vec3f>(py, px));
	limitation = 0;
	while ( nP<init_n/2 ){
		if (limitation++ > 1000) break;

		if ( x<0||x>=w||y<0||y>=h )	break;	//image ��

		//1. �ش� ������ etf
		V3DF_interpolate(src_etf, etf_, x, y, w, h);
		//2. normalize
		vnorm(etf_);
		if ( is_zero_vector(etf_) )	break;

		angle = angle_r(etf_,etf__)*180.0f/PI;
		if ( angle > threshold_Angle2 ){
			vnegate(etf_);
		}
		if ( angle > threshold_Angle )	//������ etf�� �̷�� ���� �Ӱ�ġ �̻��̸�
			break;

		x_ = x+etf_[0];		y_ = y+etf_[1];
		if ( x==x_ && y==y_ )	break;	// ���ѷ��� ����
		if ( ((int) x != (int) x_) || ((int) y != (int) y_) )	// �ʹ� �����ϰ� ���� �ʱ� ����
			vector( alpha[POSITIVE][nP++], x_, y_, 0.0f );
		x=x_;	y=y_;

		vcopy(etf__,etf_);
	}

	dst_fpath->sn = nN+nP+1;
	if ( dst_fpath->sn < 10 ){
		dst_fpath->sn = 0;
		for (i=0; i<2; i++)		delete [] alpha[i];	delete [] alpha;
		return;
	}

	//���� alpha(2D) -> Alpha(1D)(-S ~ +S)
	for (i=0; i<init_n; i++)	vzero( dst_fpath->Alpha[i] );	//0. �ʱ�ȭ
	for (i=0; i<nN; i++)	
		vcopy( dst_fpath->Alpha[i], alpha[NEGATIVE][nN-1-i] );	//1. negative ��
	vector(dst_fpath->Alpha[i], px,py,0.0f);	//2. �߽ɰ�
	int n=nN+1;	//���ݱ��� ���ǵ� Alpha ����
	for (i=0; i<nP; i++)	
		vcopy( dst_fpath->Alpha[n+i], alpha[POSITIVE][i] );	//3. positive ��

	for (i=0; i<my_n; i++)		delete [] alpha[i];	delete [] alpha;
}

void cl_set_flow_at_point_T(
	V3DF** src_grad, FlowPath* dst_fpath,
	int px, int py, int w, int h, int threshold_T /*int *nPoints, V3DF *Beta*/
) {
	float x,y, x_,y_;
	V3DF grad_;
	int nN, nP, i;
	int init_n = threshold_T;
	V3DF **beta = new V3DF*[2];
	for (i=0; i<2; i++)
		beta[i] = new V3DF[init_n];

	nP=nN=0;
	x=px;	y=py;
	vcopy(grad_, src_grad[px][py]);
	vnorm(grad_);	//normalize

	if ( is_zero_vector(grad_) ){
		dst_fpath->tn = 1;
		vector(dst_fpath->Beta[0], px,py,0.0f);
		for (i=0; i<2; i++)
			delete [] beta[i];
		delete [] beta;
		return;
	}

	//-T
	while ( nN < (init_n/2) ){
		if ( x<0||x>=w||y<0||y>=h )	break;	
		x_ = x-grad_[0];
		y_ = y-grad_[1];
		vector( beta[NEGATIVE][nN++], x_,y_,grad_[2]);
		x=x_;	y=y_;
	}

	//+T
	x=px;	y=py;
	while ( nP < (init_n/2) ){
		if ( x<0||x>=w||y<0||y>=h )	break;	//image ��

		x_ = x+grad_[0];
		y_ = y+grad_[1];
		vector( beta[POSITIVE][nP++], x_,y_,grad_[2]);
		x=x_;	y=y_;
	}

	//���� beta(2D) -> Beta(1D)(-T ~ +T)
	for (i=0; i<init_n; i++)	vzero( dst_fpath->Beta[i] );	//0. �ʱ�ȭ
	for (i=0; i<nN; i++)	vcopy( dst_fpath->Beta[i], beta[NEGATIVE][nN-1-i] );		//1. negative ��
	vector(dst_fpath->Beta[i], px,py,0.0f);	//2. �߽ɰ�
	int n=nN+1;	//���ݱ��� ���ǵ� Beta ����
	for (i=0; i<nP; i++)	vcopy( dst_fpath->Beta[n+i], beta[POSITIVE][i] );	//3. positive ��
	dst_fpath->tn = nN+nP+1;

	for (i=0; i<2; i++)
		delete [] beta[i];
	delete [] beta;
}
void cl_set_flow_at_point_T(
	cv::Mat &src_grad, FlowPath* dst_fpath,
	int px, int py, int w, int h, int threshold_T /*int *nPoints, V3DF *Beta*/
) {
	float x,y, x_,y_;
	V3DF grad_;
	int nN, nP, i;
	int init_n = threshold_T;
	V3DF **beta = new V3DF*[2];
	for (i=0; i<2; i++)
		beta[i] = new V3DF[init_n];

	nP=nN=0;
	x=px;	y=py;
	vcopy(grad_, src_grad.at<cv::Vec3f>(py, px));
	vnorm(grad_);	//normalize

	if ( is_zero_vector(grad_) ){
		dst_fpath->tn = 1;
		vector(dst_fpath->Beta[0], px,py,0.0f);
		for (i=0; i<2; i++)
			delete [] beta[i];
		delete [] beta;
		return;
	}

	//-T
	while ( nN < (init_n/2) ){
		if ( x<0||x>=w||y<0||y>=h )	break;	
		x_ = x-grad_[0];
		y_ = y-grad_[1];
		vector( beta[NEGATIVE][nN++], x_,y_,grad_[2]);
		x=x_;	y=y_;
	}

	//+T
	x=px;	y=py;
	while ( nP < (init_n/2) ){
		if ( x<0||x>=w||y<0||y>=h )	break;	//image ��

		x_ = x+grad_[0];
		y_ = y+grad_[1];
		vector( beta[POSITIVE][nP++], x_,y_,grad_[2]);
		x=x_;	y=y_;
	}

	//���� beta(2D) -> Beta(1D)(-T ~ +T)
	for (i=0; i<init_n; i++)	vzero( dst_fpath->Beta[i] );	//0. �ʱ�ȭ
	for (i=0; i<nN; i++)	vcopy( dst_fpath->Beta[i], beta[NEGATIVE][nN-1-i] );		//1. negative ��
	vector(dst_fpath->Beta[i], px,py,0.0f);	//2. �߽ɰ�
	int n=nN+1;	//���ݱ��� ���ǵ� Beta ����
	for (i=0; i<nP; i++)	vcopy( dst_fpath->Beta[n+i], beta[POSITIVE][i] );	//3. positive ��
	dst_fpath->tn = nN+nP+1;

	for (i=0; i<2; i++)
		delete [] beta[i];
	delete [] beta;
}

void V3DF_interpolate(
	V3DF** src, V3DF dst, 
	float px, float py, int w, int h
) {
	int ulx,uly,lrx,lry;		//	upper-left (x, y), lower-right (x, y)
	V3DF etfA,etfB,etfC,etfD;	//	etf at A,B,C,D
	V3DF e1,e2,t1,t2;

	if (px+1>=w || py+1>=h){
		vector(dst, 0.0f,0.0f,0.0f);
		return;
	}

	ulx = (int)(px);
	uly = (int)(py+1);
	lrx = (int)(px+1);
	lry = (int)(py);
	vcopy(etfA, src[ulx][lry]);
	vcopy(etfB, src[lrx][lry]);
	vcopy(etfC, src[ulx][uly]);
	vcopy(etfD, src[lrx][uly]);

	//distA: A~(px,py) �Ÿ�
	//eA: etf at A
	//1. distB*etfA + distA*etfB
	vcopy(t1,etfA);	vcopy(t2,etfB);
	vscale(t1,lrx-px);	vscale(t2,px-ulx);
	vadd(e1,t1,t2);
	//2. distD*etfC + distC*etfD
	vcopy(t1,etfC);	vcopy(t2,etfD);
	vscale(t1,lrx-px);	vscale(t2,px-ulx);
	vadd(e2,t1,t2);
	//3. distC*1 + distA*2
	vcopy(t1,e1);		vcopy(t2,e2);
	vscale(t1,uly-py);	vscale(t2,py-lry);
	vadd(dst,t1,t2);

	vnorm(dst);
}
void V3DF_interpolate(
	cv::Mat &src, V3DF dst, 
	float px, float py, int w, int h
) {
	int ulx,uly,lrx,lry;		//	upper-left (x, y), lower-right (x, y)
	V3DF etfA,etfB,etfC,etfD;	//	etf at A,B,C,D
	V3DF e1,e2,t1,t2;

	if (px+1>=w || py+1>=h){
		vector(dst, 0.0f,0.0f,0.0f);
		return;
	}

	ulx = (int)(px);
	uly = (int)(py+1);
	lrx = (int)(px+1);
	lry = (int)(py);
	// vcopy(etfA, src[ulx][lry]);
	vcopy(etfA, src.at<cv::Vec3f>(lry, ulx));
	// vcopy(etfB, src[lrx][lry]);
	vcopy(etfB, src.at<cv::Vec3f>(lry, lrx));
	// vcopy(etfC, src[ulx][uly]);
	vcopy(etfC, src.at<cv::Vec3f>(uly, ulx));
	// vcopy(etfD, src[lrx][uly]);
	vcopy(etfD, src.at<cv::Vec3f>(uly, lrx));

	//distA: A~(px,py) �Ÿ�
	//eA: etf at A
	//1. distB*etfA + distA*etfB
	vcopy(t1,etfA);	vcopy(t2,etfB);
	vscale(t1,lrx-px);	vscale(t2,px-ulx);
	vadd(e1,t1,t2);
	//2. distD*etfC + distC*etfD
	vcopy(t1,etfC);	vcopy(t2,etfD);
	vscale(t1,lrx-px);	vscale(t2,px-ulx);
	vadd(e2,t1,t2);
	//3. distC*1 + distA*2
	vcopy(t1,e1);		vcopy(t2,e2);
	vscale(t1,uly-py);	vscale(t2,py-lry);
	vadd(dst,t1,t2);

	vnorm(dst);
}

// void get_vector_path(
// 	V3DF** src_etf, float** src_imCL, 
// 	FlowPath** dst_vpath, 
// 	int w, int h, int threshold_S
// ) {
// 	int i,j,k;
// 	if ( !dst_vpath ){
// 		dst_vpath = new FlowPath *[w];
// 		for (i=0; i<w; i++)
// 			dst_vpath[i] = new FlowPath[h];
// 	}
// 	else{
// 		for (i=0; i<w; i++)
// 			for (j=0; j<h; j++)
// 				dst_vpath[i][j].Init();
// 	}

// 	int n;
// 	FlowPath tmp;
// 	int N = 1000;
// 	tmp.Alpha = new V3DF[N+1];
// 	for (i=0; i<w; i++){
// 		for (j=0; j<h; j++){
// 			// 1. �ɷ�����(not coherent line OR zero etf)
// 			if ( src_imCL[i][j]!=0.0f || is_zero_vector(src_etf[i][j]) ){
// 				dst_vpath[i][j].sn = 0;		dst_vpath[i][j].tn = 0;
// 				continue;
// 			}

// 			tmp.InitShallow();

// 			// 2. �ӽ����� ������ ���� �ֱ� & ���� ���� (S ����)
// 			tmp.sn = N;
// 			vpath_set_flow_at_point_S(src_etf, src_imCL, &(tmp), i, j, w, h, threshold_S);

// 			// SetFlowAtPointS4Vectorization(i,j, &(tmp.sn), tmp.Alpha);	// flow �����ٰ� coherent line�� �ƴϸ� �����

// 			// 3. �޸� �Ҵ� �� ���� (S ����)
// 			n = dst_vpath[i][j].sn = tmp.sn;
// 			if ( n!=0 ){
// 				if ( dst_vpath[i][j].Alpha )	delete [] dst_vpath[i][j].Alpha;
// 				dst_vpath[i][j].Alpha = new V3DF[n];
// 				for (k=0; k<n; k++)
// 					vcopy(dst_vpath[i][j].Alpha[k], tmp.Alpha[k]);
// 			}
// 		}
// 	}
// }

// void vpath_set_flow_at_point_S(
// 	V3DF** src_etf, float** src_imCL,
// 	FlowPath* dst_vpath,
// 	int px, int py, int w, int h, int threshold_S
// ) {
// 	// flow �����ٰ� coherent line�� �ƴϸ� �����
// 	float x,y, x_,y_;	//ù��, ����
// 	V3DF etf_/*,etf__/*��etf*/;
// 	int nN, nP, i;
// 	int init_n = dst_vpath->sn;
// 	V3DF **alpha = new V3DF*[2];
// 	for (i=0; i<2; i++)
// 		alpha[i] = new V3DF[init_n];

// 	nP=nN=0;
// 	x=px;	y=py;

// 	float angle;
// 	int threshold_Angle = 15.0f;
// 	int threshold_Angle2 = 135.0f;
// 	V3DF etf__;	vcopy(etf__, src_etf[px][py]);

// 	float sqrt2 = sqrt(2.0f);	//sqrt(2)��ŭ ����.
// 	//-S
// 	int limitation = 0;
// 	while ( nN<init_n/2 ){
// 		if (limitation++ > 100) break;

// 		if ( x<0||x>=w||y<0||y>=h )	break;	//image ��

// 		//1. �ش� ������ etf
// 		V3DF_interpolate(src_etf, etf_, x, y, w, h);
// 		//2. normalize
// 		vnorm(etf_);
// 		if ( is_zero_vector(etf_) )	break;

// 		angle = angle_r(etf_,etf__)*180.0f/PI;
// 		if ( angle > threshold_Angle2 ){
// 			vnegate(etf_);
// 		}
// 		if ( angle > threshold_Angle )	//������ etf�� �̷�� ���� �Ӱ�ġ �̻��̸�
// 			break;

// 		x_ = x-etf_[0];
// 		y_ = y-etf_[1];
// 		if ( x_<0 || x_>=w || y_<0 || y_>=h || src_imCL[(int)x_][(int)y_]!=0.0f )	break;				// coherent line �� �ƴϸ� stop
// 		if ( ((int) x != (int) x_) || ((int) y != (int) y_) )	// �ʹ� �����ϰ� ���� �ʱ� ����
// 			vector( alpha[NEGATIVE][nN++], x_, y_, 0.0f );
// 		x=x_;	y=y_;

// 		vcopy(etf__,etf_);
// 	}

// 	//+S
// 	x=px;	y=py;
// 	vcopy(etf__, src_etf[px][py]);
// 	limitation = 0;
// 	while ( nP<init_n/2 ){
// 		if (limitation++ > 100) break;

// 		if ( x<0||x>=w||y<0||y>=h )	break;	//image ��

// 		//1. �ش� ������ etf
// 		V3DF_interpolate(src_etf, etf_, x, y, w, h);
// 		//2. normalize
// 		vnorm(etf_);
// 		if ( is_zero_vector(etf_) )	break;

// 		angle = angle_r(etf_,etf__)*180.0f/PI;
// 		if ( angle > threshold_Angle2 ){
// 			vnegate(etf_);
// 		}
// 		if ( angle > threshold_Angle )	//������ etf�� �̷�� ���� �Ӱ�ġ �̻��̸�
// 			break;

// 		x_ = x+etf_[0];
// 		y_ = y+etf_[1];
// 		if ( x_<0 || x_>=w || y_<0 || y_>=h || src_imCL[(int)x_][(int)y_]!=0.0f )	break;				// coherent line �� �ƴϸ� stop
// 		if ( ((int) x != (int) x_) || ((int) y != (int) y_) )
// 			vector( alpha[POSITIVE][nP++], x_, y_, 0.0f );
// 		x=x_;	y=y_;

// 		vcopy(etf__,etf_);
// 	}
// 	dst_vpath->sn = nN+nP+1;
// 	if ( dst_vpath->sn < 10 ){
// 		dst_vpath->sn = 0;
// 		for (i=0; i<2; i++)		delete [] alpha[i];	delete [] alpha;
// 		return;
// 	}

// 	//���� alpha(2D) -> Alpha(1D)(-S ~ +S)
// 	for (i=0; i<init_n; i++)	vzero( dst_vpath->Alpha[i] );	//0. �ʱ�ȭ
// 	for (i=0; i<nN; i++)	
// 		vcopy( dst_vpath->Alpha[i], alpha[NEGATIVE][nN-1-i] );	//1. negative ��
// 	vector(dst_vpath->Alpha[i], px,py,0.0f);	//2. �߽ɰ�
// 	int n=nN+1;	//���ݱ��� ���ǵ� Alpha ����
// 	for (i=0; i<nP; i++)	
// 		vcopy( dst_vpath->Alpha[n+i], alpha[POSITIVE][i] );	//3. positive ��

// 	for (i=0; i<2; i++)		delete [] alpha[i];	delete [] alpha;
// }

// void vectorize_vector_path(
// 	FlowPath** src_vpath,
// 	int** dst_skltn,
// 	int w, int h
// ) {
// 	int i,j;
// 	if ( !dst_skltn ){
// 		dst_skltn = new int *[w];
// 		for (i=0; i<w; i++){
// 			dst_skltn[i] = new int[h];
// 			for (j=0; j<h; j++)
// 				dst_skltn[i][j] = 0;
// 		}
// 	}
// 	else{
// 		for (i=0; i<w; i++)
// 			for (j=0; j<h; j++)
// 				dst_skltn[i][j] = 0;
// 	}

// 	_vectorize_build_skltn(src_vpath, dst_skltn, w, h);			//	��� pixel���� skeleton�� ��� -->	0: outside
// 	//										1: border
// 	//										2: border���� �� �ȼ� ��������
// 	//										3: border���� �� �ȼ� ��������
// 	_vectorize_compute_center_line(src_vpath, dst_skltn, w, h);	//	Streams[i][j]���� stream ���
// 	//	stream�� ���� pixel�� index (i, j)�� slist�� ����
// 	// if ( !dparam.StreamsExisting ){	// Stream�� �������� ������
// 	// 	return;
// 	// }
// 	build_center_line1();	//	stream�� �߿��� max stream ���
// 	//	max stream�� plist�� ����
// 	build_center_line2();	//	�̿��� stream�� ã�Ƽ� ������
// 	//	�̿��� stream�� ������ stream���� qlist�� ����

// 	sort_length();			//	qlist�� stream���� ���̼����� ����
// 	//	���̼����� ������ stream���� rlist�� ����
// 	curv_line ( );			//	compute the curvatures at each point

// 	SetAlphaMax();	// [����߰�] vectorlized line(vpath) �̵���Ű�� ���� Alpha -> AlphaMax
// 	// [��� �߰�] importance ���

// 	if ( dparam.is_shifted_line_selection ){
// 		// 1. ���� ���� AlphaMax -> Alpha
// 		SwapAlphaNAlphaMax();

// 		// 2. ��ȣ��: sort_length(), curv_line()
// 		sort_length();	curv_line ( );	// �׳� �ҷ��� �ɱ�?
// 	}
// }

// void _vectorize_build_skltn(
// 	FlowPath** src_vpath,
// 	int** dst_skltn,
// 	int w, int h
// ) {
// 	int i, j, k, l, px,py;

// 	//	1. �ʱ�ȭ: Streams�� �������� pixel�� 1, �׷��� ���� pixel�� 0���� setting
// 	for ( i = 0; i < w; i++ )  {
// 		for ( j = 0; j < h; j++ ) {
// 			if ( src_vpath[i][j].sn==0 )	continue;

// 			vector(src_vpath[i][j].BBox[0], 10000.0f,  10000.0f,  10000.0f);
// 			vector(src_vpath[i][j].BBox[1], -10000.0f,  -10000.0f,  -10000.0f);
// 			for ( k = 0; k < src_vpath[i][j].sn; k++ ) {
// 				px = (int) src_vpath[i][j].Alpha[k][0];	if (px<0) px = 0;	else if (px>=w) px = w-1;
// 				py = (int) src_vpath[i][j].Alpha[k][1];	if (py<0) py = 0;	else if (py>=h) py = h-1;
// 				dst_skltn[px][py] = 1;
// 				if ( src_vpath[i][j].BBox[0][0] > px )	src_vpath[i][j].BBox[0][0] = px;
// 				if ( src_vpath[i][j].BBox[0][1] > py )	src_vpath[i][j].BBox[0][1] = py;

// 				if ( src_vpath[i][j].BBox[1][0] < px )	src_vpath[i][j].BBox[1][0] = px;
// 				if ( src_vpath[i][j].BBox[1][1] < py )	src_vpath[i][j].BBox[1][1] = py;
// 			}
// 		}
// 	}

// 	int cnt, scnt;
// 	int flag = 1;
// 	int limitation = 0;
// 	do {
// 		for ( i = 1, cnt = 0; i < w-1; i++ ) {
// 			for ( j = 1; j < h-1; j++ ) {
// 				if ( dst_skltn[i][j] != flag )
// 					continue;
// 				for ( k = i-1, scnt = 0; k <= i+1; k++ ) {
// 					for ( l = j-1; l <= j+1; l++ ) {
// 						scnt += (dst_skltn[k][l] >= flag);
// 					}
// 				}
// 				if ( scnt == 9 ) {
// 					dst_skltn[i][j] = flag+1;
// 					cnt++;
// 				}
// 			}
// 		}
// 		flag++;

// 		if (limitation++ > 100) break;
// 	} while ( cnt > 0 );

// 	printf("Maximum skeleton: %d\n", flag);

// 	for ( i = 0; i < w; i++ ) {
// 		for ( j = 0; j < h; j++ ) {
// 			src_vpath[i][j].skltn = 0.0f;
// 			if ( src_vpath[i][j].sn == 0 )
// 				continue;

// 			for ( k = 0; k < src_vpath[i][j].sn; k++ ) {
// 				px = (int) src_vpath[i][j].Alpha[k][0];	if (px<0) px = 0;	else if (px>=w) px = w-1;
// 				py = (int) src_vpath[i][j].Alpha[k][1];	if (py<0) py = 0;	else if (py>=h) py = h-1;
// 				src_vpath[i][j].skltn += dst_skltn[px][py];
// 			}
// 			src_vpath[i][j].skltn /= src_vpath[i][j].sn;
// 		}
// 	}
// }

// void _vectorize_compute_center_line(
// 	FlowPath** src_vpath,
// 	int* dst_ns_list, int* dst_slist, 
// 	int w, int h
// ) {
// 	int i, j;
// 	float max_cl;
// 	int cl_cnt;

// 	for ( i = 0, max_cl = 0.0f, cl_cnt = 0; i < w; i++ ) {
// 		for ( j = 0; j < h; j++ ) {
// 			if ( src_vpath[i][j].sn == 0 )
// 				continue;

// 			//	1. estimate center_line
// 			compute_center_line ( i, j );
// 			cl_cnt++;
// 			if ( src_vpath[i][j].cl > max_cl )
// 				max_cl = src_vpath[i][j].cl;
// 		}
// 	}
// 	dst_slist = cl_cnt;

// 	if ( dst_slist )	delete [] dst_slist;
// 	dst_slist = new V3DI[dst_ns_list];
// 	for ( i = 0, cl_cnt = 0; i < w; i++ ) {
// 		for ( j = 0; j < h; j++ ) {
// 			if ( src_vpath[i][j].cl > 0 ) {
// 				slist[cl_cnt][0] = i;
// 				slist[cl_cnt][1] = j;
// 				slist[cl_cnt++][2] = (int) (vpath[i][j].cl * 10000.0f);
// 			}
// 		}
// 	}

// 	if ( n_slist==0 ){	// There is no Stream
// 		nVpath = n_qlist = n_plist = n_slist;
// 		dparam.StreamsExisting = false;
// 		return;
// 	}

// 	qs ( 0, n_slist-1, slist );
// }

// void __vectorize_compute_center_line_at_point(
// 	int** src_skltn,
// 	FlowPath** dst_vpath,
// 	int px, int py, int w, int h
// ) {
// 	int i, j,r,c;
// 	V3DF vec;
// 	V3DF nvec;
// 	V3DF Zaxis;
// 	float tcl;
// 	V3DF tpt1, tpt2;
// 	V3DF weights = { 0.5f, 0.3f, 0.2f };

// 	vector ( Zaxis, 0.0f, 0.0f, 1.0f );
// 	dst_vpath[px][py].cl = 0;
// 	for ( i = 0; i < dst_vpath[px][py].sn; i++ ) {
// 		if ( i == 0 ) {
// 			nvector ( vec, dst_vpath[px][py].Alpha[i+1], dst_vpath[px][py].Alpha[i] );
// 		}
// 		else if ( i == dst_vpath[px][py].sn - 1 ) {
// 			nvector ( vec, dst_vpath[px][py].Alpha[i], dst_vpath[px][py].Alpha[i-1] );
// 		}
// 		else {
// 			nvector ( vec, dst_vpath[px][py].Alpha[i+1], dst_vpath[px][py].Alpha[i-1] );
// 		}

// 		nvcross ( nvec, vec, Zaxis );

// 		r = (int) dst_vpath[px][py].Alpha[i][0];	if (r<0) r = 0;	else if (r>=w) r = w-1;
// 		c = (int) dst_vpath[px][py].Alpha[i][1];	if (c<0) c = 0;	else if (c>=h) c = h-1;
// 		tcl = src_skltn[r][c];

// 		vcopy( tpt1, dst_vpath[px][py].Alpha[i] );
// 		for ( j = 0; j < 3; j++ ) {
// 			get_point ( tpt2, tpt1, 1.0f, nvec );
// 			if ( tpt2[0] < 0.0f )
// 				tpt2[0] = 0.0f;
// 			if ( tpt2[1] < 0.0f )
// 				tpt2[1] = 0.0f;
// 			if ( tpt2[0] > w-1 )
// 				tpt2[0] = w-1;
// 			if ( tpt2[1] > h-1 )
// 				tpt2[1] = h-1;
// 			tcl += weights[j] * src_skltn[(int) tpt2[0]][(int) tpt2[1]];
// 			vcopy( tpt1, tpt2 );
// 		}

// 		vnegate ( nvec );
// 		vcopy( tpt1, dst_vpath[px][py].Alpha[i] );
// 		for ( j = 0; j < 3; j++ ) {
// 			get_point ( tpt2, tpt1, 1.0f, nvec );
// 			if ( tpt2[0] < 0.0f )
// 				tpt2[0] = 0.0f;
// 			if ( tpt2[1] < 0.0f )
// 				tpt2[1] = 0.0f;
// 			if ( tpt2[0] > w-1 )
// 				tpt2[0] = w-1;
// 			if ( tpt2[1] > h-1 )
// 				tpt2[1] = h-1;
// 			tcl += weights[j] * src_skltn[(int) tpt2[0]][(int) tpt2[1]];
// 			vcopy( tpt1, tpt2 );
// 		}

// 		dst_vpath[px][py].cl += tcl;
// 	}
// }

