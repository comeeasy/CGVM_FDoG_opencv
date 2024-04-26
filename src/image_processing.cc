#include "image_processing.h"


void apply_FBL_filter( 
    V3DF** src_cim, 
	FlowPath** src_fpath,
    V3DF** cimFBL,
    int w, int h,
    float sigma_e, float gamma_e, float sigma_g, float gamma_g,	int threshold_T,
	int iteration
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
    cv::Mat &cimFBL,
    int w, int h,
    float sigma_e, float gamma_e, 
	float sigma_g, float gamma_g, 
	int threshold_T,
	int iteration
)
{
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
			vcopy(InImg[i][j], src_cim.at<cv::Vec3f>(j,i));
			vcopy(OImg[i][j], src_cim.at<cv::Vec3f>(j,i));
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
			vcopy(cimFBL.at<cv::Vec3f>(j,i), OImg[i][j]);

	for (i=0; i<w; i++){
		delete [] InImg[i];
		delete [] OImg[i];
	}
	delete [] InImg;
	delete [] OImg;

	printf("done\n");
}
cv::Mat apply_FBL_filter(
	cv::Mat &src_cim, 
	FlowPath** src_fpath,
    float sigma_e, float gamma_e, 
	float sigma_g, float gamma_g,
	int threshold_T,
	int iteration
)
{
	int w = src_cim.cols, h = src_cim.rows;

	cv::Mat cimFBL = cv::Mat::zeros(src_cim.size(), CV_32FC3);
	apply_FBL_filter(
		src_cim, src_fpath, cimFBL, w, h, 
		sigma_e, gamma_e, sigma_g, gamma_g, 
		threshold_T, iteration
	);
	return cimFBL;
}


bool IsSameCluster(int sr, int sc, int rr, int rc, int m, int n, int hs, int hr, V3DF **img)
{
	int hr2 = hr*hr;
	int hs2 = hs*hs;
	int dx, dy;		//공간차이
	int dL, du, dv;	//색차이
	dx=sr-m;	dy=sc-n;
	if (dx*dx + dy*dy <= hs2){
		dL = img[rr][rc][0] - img[m][n][0];
		du = img[rr][rc][1] - img[m][n][1];
		dv = img[rr][rc][2] - img[m][n][2];
		if (dL*dL + du*du + dv*dv <= hr2){
			return true;
		}
	}

	return false;
}

int push(short *stackx, short *stacky, short vx, short vy, int *top)
{
	if (*top >= 500000) return (-1);		// 숫자가 너무 작으면, 직사각형으로 세그멘테이션 된다.	//[요기]
	(*top)++;
	stackx[*top] = vx;
	stacky[*top] = vy;
	return (1);
}

int pop(short *stackx, short *stacky, short *vx, short *vy, int *top)
{
	if (*top == 0) return (-1);
	*vx = stackx[*top];
	*vy = stacky[*top];
	(*top)--;
	return (1);
}

int m_BlobColoring(V3DF **img, int width, int height, int hs, int hr, int **idx){
	///1. 라벨링하기
	int i,j,m,n,top, BlobArea[1000];
	short r,c, area;
	int curColor=0;

	//스택으로 사용할 메모리 할당
	short* stackx=new short [width*height];
	short* stacky=new short [width*height];

	//라벨링된 픽셀을 저장하기 위해 메모리 할당
	int *coloring = new int[width*height];
	for (i = 0; i < width*height; i++)
		coloring[i]=0;		//메모리 초기화

	for (i = 0; i < width; i++){
		for (j = 0; j < height; j++){
			//이미 방문한 점이면 처리 안함
			if (coloring[i*height+j] != 0)	continue;

			r=i;	//y축 중심위치
			c=j;	//x축 중심위치
			top = 0;	area = 1;
			curColor++;

			while (1){
GRASSFIRE:
				for (m = r-1; m <= r+1; m++){
					for (n = c-1; n <= c+1; n++){
						//4근방화소
						//if ( (m==-1 && n==-1) || (m==-1 && n==1) || (m==1 && n==-1) || (m==1 && n==1) )	continue;

						//관심 픽셀이 영상경계를 벗어나면 처리 안함
						if (m<0 || m>=width || n<0 || n>=height) continue;

						//윈도우 크기 내에 있고, 방문하지 않은 점
						//2. 클러스터링
						if ( IsSameCluster(r,c,i,j,m,n, hs, hr, img) && coloring[m*height+n]==0 ){
							coloring[m*height+n] = curColor;		//현재 라벨로 마크
							if(push(stackx,stacky,(short)m,(short)n,&top)==-1) continue;

							r=m; 
							c=n;
							area++;
							goto GRASSFIRE;
						}
					}
				}
				if (pop(stackx,stacky,&r,&c,&top) == -1) break;
			}
			if(curColor<1000) BlobArea[curColor] = area;
		}
	}

	//coloring --> idx
	for (i = 0; i < width; i++){
		for (j = 0; j < height; j++){
			idx[i][j] = coloring[i*height+j]-1;
		}
	}

	delete [] coloring; delete [] stackx; delete [] stacky;

	return curColor;
}

void InitLabel2(V3DF **img, int w, int h, int n/*regions개수*/, LABEL *label, int **idx){
	int i,j,r,c;
	int indx;	//인덱스
	//라벨 초기화(id, 개수, 대표색, 대표위치)
	for (i=0; i<n; i++){
		label[i].num_pixel = 0;		//픽셀개수
		vector(label[i].clr, 0.0f, 0.0f, 0.0f);	//색
		vector2(label[i].pos,0,0);	//위치
	}

	for (i=0; i<w; i++){
		for (j=0; j<h; j++){
			indx = idx[i][j];
			label[indx].num_pixel++;			//픽셀개수
			vadd(label[indx].clr, img[i][j]);	//색
			label[indx].pos[X] += i;		label[indx].pos[Y] +=j; //위치
		}
	}

	for (i=0; i<n; i++){

		vscalar( label[i].clr, 1.0f/(float)(label[i].num_pixel) );	//평균색
		label[i].pos[X] = (float)(label[i].pos[X]) / (float)(label[i].num_pixel);
		label[i].pos[Y] = (float)(label[i].pos[Y]) / (float)(label[i].num_pixel);
		//vscalar( (V3DF)label[i].pos, 1.0f/(float)(label[i].num_pixel) );	//무게중심
		//대표위치 정하기: 무게중심에서 
		//좌->우
		for (r=0; r<w; r++){
			if (idx[r][label[i].pos[Y]] == i){
				vector2(label[i].pos, r,label[i].pos[Y]);
				break;
			}
		}
		//상->하
		for (c=0; c<h; c++){
			if (idx[label[i].pos[X]][c] == i){
				vector2(label[i].pos, label[i].pos[X],c);
				break;
			}
		}
	}
}

int delete_smallRegion3(V3DF **img, int w, int h, int rn/*regions개수*/, LABEL *label, int M, int **idx){
	//1. 라벨링하기
	int i,j,m,n,top, BlobArea[2000];
	short r,c, area;
	int curColor=0;
	int orig_n = rn;	// 원래 영역 개수
	int id, otherid, changedid;
	V3DF thisclr, diffclr;
	float diffmin, diff;	//색차

	//스택으로 사용할 메모리 할당
	short* stackx=new short [w*h];
	short* stacky=new short [w*h];

	//라벨링된 픽셀을 저장하기 위해 메모리 할당
	int *coloring = new int[w*h];
	for (i = 0; i < w*h; i++)
		coloring[i]=0;		//메모리 초기화

	for (i = 0; i < w; i++){
		for (j = 0; j < h; j++){
			// M 이상이면 처리 안함
			id = idx[i][j];
			vcopy(thisclr, label[id].clr);
			diffmin = INFINITY;
			changedid = -1;

			// 이미 방문한 점이면 처리 안함
			if (coloring[i*h+j] != 0)	continue;

			r=i;	//y축 중심위치
			c=j;	//x축 중심위치
			top = 0;	area = 1;
			curColor++;

			while (1){
GRASSFIRE:
				for (m = r-1; m <= r+1; m++){
					for (n = c-1; n <= c+1; n++){
						//4근방화소
						//if ( (m==-1 && n==-1) || (m==-1 && n==1) || (m==1 && n==-1) || (m==1 && n==1) )	continue;

						//관심 픽셀이 영상경계를 벗어나면 처리 안함
						if (m<0 || m>=w || n<0 || n>=h) continue;

						// 같은 인덱스이고, 방문하지 않은 점
						// 클러스터링
						otherid = idx[m][n];

						//if ( id==2465 )
						//	printf("hihi\n");
						//if ( m==228&&n==247 )
						//	printf("hoho\n");

						if ( id==otherid && coloring[m*h+n]==0 ){
							coloring[m*h+n] = curColor;		//현재 라벨로 마크
							if(push(stackx,stacky,(short)m,(short)n,&top)==-1) continue;

							r=m; 
							c=n;
							area++;
							goto GRASSFIRE;
						}
						//1. 주위 픽셀 중 Luv 색 차이 최소인 라벨로 검색
						else if ( id!=otherid && otherid!=changedid ){
							vsub(diffclr, thisclr, label[otherid].clr);
							diff = diffclr[0]*diffclr[0] + diffclr[1]*diffclr[1] + diffclr[2]*diffclr[2];
							if ( diffmin>diff ){
								diffmin = diff;
								changedid = otherid;
							}
						}
					}
				}
				if (pop(stackx,stacky,&r,&c,&top) == -1) break;
			}

			if ( label[id].num_pixel < M && changedid>=0){
				// 2. 최소차이 라벨로 대체
				for (m=0; m<w; m++){
					for (n=0; n<h; n++){
						if ( idx[m][n]==id ){
							if ( changedid<0 )	continue;
							idx[m][n] = changedid;
							label[changedid].num_pixel++;
						}
					}
				}
				//3. 해당 라벨 id 제거, 개수 0으로
				label[id].num_pixel = 0;
				rn--;		//총 라벨 영역 중 M개 픽셀 이하인 영역 제외
			}

			if(curColor<2000) BlobArea[curColor] = area;
		}
	}

	// 빈 곳 정리
	int *pt = new int[orig_n];	// 빈 곳
	// 1) 빈 곳 주소 찾기
	for (i=0,j=0; i<orig_n; i++){
		if ( label[i].num_pixel==0 )
			pt[j++] = i; 
	}
	// 2) 뒤의 것 --> 앞으로
	for (i=0,j=0; i<orig_n; i++){
		if ( label[i].num_pixel==0 || i<rn )	continue;

		// label[j] <-- label[i];
		for (r=0; r<w; r++){
			for (c=0; c<h; c++){
				if (idx[r][c]==i){
					idx[r][c] = pt[j];
				}
			}
		}
		//label[pt[j]].num_pixel = label[i].num_pixel;
		//vassign(label[pt[j]].clr, label[i].clr);
		//vassign(label[pt[j]].pos, label[i].pos);
		//j++;

		label[pt[j++]] = label[i];
		label[i].num_pixel = 0;
	}

	delete [] pt;
	delete [] coloring; delete [] stackx; delete [] stacky;

	return rn;
}

void get_segmentation (
	cv::Mat &InImg, 
	int** dst_sgid, cv::Mat &dst_cimSG,
	int hs, int hr, int M, int w, int h
)
{ 
	printf("Segmentation..");

	int i,j,n;
	LABEL *label;
	V3DF **img = new V3DF *[w];
	int **idx = new int*[w];
	for (i=0; i<w; i++){
		img[i] = new V3DF[h];
		idx[i] = new int[h];
	}

	// 1. mean-shift filtering
	//MeanShiftFilter(InImg, w,h, hs, hr, img);
	init_buf(img, InImg, w,h);
	conv_RGBtoLUV(img,w,h);

	// 2. labeling
	n = m_BlobColoring(img, w,h, hs,hr, idx);
	printf ("%d 개-> ", n);

	// 3. init label
	label = new LABEL[n];
	InitLabel2(img, w,h, n, label, idx);

	// 4. M개 이하 픽셀 갖는 클러스터 삭제
	//n = delete_smallRegion2(img, w, h, n, label, M, idx);
	n = delete_smallRegion3(img, w, h, n, label, M, idx);
	InitLabel2(img, w,h, n, label, idx);
	printf ("%d 개", n);

	// 5. LUV -> RGB
	for (i=0; i<n; i++){
		LUVtoRGB(label[i].clr[0],label[i].clr[1],label[i].clr[2],
			&label[i].clr[0],&label[i].clr[1],&label[i].clr[2]);
	}

	int id;
	for (i=0; i<w; i++){
		for (j=0; j<h; j++){
			id = idx[i][j];
			dst_sgid[i][j] = id;
			vcopy(dst_cimSG.at<cv::Vec3f>(j, i), label[id].clr);
		}
	}

	for (i=0; i<w; i++){
		delete [] img[i];
		delete [] idx[i];
	}
	delete [] img;	delete [] idx;
	delete [] label;

	printf("done\n");
	// return n;
}



void get_segmentation (
	V3DF** InImg, 
	int** dst_sgid, V3DF** dst_cimSG,
	int hs, int hr, int M, int w, int h
)
{ 
	printf("Segmentation..");

	int i,j,n;
	LABEL *label;
	V3DF **img = new V3DF *[w];
	int **idx = new int*[w];
	for (i=0; i<w; i++){
		img[i] = new V3DF[h];
		idx[i] = new int[h];
	}

	// 1. mean-shift filtering
	//MeanShiftFilter(InImg, w,h, hs, hr, img);
	init_buf(img, InImg, w,h);
	conv_RGBtoLUV(img,w,h);

	// 2. labeling
	n = m_BlobColoring(img, w,h, hs,hr, idx);
	printf ("%d 개-> ", n);

	// 3. init label
	label = new LABEL[n];
	InitLabel2(img, w,h, n, label, idx);

	// 4. M개 이하 픽셀 갖는 클러스터 삭제
	//n = delete_smallRegion2(img, w, h, n, label, M, idx);
	n = delete_smallRegion3(img, w, h, n, label, M, idx);
	InitLabel2(img, w,h, n, label, idx);
	printf ("%d 개", n);

	// 5. LUV -> RGB
	for (i=0; i<n; i++){
		LUVtoRGB(label[i].clr[0],label[i].clr[1],label[i].clr[2],
			&label[i].clr[0],&label[i].clr[1],&label[i].clr[2]);
	}

	int id;
	for (i=0; i<w; i++){
		for (j=0; j<h; j++){
			id = idx[i][j];
			dst_sgid[i][j] = id;
			vcopy(dst_cimSG[i][j], label[id].clr);
		}
	}

	for (i=0; i<w; i++){
		delete [] img[i];
		delete [] idx[i];
	}
	delete [] img;	delete [] idx;
	delete [] label;

	printf("done\n");
	// return n;
}

void sum_FDoG_FBL(
	cv::Mat &src_imCL, cv::Mat &src_cimSG,
	cv::Mat &dst_cimsum,
	int w, int h
){
	int i,j;

	for (i=0; i<w; i++){
		for (j=0; j<h; j++){
			if (src_imCL.at<float>(j, i) == 0.0f) {
				vzero(dst_cimsum.at<cv::Vec3f>(j, i));
			} else {
				// vcopy(cimsum[i][j], cimSG[i][j]);
				dst_cimsum.at<cv::Vec3f>(j, i) = src_cimSG.at<cv::Vec3f>(j, i);
			}
		}
	}
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
void get_gradient(const cv::Mat& src, cv::Mat& dst, float grad_thr) {
    dst.create(src.size(), CV_32FC3); // make sure it's the correct type
    GradientOperator op(src, dst, grad_thr);
    cv::parallel_for_(cv::Range(0, src.rows), op);
}
cv::Mat get_gradient(
	cv::Mat &src_gray_img,
	float grad_thr
)
{
	int w = src_gray_img.cols, h = src_gray_img.rows;
	cv::Mat grad = cv::Mat::zeros(h, w, CV_32FC3);
	get_gradient(src_gray_img, grad, w, h, grad_thr);
	return grad;
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
cv::Mat get_tangent(
	cv::Mat &src_grad
)
{
	int w = src_grad.cols, h = src_grad.rows;
	cv::Mat tangent = cv::Mat::zeros(h, w, CV_32FC3);
	get_tangent(src_grad, tangent, w, h);
	return tangent;
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
void __get_ETF_partial(int start, int end, int w, int h, int nbhd, float r, cv::Mat& src_grad, cv::Mat& dst_etf, V3DF** oldetf) {
	V3DI xvec,yvec;
	V3DF tcur,tsum;
	float W;

    for (int i = start; i < end; i++) {
        for (int j = 0; j < h; j++) {
            cv::Vec3f &src_grad_pixel = src_grad.at<cv::Vec3f>(j, i);
            cv::Vec3f &dst_etf_pixel = dst_etf.at<cv::Vec3f>(j, i);
            if (src_grad_pixel[2] == 0.0f) {
                vzero(dst_etf_pixel);
                continue;
            }

            vectori(xvec, i,j,0);
			vzero(tsum);
            for (int a = -nbhd; a <= nbhd; a++) {
                for (int b = -nbhd; b <= nbhd; b++) {
                    vectori(yvec, i+a, j+b, 0);
                    // Boundary checks and corrections for yvec
                    if (yvec[0] < 0) yvec[0] *= -1;
                    else if (yvec[0] >= w) yvec[0] = w - (yvec[0] - (w - 1));
                    if (yvec[1] < 0) yvec[1] *= -1;
                    else if (yvec[1] >= h) yvec[1] = h - (yvec[1] - (h - 1));

                    cv::Vec3f &src_grad_pixel_yvec = src_grad.at<cv::Vec3f>(yvec[1], yvec[0]);
                    cv::Vec3f &src_grad_pixel_xvec = src_grad.at<cv::Vec3f>(xvec[1], xvec[0]);

                    if (src_grad_pixel_yvec[2] == 0.0f) continue;
                    if (pdist2(xvec, yvec) > r) continue;

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
            vnorm(tsum);
            vcopy(dst_etf_pixel, tsum);
        }
    }
}
void get_ETF(
	cv::Mat &src_grad, cv::Mat &dst_etf,
	int nbhd, int w, int h, int num_workers
)
{
	printf("GetETF.. ");
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

    std::vector<std::thread> threads(num_workers);
    int slice = w / num_workers;

	for (int t = 0; t < num_workers; t++) {
        int start = t * slice;
        int end = (t == num_workers - 1) ? w : (t + 1) * slice;
        threads[t] = std::thread(__get_ETF_partial, start, end, w, h, nbhd, r, std::ref(src_grad), std::ref(dst_etf), std::ref(oldetf));
    }

    for (auto &thread : threads) {
        thread.join();
    }

	for (i=0; i<w; i++)		delete [] oldetf[i];	delete [] oldetf;
}
cv::Mat get_ETF(
	cv::Mat &src_grad, cv::Mat &src_tangent,
	int nbhd, int iteration
)
{
	if (src_grad.cols != src_tangent.cols || src_grad.rows != src_tangent.rows) {
		printf("src_grad and src_tangent must have same shape. (width and height)\n");
		exit(1);
	}

	int w = src_grad.cols, h = src_tangent.rows;
	cv::Mat etf = src_tangent.clone();
	for(int i=0; i<iteration; ++i) {
		get_ETF(src_grad, etf, nbhd, w, h);
	}
	return etf;
}
cv::Mat get_ETF(
	cv::Mat &src_grad, cv::Mat &src_tangent,
	int nbhd, int iteration, size_t num_workers
)
{
	if (src_grad.cols != src_tangent.cols || src_grad.rows != src_tangent.rows) {
		printf("src_grad and src_tangent must have same shape. (width and height)\n");
		exit(1);
	}

	int w = src_grad.cols, h = src_tangent.rows;
	cv::Mat etf = src_tangent.clone();
	for(int i=0; i<iteration; ++i) {
		float start = clock();
		get_ETF(src_grad, etf, nbhd, w, h, num_workers);
		float finish = clock();	
		printf("Took %.2f seconds.\n", ((finish-start)/CLOCKS_PER_SEC ));
	}
	
	return etf;
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
void __get_flow_path_partial(int start, int end, int h, int threshold_S, int threshold_T, cv::Mat &src_etf, cv::Mat &src_grad, FlowPath** dst_fpath, int num_workers) {
    V3DF zerov = {0.0f, 0.0f, 0.0f};
    V3DF tmp = {0.0f, 0.0f, 0.0f};

    for (int i = start; i < end; i++) {
        for (int j = 0; j < h; j++) {
            cv::Vec3f &src_etf_pixel = src_etf.at<cv::Vec3f>(j, i);
            vcopy(tmp, src_etf_pixel);

            if (is_similar_vector(tmp, zerov, 0.000001f)) {
                dst_fpath[i][j].sn = dst_fpath[i][j].tn = 0;
                continue;
            }

            dst_fpath[i][j].sn = threshold_S;
            if (dst_fpath[i][j].Alpha) delete[] dst_fpath[i][j].Alpha;
            dst_fpath[i][j].Alpha = new V3DF[threshold_S];
            cl_set_flow_at_point_S(src_etf, &(dst_fpath[i][j]), i, j, src_etf.cols, src_etf.rows, threshold_S);

            dst_fpath[i][j].tn = threshold_T;
            if (dst_fpath[i][j].Beta) delete[] dst_fpath[i][j].Beta;
            dst_fpath[i][j].Beta = new V3DF[threshold_T];
            cl_set_flow_at_point_T(src_grad, &(dst_fpath[i][j]), i, j, src_grad.cols, src_grad.rows, threshold_T);
        }
    }
}
void get_flow_path(
	cv::Mat &src_etf, 
	cv::Mat &src_grad, 
	FlowPath** dst_fpath, 
	int w, int h, int threshold_S, int threshold_T,
	int num_workers
)
{
	printf("SetFlowPath.. ");
    float start = clock();

    if (dst_fpath == nullptr) {
        dst_fpath = new FlowPath *[w];
        for (int i = 0; i < w; i++) {
            dst_fpath[i] = new FlowPath[h];
        }
    } else {
        for (int i = 0; i < w; i++)
            for (int j = 0; j < h; j++)
                dst_fpath[i][j].Init();
    }

    std::vector<std::thread> threads(num_workers);
    int slice = w / num_workers;

    for (int t = 0; t < num_workers; t++) {
        int start = t * slice;
        int end = (t == num_workers - 1) ? w : (t + 1) * slice;
        threads[t] = std::thread(__get_flow_path_partial, start, end, h, threshold_S, threshold_T, std::ref(src_etf), std::ref(src_grad), dst_fpath, num_workers);
    }

    for (auto &thread : threads) {
        thread.join();
    }

    float finish = clock();
    printf("Took %.2f seconds.\n", ((finish - start) / CLOCKS_PER_SEC));
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
void __get_coherent_line_partial(int start, int end, int h, float** oldim, FlowPath** src_fpath, float *DogMask, float sigmaM, float threshold_T, float CL_tanh_he_thr, cv::Mat &dst_imCL) {
    const float div = sqrt(2 * M_PI) * sigmaM;
    const float div2 = 2 * sigmaM * sigmaM;

    for (int i = start; i < end; i++) {
        for (int j = 0; j < h; j++) {
            float Hg = 0.0f, He = 0.0f;

            // Calculate Hg
            for (int k = 0; k < std::min((int)threshold_T, src_fpath[i][j].tn); k++) {
                float px = src_fpath[i][j].Beta[k][0];
                float py = src_fpath[i][j].Beta[k][1];
                if (px < 0 || px >= end || py < 0 || py >= h) continue;
                float intensity;
                GetValAtPoint(oldim, end, h, px, py, &intensity);
                Hg += (intensity * DogMask[k]);
            }

            // Calculate He
            for (int k = 0; k < src_fpath[i][j].sn; k++) {
                float px = src_fpath[i][j].Alpha[k][0];
                float py = src_fpath[i][j].Alpha[k][1];
                if (px < 0 || px >= end || py < 0 || py >= h) continue;
                float dist = sqrt((i - px) * (i - px) + (j - py) * (j - py));
                He += (exp(-(dist * dist) / div2) / div) * Hg;
            }

            // Set line
            if (He < 0 && 1.0f + tanh(He) < CL_tanh_he_thr)
                dst_imCL.at<float>(j, i) = 0.0f;
            else
                dst_imCL.at<float>(j, i) = 1.0f;

            // Update old image
            if (dst_imCL.at<float>(j, i) == 0.0f)
                oldim[i][j] = 0.0f;
        }
    }
}
void get_coherent_line(
	cv::Mat &src_gray_im, cv::Mat &src_etf, FlowPath** src_fpath, cv::Mat &dst_imCL,
 	int w, int h, float threshold_T, float CL_tanh_he_thr,
	float sigmaC, float sigmaM, float P, int iterations,
	int num_workers
)
{
	printf("SetCoherentLine.. ");
    float** oldim = new float*[w];
    for (int i = 0; i < w; i++) {
        oldim[i] = new float[h];
        for (int j = 0; j < h; j++) {
            oldim[i][j] = src_gray_im.at<float>(j, i);
        }
    }

    // Dog mask calculation
    int MaskSize = threshold_T;
    float *Mc = new float[MaskSize];
    float *Ms = new float[MaskSize];
    float *DogMask = new float[MaskSize];
    float sigmaS = sigmaC * 1.6f;
    BuildGaussianMask(MaskSize, sigmaC, Mc);
    BuildGaussianMask(MaskSize, sigmaS, Ms);
    ScalarBuf(Ms, P, MaskSize);
    DifferenceBuf(DogMask, Mc, Ms, MaskSize);

    std::vector<std::thread> threads(num_workers);
    int slice = w / num_workers;

    for (int iter = 0; iter < iterations; ++iter) {
        for (int t = 0; t < num_workers; t++) {
            int start = t * slice;
            int end = (t == num_workers - 1) ? w : (t + 1) * slice;
            threads[t] = std::thread(__get_coherent_line_partial, start, end, h, oldim, src_fpath, DogMask, sigmaM, threshold_T, CL_tanh_he_thr, std::ref(dst_imCL));
        }

        for (auto &thread : threads) {
            thread.join();
        }
    }

    delete [] Mc;
    delete [] Ms;
    delete [] DogMask;
    for (int i = 0; i < w; i++) {
        delete [] oldim[i];
    }
    delete [] oldim;

    printf("done\n");
}
cv::Mat get_coherent_line(
	cv::Mat &src_gray_im, cv::Mat &src_etf, FlowPath** src_fpath,
 	float threshold_T, float CL_tanh_he_thr,
	float sigmaC, float sigmaM, float P, int iterations
)
{
	if (src_gray_im.rows != src_etf.rows || src_gray_im.cols != src_etf.cols) {
		printf("src_gray and src_etf must have shape shape (width & height)\n");
		exit(1);
	}
	int w = src_gray_im.cols, h = src_gray_im.rows;
	cv::Mat imCL = cv::Mat::zeros(h, w, CV_32F);
	get_coherent_line(src_gray_im, src_etf, src_fpath, imCL, w, h, threshold_T, CL_tanh_he_thr, sigmaC, sigmaM, P, iterations);

	return imCL;
}
cv::Mat get_coherent_line(
	cv::Mat &src_gray_im, cv::Mat &src_etf, FlowPath** src_fpath,
 	float threshold_T, float CL_tanh_he_thr,
	float sigmaC, float sigmaM, float P, int iterations,
	int num_workers
)
{
	if (src_gray_im.rows != src_etf.rows || src_gray_im.cols != src_etf.cols) {
		printf("src_gray and src_etf must have shape shape (width & height)\n");
		exit(1);
	}
	int w = src_gray_im.cols, h = src_gray_im.rows;
	cv::Mat imCL = cv::Mat::zeros(h, w, CV_32F);
	get_coherent_line(src_gray_im, src_etf, src_fpath, imCL, w, h, threshold_T, CL_tanh_he_thr, sigmaC, sigmaM, P, iterations, num_workers);

	return imCL;
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