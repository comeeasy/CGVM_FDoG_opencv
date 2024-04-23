#include "common_func.h"


typedef float V3DF[3];

void init_buf( V3DF **OutImg, V3DF **InImg, int width, int height) {
	int i,j;
	//OutImg �ʱ�ȭ
	for (i = 0; i < width; i++){
		for (j = 0; j < height; j++){
			vcopy(OutImg[i][j], InImg[i][j]);
		}
	}
}

void GetValAtPoint_V3DF(V3DF **ival, int w, int h, float px, float py, V3DF oval)
{
	int ulx,uly,lrx,lry;		//	upper-left (x, y), lower-right (x, y)
	V3DF ivalA,ivalB,ivalC,ivalD;	//	etf at A,B,C,D
	V3DF e1,e2,t1,t2;

	if (px+1>=w || py+1>=h){
		vector(oval, 0.0f,0.0f,0.0f);
		return;
	}

	ulx = (int)(px);
	uly = (int)(py+1);
	lrx = (int)(px+1);
	lry = (int)(py);
	vcopy(ivalA, ival[ulx][lry]);
	vcopy(ivalB, ival[lrx][lry]);
	vcopy(ivalC, ival[ulx][uly]);
	vcopy(ivalD, ival[lrx][uly]);

	//distA: A~(px,py) �Ÿ�
	//eA: etf at A
	//1. distB*etfA + distA*etfB
	vcopy(t1,ivalA);	vcopy(t2,ivalB);
	vscale(t1,lrx-px);	vscale(t2,px-ulx);
	vadd(e1,t1,t2);
	//2. distD*etfC + distC*etfD
	vcopy(t1,ivalC);	vcopy(t2,ivalD);
	vscale(t1,lrx-px);	vscale(t2,px-ulx);
	vadd(e2,t1,t2);
	//3. distC*1 + distA*2
	vcopy(t1,e1);		vcopy(t2,e2);
	vscale(t1,uly-py);	vscale(t2,py-lry);
	vadd(oval,t1,t2);

	vnorm(oval);
	//return vlength(resETF);
}

void GetValAtPoint(float **ival, int w, int h, float px, float py, float *oval)
{
	int ulx,uly,lrx,lry;		//	upper-left (x, y), lower-right (x, y)
	float ivalA,ivalB,ivalC,ivalD;	//	etf at A,B,C,D
	float e1,e2,t1,t2;

	if (px+1>=w || py+1>=h){
		*oval = 0.0f;
		return;
	}

	ulx = (int)(px);
	uly = (int)(py+1);
	lrx = (int)(px+1);
	lry = (int)(py);
	ivalA = ival[ulx][lry];
	ivalB = ival[lrx][lry];
	ivalC = ival[ulx][uly];
	ivalD = ival[lrx][uly];

	//distA: A~(px,py) �Ÿ�
	//eA: etf at A
	//1. distB*etfA + distA*etfB
	t1 = ivalA;		t2 = ivalB;
	t1 *= (lrx-px);	t2 *= (px-ulx);
	e1 = t1+t2;
	//2. distD*etfC + distC*etfD
	t1 = ivalC;		t2 = ivalD;
	t1 *= (lrx-px);	t2 *= (px-ulx);
	e2 = t1+t2;
	//3. distC*1 + distA*2
	t1 = e1;		t2 = e2;
	t1 *= (uly-py);	t2 *= (py-lry);
	*oval = t1+t2;
}

void GaussianFilter(float **InImg, float **OutImg, int width, int height)
{
	int i,j;
	int row, col;
	int rowOffset;			// Row offset from the current pixel
	int colOffset;			// Col offset from the current pixel
	int rowTotal = 0;		// Row position of offset pixel
	int colTotal = 0;		// Col position of offset pixel
	float newPixel;
	float gaussianMask[5][5];		// Gaussian mask
	float **img_temp = new float*[height];
	for (i=0; i<height; i++)
		img_temp[i] = new float[width];

	/* Declare Gaussian mask */
	gaussianMask[0][0] = 2.0f;		gaussianMask[0][1] = 4.0f;		gaussianMask[0][2] = 5.0f;		gaussianMask[0][3] = 4.0f;		gaussianMask[0][4] = 2.0f;	
	gaussianMask[1][0] = 4.0f;		gaussianMask[1][1] = 9.0f;		gaussianMask[1][2] = 12.0f;		gaussianMask[1][3] = 9.0f;		gaussianMask[1][4] = 4.0f;	
	gaussianMask[2][0] = 5.0f;		gaussianMask[2][1] = 12.0f;		gaussianMask[2][2] = 15.0f;		gaussianMask[2][3] = 12.0f;		gaussianMask[2][4] = 2.0f;	
	gaussianMask[3][0] = 4.0f;		gaussianMask[3][1] = 9.0f;		gaussianMask[3][2] = 12.0f;		gaussianMask[3][3] = 9.0f;		gaussianMask[3][4] = 4.0f;	
	gaussianMask[4][0] = 2.0f;		gaussianMask[4][1] = 4.0f;		gaussianMask[4][2] = 5.0f;		gaussianMask[4][3] = 4.0f;		gaussianMask[4][4] = 2.0f;

	/* Gaussian Blur */
	for(row = 0; row < height; row++){
		for (col = 0; col < width; col++){
			newPixel = 0.0f;
			for (rowOffset = -2; rowOffset <= 2; rowOffset++){
				for (colOffset = -2; colOffset <= 2; colOffset++){
					rowTotal = row + rowOffset;
					colTotal = col + colOffset;
					//���ó��(��Ī)
					if (rowTotal<0) rowTotal *= (-1.0f);
					else if (rowTotal>=height)	rowTotal = height-(rowTotal-(height-1));
					if (colTotal<0) colTotal *= (-1.0f);
					else if (colTotal>=width)	colTotal = width-(colTotal-(width-1));
					//iOffset = (unsigned long)(rowTotal*3*W + colTotal*3);	//BMP�� 1���� �迭�� �������� ��� ����.
					newPixel += (InImg[colTotal][rowTotal])*gaussianMask[2 + rowOffset][2 + colOffset];
				}
			}
			img_temp[row][col] = (newPixel/159.0f);
		}
	}

	for (i=0; i<width; i++)
		for (j=0; j<height; j++)
			OutImg[i][j] = img_temp[j][i];

	for (i=0; i<height; i++)
		delete []img_temp[i];
	delete []img_temp;
}

void BuildGaussianMask(int n, float s, float *m)
{
	if (n%2==0){
		printf("Gaussian mask ���� ����; ����ũ ũ��� Ȧ������ �մϴ�.\n");
		exit(0);
	}
	int i;
	float div = pow(2.0f*PIF*s*s, 0.5f);
	int ct = n/2;	//center
	float sum=0.0f;
	float x2;
	
	for (i=0; i<n; i++){
		x2 = (i-ct)*(i-ct);
		m[i] = pow(EXPF, -x2/(2.0f*s*s)) / div;
		sum += m[i];
	}
	
	for (i=0; i<n; i++){
		m[i] /= sum;
	}
}

void ScalarBuf(float *OutBuf, float scalar, int length)
{
	for (int i=0; i<length; i++)
		OutBuf[i] = OutBuf[i] * scalar;
}

void DifferenceBuf(float *OutBuf, float *InBuf1, float *InBuf2, int length)
{
	for (int i=0; i<length; i++){
		OutBuf[i] = InBuf1[i] - InBuf2[i];
	}
}