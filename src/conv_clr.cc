#include "conv_clr.h"


void conv_RGBtoLAB(V3DF** buf, int width, int height)		//RGB -> LAB
{
	int i,j;

	for (i=0; i<width; i++){
		for (j=0; j<height; j++){
			RGBtoLAB(buf[i][j][0], buf[i][j][1], buf[i][j][2], &buf[i][j][0], &buf[i][j][1], &buf[i][j][2]);
		}
	}
}

void conv_LABtoRGB(V3DF** buf, int width, int height)		//LAB -> RGB
{
	int i,j;

	for (i=0; i<width; i++){
		for (j=0; j<height; j++){
			LABtoRGB(buf[i][j][0], buf[i][j][1], buf[i][j][2], &buf[i][j][0], &buf[i][j][1], &buf[i][j][2]);
		}
	}
}

void RGBtoLAB(float iR, float iG, float iB, float *oL, float *oa, float *ob){
	double r, g, b;
	double x, y, z;
	double l, a, b1;
	const double REF_X = 95.047; // Observer= 2��, Illuminant= D65
	const double REF_Y = 100.000; 
	const double REF_Z = 108.883; 

	b = iR;
	g = iG;
	r = iB;

	if (r > 0.04045)
	{
		r = (r+0.055) / 1.055;
		r = pow(r, 2.4); 
	}
	else { r = r / 12.92; }

	if ( g > 0.04045)
	{ 
		g = (g+0.055) / 1.055;
		g = pow(g, 2.4); 
	}
	else { g = g / 12.92; }

	if (b > 0.04045)
	{
		b = (b+0.055) / 1.055;
		b = pow(b, 2.4); 
	}
	else {	b = b / 12.92; }

	r = r * 100;
	g = g * 100;
	b = b * 100;

	//Observer. = 2°, Illuminant = D65
	x = r * 0.4124 + g * 0.3576 + b * 0.1805;
	y = r * 0.2126 + g * 0.7152 + b * 0.0722;
	z = r * 0.0193 + g * 0.1192 + b * 0.9505;

	x /= REF_X;   
	y /= REF_Y;  
	z /= REF_Z;

	if ( x > 0.008856 ) { x = pow( x , 0.33333 ); } // 1/3
	else { x = ( 7.787 * x ) + ( 0.137931 ); } //16/116

	if ( y > 0.008856 ) { y = pow( y , 0.33333 ); }
	else { y = ( 7.787 * y ) + ( 0.137931 ); }

	if ( z > 0.008856 ) { z = pow( z , 0.33333 ); }
	else { z = ( 7.787 * z ) + ( 0.137931 ); }

	l = ( 116.0 * y ) - 16.0;
	a = 504.3 * ( x - y );
	b1 = 201.7 * ( y - z );

	if(l > 100)
		(*oL) = 100;
	else if(l < 0)
		(*oL) = 0;
	else
		(*oL) = (char) ( (116.0 * y) - 16.0 ); // l
	
	if(a > 127)
		(*oa) = 127;
	else if(a < -127)
		(*oa) = -127;
	else
		(*oa) = (char) ( 504.3 * (x-y) ); // a
	
	if(b > 127)
		(*ob) = 127;
	else if( b < -127)
		(*ob) = -127;
	else
		(*ob) = (char) ( 201.7 * (y-z) );  // b
}

void LABtoRGB(float iL, float ia, float ib, float *oR, float *oG, float *oB){
	const double REF_X = 95.047; // Observer= 2��, Illuminant= D65
	const double REF_Y = 100.000; 
	const double REF_Z = 108.883; 

	char l, a, b1;
	double x, y, z;
	double r, g, b;

	l = iL;
	a = ia;
	b1 = ib;

	y = (l + 16) / 116.0;
	x = a / 504.3 + y; // 500
	z = y - (b1 / 201.7); // 200

	if ( pow( y , 3 ) > 0.008856 ) { y = pow( y , 3 ); }
	else { y = ( y - 0.137931 ) / 7.787; }
	if ( pow( x , 3 ) > 0.008856 ) { x = pow( x , 3 ); }
	else { x = ( x - 0.137931 ) / 7.787; }
	if ( pow( z , 3 ) > 0.008856 ) { z = pow( z , 3 ); }
	else { z = ( z - 0.137931 ) / 7.787; }

	x = REF_X * x;     
	y = REF_Y * y; 
	z = REF_Z * z; 

	x /= 100.0;
	y /= 100.0;
	z /= 100.0;

	r = x * 3.2406 + y * -1.5372 + z * -0.4986;
	g = x * -0.9689 + y * 1.8758 + z * 0.0415;
	b = x * 0.0557 + y * -0.2040 + z * 1.0570;

	if ( r > 0.0031308 ) { r = 1.055 * pow( r , 0.416667 ) - 0.055; } //1/2.4
	else { r = 12.92 * r; }
	if ( g > 0.0031308 ) { g = 1.055 * pow( g , 0.416667 ) - 0.055; }
	else { g = 12.92 * g; }
	if ( b > 0.0031308 ) { b = 1.055 * pow( b , 0.416667 ) - 0.055; }
	else { b = 12.92 * b; }

	if(r > 1.0) 
	{
		(*oB) = 1.0f;
	}
	else if(r < 0.0)
	{
		(*oB) = 0.0f;
	}
	else
		(*oB) = r;
	
	if(g > 1.0) 
	{
		(*oG) = 1.0f;
	}
	else if(g < 0.0)
	{
		(*oG) = 0.0f;
	}
	else
		(*oG) = g;
	
	if(b > 1.0) 
	{
		(*oR) = 1.0f;
	}
	else if(b < 0.0)
	{
		(*oR) = 0.0f;
	}
	else
		(*oR) = b;
}

void conv_RGBtoLUV(V3DF** buf, int width, int height)		//RGB -> LUV
{
	conv_RGBtoXYZ(buf, width, height);	//RGB->XYZ
	conv_XYZtoLUV(buf, width, height);	//XYZ->LUV
}

void conv_LUVtoRGB(V3DF** buf, int width, int height)		//LUV -> RGB
{
	conv_LUVtoXYZ(buf, width, height);	//LUV->XYZ
	conv_XYZtoRGB(buf, width, height);	//XYZ->RGB
}

void conv_RGBtoXYZ(V3DF** buf, int width, int height)		//RGB -> XYZ
{
	int i,j;

	for (i=0; i<width; i++){
		for (j=0; j<height; j++){
			RGBtoXYZ(buf[i][j][R_], buf[i][j][G_], buf[i][j][B_], &buf[i][j][X], &buf[i][j][Y], &buf[i][j][Z]);			
		}
	}
}

void conv_XYZtoRGB(V3DF** buf, int width, int height)		//XYZ -> RGB
{
	int i,j;

	for (i=0; i<width; i++){
		for (j=0; j<height; j++){
			XYZtoRGB(buf[i][j][X], buf[i][j][Y], buf[i][j][Z], &buf[i][j][R_], &buf[i][j][G_], &buf[i][j][B_]);
		}
	}
}

void conv_XYZtoLUV(V3DF** buf, int width, int height)		//XYZ -> LUV
{
	int i,j;

	for (i=0; i<width; i++){
		for (j=0; j<height; j++){
			XYZtoLUV(buf[i][j][X], buf[i][j][Y], buf[i][j][Z], &buf[i][j][X], &buf[i][j][Y], &buf[i][j][Z]);
		}
	}
}

void conv_LUVtoXYZ(V3DF** buf, int width, int height)		//LUV -> XYZ
{
	int i,j;

	for (i=0; i<width; i++){
		for (j=0; j<height; j++){
			LUVtoXYZ(buf[i][j][X], buf[i][j][Y], buf[i][j][Z], &buf[i][j][X], &buf[i][j][Y], &buf[i][j][Z]);
		}
	}
}

void RGBtoXYZ(float iR, float iG, float iB, float *oX, float *oY, float *oZ)
{
	V3DF res;
	V3DF iRGB;		//vector3(iRGB, iR, iG, iB);
	vector(iRGB, iR, iG, iB);

	Matrix3DF m_RGBtoXYZ = {0.607f, 0.174f, 0.200f, 0.299f, 0.587f, 0.114f, 0.000f, 0.066f, 1.116f};

	// m_multi_v(res, m_RGBtoXYZ, iRGB);
	mult_mat_pt_3(res, m_RGBtoXYZ, iRGB);

	(*oX) = res[X];
	(*oY) = res[Y];
	(*oZ) = res[Z];
}

void XYZtoRGB(float iX, float iY, float iZ, float *oR, float *oG, float *oB)
{
	V3DF res;
	V3DF iXYZ;		vector(iXYZ, iX, iY, iZ);

	Matrix3DF m_XYZtoRGB = {1.91046f, -0.53394f, -0.287834f, -0.984436f, 1.9985f, -0.027726f, 0.0582193f, -0.118191f, 0.897697f};

	// m_multi_v(res, m_XYZtoRGB, iXYZ);
	mult_mat_pt_3(res, m_XYZtoRGB, iXYZ);

	(*oR) = res[X];
	(*oG) = res[Y];
	(*oB) = res[Z];
}

void XYZtoLUV(float iX, float iY, float iZ, float *oL, float *oU, float *oV)
{
	float x0,y0,z0;
	float u0,v0;

	x0 = 0.98072;
	y0 = 1.00000;
	z0 = 1.18225;
	/*x0 = 1.00000;
	y0 = 1.00000;
	z0 = 1.00000;*/

	u0 = 4.0f*x0 / (x0 + 15.0f*y0 + 3.0f*z0);
	v0 = 9.0f*y0 / (x0 + 15.0f*y0 + 3.0f*z0);

	float u_,v_;
	if ( (iX + 15.0f*iY + 3.0f*iZ)==0.0f )
		u_ = 0.0f;
	else
		u_ = 4.0f*iX / (iX + 15.0f*iY + 3.0f*iZ);

	if ( (iX + 15.0f*iY + 3.0f*iZ)==0.0f )
		v_ = 0.0f;
	else
		v_ = 9.0f*iY / (iX + 15.0f*iY + 3.0f*iZ);

	if (iY/y0 >= 0.008856)
		(*oL) = 25.0f * pow ( 100.0f * (iY/y0), 1.0f/3.0f ) - 16.0f;
	else if (iY/y0 < 0.008856)
		(*oL) = 903.3f * (iY/y0);

	(*oU) = 13.0f*(*oL)*(u_ - u0);
	(*oV) = 13.0f*(*oL)*(v_ - v0);
//	(*oU) = (*oL);
//	(*oV) = (*oL);
}

void LUVtoXYZ(float iL, float iU, float iV, float *oX, float *oY, float *oZ)
{
	float x0,y0,z0;
	float u0,v0;

	x0 = 0.98072;
	y0 = 1.00000;
	z0 = 1.18225;

	u0 = 4.0f*x0 / (x0 + 15.0f*y0 + 3.0f*z0);
	v0 = 9.0f*y0 / (x0 + 15.0f*y0 + 3.0f*z0);

	float u_ = iU/(13.0f*iL) + u0;
	float v_ = iV/(13.0f*iL) + v0;

	(*oY) = (y0/100.0f) * pow( (iL+16.0f) / 25.0f, 3.0f );
	(*oX) = ( (9.0f*u_)/(4.0f*v_) ) * (*oY);
	(*oZ) = (*oY) * ( (12.0f - 3.0f*u_ - 20.0f*v_) / (4.0f*v_) );
}

void LUVtoRGB(float iL, float iU, float iV, float *oR, float *oG, float *oB)
{
	LUVtoXYZ(iL,iU,iV, oR,oG,oB);
	XYZtoRGB((*oR),(*oG),(*oB), oR,oG,oB);
}

void RGBtoLUV(float iR, float iG, float iB, float *oL, float *oU, float *oV)
{
	RGBtoXYZ(iR,iG,iB, oL,oU,oV);
	XYZtoLUV((*oL),(*oU),(*oV), oL,oU,oV);
}