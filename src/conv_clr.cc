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

	//Observer. = 2��, Illuminant = D65
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