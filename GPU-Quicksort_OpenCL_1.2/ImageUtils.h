/* ************************************************************************* *\
               INTEL CORPORATION PROPRIETARY INFORMATION
     This software is supplied under the terms of a license agreement or 
     nondisclosure agreement with Intel Corporation and may not be copied 
     or disclosed except in accordance with the terms of that agreement. 
        Copyright (C) 2014 Intel Corporation. All Rights Reserved.
\* ************************************************************************* */
// junkins image format and ppm file utilities


#include <stdio.h>
#include <stdlib.h>

#define U8	unsigned char
#define U16	wchar_t
#define F32	float

void convertU8RGB_to_U8Gray (int w, int h, U8* pSrcU8_RGB, U8* pDstU8_G)
{
	// We’re more sensitive to green than other colors, so green is weighted most heavily. 
	// The formula for luminosity is 0.21 R + 0.71 G + 0.07 B

	unsigned int strideU8_RGB	= w * 3;	// 3 bytes for each RGB
	unsigned int strideU8_G		= w;		// 1 byte for each grayscale

	for (int y=0; y < h; y++)
	{
		for (int x=0; x < w; x++)
		{
			unsigned int locU8_RGB	= y*strideU8_RGB	+ x*3;
			unsigned int locU8_G	= y*strideU8_G		+ x;

			pDstU8_G[locU8_G]	= (U8) (	0.21f * pSrcU8_RGB[locU8_RGB]
										+	0.71f * pSrcU8_RGB[locU8_RGB+ 1]
										+	0.07f * pSrcU8_RGB[locU8_RGB+ 1]);
		}
	}
}


void convertU8Gray_to_U8RGB (int w, int h, U8* pSrcU8_G,  U8* pDstU8_RGB)
{


	unsigned int strideU8_RGB	= w * 3;	// 3 bytes for each RGB
	unsigned int strideU8_G		= w;		// 1 byte for each grayscale

	for (int y=0; y < h; y++)
	{
		for (int x=0; x < w; x++)
		{
			unsigned int locU8_RGB	= y*strideU8_RGB	+ x*3;
			unsigned int locU8_G	= y*strideU8_G		+ x;

			pDstU8_RGB[locU8_RGB]		= (U8) pSrcU8_G[locU8_G]; // Write gray value as Red
			pDstU8_RGB[locU8_RGB + 1]	= (U8) pSrcU8_G[locU8_G]; // Write gray value as Green
			pDstU8_RGB[locU8_RGB + 2]	= (U8) pSrcU8_G[locU8_G]; // Write gray value as Blue
		}
	}
}


void convertU16Gray_to_U8Gray (int w, int h, U16* pSrcU16, U8* pDstU8)
{
	unsigned int strideU8	= w;
	unsigned int locU8;

	for (int y=0; y < h; y++)
	{
		for (int x=0; x < w; x++)
		{
			locU8	= y*strideU8	+ x;

			pDstU8[locU8]	= (U8) pSrcU16[locU8];
		}
	}
}

void convertU8Gray_to_F32Gray (int w, int h, U8* pSrcU8, F32* pDstF32)
{
	unsigned int strideU8 = w;
	unsigned int locU8;

	for (int y=0; y < h; y++)
	{
		for (int x=0; x < w; x++)
		{
			locU8	= y*strideU8	+ x;

			pDstF32[locU8]	= (F32) pSrcU8[locU8];
		}
	}
}



void convertF32Gray_to_U8Gray (int w, int h, F32* pSrcF32, U8* pDstU8)
{
	unsigned int strideU8 = w;
	unsigned int locU8;

	for (int y=0; y < h; y++)
	{
		for (int x=0; x < w; x++)
		{
			locU8	= y*strideU8	+ x;

			pDstU8[locU8]	= (U8) pSrcF32[locU8];
		}
	}
}


void loadppm_toU8RGB(const char *fname, unsigned int* pW, unsigned int* pH, U8 **img, unsigned int *upImgSize)
{
	char			sPPMHeader[2];
	FILE			*fp;
	char            c;
	unsigned int    uWidth = 0, uHeight = 0, uImgSize = 0, uPPMRange = 0;
	size_t          szRecords;
	U8              *pImg = 0;
	int             iFields = 0;

	// initialize return values to zeros
	*pW = *pH = *upImgSize = 0;
	*img = 0;

	errno_t err = fopen_s(&fp, fname, "rb");
	if( err != 0 )
	{
		printf( "loadppm_toU8RGB: The file %s was not opened: %s\n", fname, strerror(err) );
		exit(-1);
	}

	// read ppm image header info:
	c = getc(fp);
	if (c != EOF) {
		sPPMHeader[0] = c;
	} else {
		printf( "loadppm_toU8RGB: Couldn't read the ppm image header in file %s: %s\n", fname, strerror(err) );
		exit(-1);
	}
	c = getc(fp);
	if (c != EOF) {
		sPPMHeader[1] = c;
	} else {
		printf( "loadppm_toU8RGB: Couldn't read the ppm image header in file %s: %s\n", fname, strerror(err) );
		exit(-1);
	}

	//check for comments
	c = getc(fp);
	if (c == EOF) {
		printf( "loadppm_toU8RGB: Improper ppm image format %s: %s\n", fname, strerror(err) );
		exit(-1);
	}

	while (c == '#' || c == '\n') 
	{
		do 
		{
			c = getc(fp);
			if (c == EOF) {
				printf( "loadppm_toU8RGB: Improper ppm image format %s: %s\n", fname, strerror(err) );
				exit(-1);
			}
		} 
		while (c != '\n');
		c = getc(fp);
		if (c == EOF) {
			printf( "loadppm_toU8RGB: Improper ppm image format %s: %s\n", fname, strerror(err) );
			exit(-1);
		}
	}

	if (ungetc(c, fp) == EOF) {
		printf( "loadppm_toU8RGB: Cannot put the character %c back %s: %s\n", c, fname, strerror(err) );
		exit(-1);
	}

	iFields = fscanf_s(fp, "%d %d", &uWidth, &uHeight);
	if( iFields != 2 )
	{
		printf( "loadppm_toU8RGB: Couldn't read width and height in file %s: %s\n", fname, strerror(err) );
		exit(-1);
	}

	iFields = fscanf_s(fp, "%d\n", &uPPMRange);
	if( iFields != 1 )
	{
		printf( "loadppm_toU8RGB: Couldn't read ppm range in file %s: %s\n", fname, strerror(err) );
		exit(-1);
	}

	uImgSize = 3*uWidth*uHeight;;
	pImg = (U8 *) malloc (uImgSize);
	if (pImg == NULL)
	{
		printf("loadppm_toU8RGB: Couldn't allocate memory of size: %d\n", pImg);
		exit(-1);
	}

	//Read ppm image data:
	szRecords = fread_s(pImg, uImgSize, uImgSize, 1, fp);
	if (szRecords != 1) 
	{
		printf( "Couldn't read image content from file %s: %s\n", fname, strerror(err) );
		exit(-1);
	}

	err = fclose(fp);
	if( err != 0 )
	{
		printf( "Couldn't close file %s: %s\n", fname, strerror(err) );
		exit(-1);
	}

	// return image and its parameters
	*img = pImg;
	*pW  = uWidth;
	*pH  = uHeight;
	*upImgSize = uImgSize;
}


void saveppm_fromU8RGB(const char *fname, unsigned int w, unsigned int h, U8 *img)
{
	FILE *fp;

	fopen_s(&fp, fname, "wb");
	// assert(fp);

	fprintf(fp, "P6\n");
	fprintf(fp, "%d %d\n", w, h);
	fprintf(fp, "255\n");
	fwrite(img, w * h * 3, 1, fp);
	fclose(fp);
}

void resizeU8Gray (unsigned int w, unsigned int h, unsigned int wReSz, unsigned int hReSz, U8* pImgData_orig, U8* pImgData_resz)
{
	unsigned int strideU8_orig = w;
	unsigned int strideU8_resz = wReSz;
	unsigned int locU8_orig, locU8_resz;

	for (unsigned int y=1; y < hReSz-1; y++)
	{
		for (unsigned int x=16; x < wReSz-16; x++)
		{
			locU8_orig	= ((y-1) % h) * strideU8_orig + ((x-16) % w);
			locU8_resz	= y*strideU8_resz + x ;

			pImgData_resz[locU8_resz] = pImgData_orig[locU8_orig]; 
		}
	}

	// Bottom row y = 0
	for (unsigned int x=16; x < wReSz-16; x++)
	{
		locU8_orig	= ((x-16) % w);
		locU8_resz	= x ;

		pImgData_resz[locU8_resz] = pImgData_orig[locU8_orig]; 
	}
	// Last row y = hReSz-1
	for (unsigned int x=16; x < wReSz-16; x++)
	{
		locU8_orig	= ((hReSz-3) % h) * strideU8_orig + ((x-16) % w);
		locU8_resz	= (hReSz-1)*strideU8_resz + x ;

		pImgData_resz[locU8_resz] = pImgData_orig[locU8_orig]; 
	}
	
	// Now for columns
	// x = 15
	for (unsigned int y=0; y < hReSz; y++)
	{
		locU8_resz	= y*strideU8_resz+15;

		pImgData_resz[locU8_resz] = pImgData_resz[locU8_resz+1]; 
	}
	
	// x = wReSz-16
	for (unsigned int y=0; y < hReSz; y++)
	{
		locU8_orig	= y*strideU8_resz + wReSz-16 ;
		locU8_resz	= y*wReSz + wReSz-16 ;

		pImgData_resz[locU8_resz] = pImgData_resz[locU8_orig-1]; 
	}
}