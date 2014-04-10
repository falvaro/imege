/*
* Copyright (C) 2012 Francisco Álvaro <falvaro@dsic.upv.es>.
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 3 of the License, or (at
* your option) any later version.
*
* This program is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
*/

#include <limits.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pgmtextfea.h"


/* Write image to file */
void write_pgm(int **img, int W, int H, char *pfile) {
  int i, j;
  FILE *fd;

  fd=fopen(pfile, "w");
  if( !fd ) {
    fprintf(stderr, "Error creando archivo '%s'\n", pfile);
    exit(-1);
  }

  fprintf(fd, "P2\n%d %d\n255\n", W, H);
  for(i=0; i<H; i++) {
    for(j=0; j<W; j++)
      fprintf(fd, " %3d", img[i][j]);
    fprintf(fd, "\n");
  }

  fclose(fd);
}


/*
 * Otsu's method code based on the implementation in file
 * 'ocr-binarize-otsu.cc' from the project OCRopus
 * 
 * http://code.google.com/p/ocropus/
 * The OCRopus(tm) open source document analysis and OCR system
 *
 */

void otsu(float **img, int X, int Y, int maxv) {
  int x, y, i, j, N=X*Y;
  int threshold;
  float pdf[256];   /* probability distribution */
  float cdf[256];   /* cumulative probability distribution */
  float myu[256];   /* mean value for separation */
  float sigma[256]; /* inter-class variance */
  float p1p2, mu1mu2diff, max_sigma;

  /* Compute the histogram */
  for(i=0; i<256; i++)
    pdf[i] = 0.0;

  for(y=0; y<Y; y++)
    for(x=0; x<X; x++)
      pdf[ (int)(255*(img[y][x]/maxv)) ] += 1.0;

  /* Histogram normalization (probability distribution) */
  for(i=0; i<256; i++)
    pdf[i] /= N;

  /* cdf & myu generation */
  cdf[0] = pdf[0];
  myu[0] = 0.0;       /* 0.0 times prob[0] equals zero */
  for(i=1; i<256; i++){
    cdf[i] = cdf[i-1] + pdf[i];
    myu[i] = myu[i-1] + i*pdf[i];
  }

  /* sigma maximization
     sigma stands for inter-class variance
     and determines optimal threshold value */
  threshold = 0;
  max_sigma = 0.0;
  for(i=0; i<255; i++){
    if(cdf[i] != 0.0 && cdf[i] != 1.0){
      p1p2 = cdf[i]*(1.0 - cdf[i]);
      mu1mu2diff = myu[256-1]*cdf[i]-myu[i];
      sigma[i] = mu1mu2diff * mu1mu2diff / p1p2;
    }
    else
      sigma[i] = 0.0;
    if(sigma[i] > max_sigma){
      max_sigma = sigma[i];
      threshold = i;
    }
  }

  /* Binarize the image using the computed threshold */
  for(x=0; x<X; x++)
    for(y=0; y<Y; y++)
      img[y][x] = ( 255*(img[y][x]/maxv) > threshold ) ? 1.0 : 0.0;

  /* printf("Otsu method, threshold=%d\n", threshold); */
}

void printInfo(char *str) {
  fprintf(stderr, "Compute BIDM(w,c) between 2 images\n\n");
  fprintf(stderr, "Usage: %s img.pgm ref.pgm w c [options]\n\n", str);
  fprintf(stderr, "    img   - Test image\n");
  fprintf(stderr, "    ref   - Reference image\n");
  fprintf(stderr, "    w     - Warp range (w x w)\n");
  fprintf(stderr, "    c     - Context-window size (c x c)\n\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  -m file - Save mapping values in 'file' (gray)\n");
  fprintf(stderr, "  -b file - Save mapping values in 'file' (binary)\n");
  fprintf(stderr, "  -e file - Save test image pixels evaluation in 'file' (gray)\n");
  fprintf(stderr, "  -f      - Save image derivatives in '*.dh.pgm' and '*.dv.pgm'\n");
  fprintf(stderr, "  -d size - Set derivative computation window size (default 3)\n");
  fprintf(stderr, "  -i lev  - Set inverpolation level {0,1,2} (default 0, no interpolation)\n");
}

float bidme(int **imgh, int **imgv, int I, int J,
	    int **refh, int **refv, int X, int Y,
	    float **warp, int w, int c, FILE *fm, int itp) {

  float min, d, MAX, XI, YJ;
  float v1, v2, aux;
  int i, j, m, n, ip, jp;
  int ha, hb, va, vb, x, y;

  XI = (float)X/I;
  YJ = (float)Y/J;
  MAX = 0;
  
  for(i=0; i<I; i+=1+itp)
    for(j=0; j<J; j+=1+itp) {
      ip = (int)(i*XI+0.5);
      jp = (int)(j*YJ+0.5);

      min=FLT_MAX;

      /* First direct mapping, then the remaining posibilities */
      d = 0;

      /* For each pixel in context window (C x C) */
      for(m=-c/2; m<=c/2 && d<min; m++)
	for(n=-c/2; n<=c/2 && d<min; n++) {
	  /* Compute the difference */

	  ha = hb = va = vb = 127; /* "Background" value */
	  
	  if( i+n >= 0 && i+n < I && j+m >= 0 && j+m < J ) {
	    ha = imgh[j+m][i+n];
	    va = imgv[j+m][i+n];
	  }
	  if( ip+n >= 0 && ip+n < X && jp+m >= 0 && jp+m < Y ) {
	    hb = refh[jp+m][ip+n];
	    vb = refv[jp+m][ip+n];
	  }
	  
	  d += ((ha-hb)*(ha-hb) + (va-vb)*(va-vb))*(1.0/(255*255));
	}
      
      /* Store the minimum value (pruning) */
      if( d < min )
	min = d;
      
      /* Remaining posibilities in warping window W x W */
      for(x=ip-w; x<=ip+w && min>0.0; x++)
	if( x >= 0 && x < X ) {
	  for(y=jp-w; y<=jp+w && min>0.0; y++) 
	    if( y >= 0 && y < Y ) {
	      d = 0;

	      /* For each pixel in context window (C x C) */
	      for(m=-c/2; m<=c/2 && d<min; m++)
		for(n=-c/2; n<=c/2 && d<min; n++) {
		  /* Compute the difference */

		  ha = hb = va = vb = 127; /* "Background" value */
		  
		  if( i+n >= 0 && i+n < I && j+m >= 0 && j+m < J ) {
		    ha = imgh[j+m][i+n];
		    va = imgv[j+m][i+n];
		  }
		  if( x+n >= 0 && x+n < X && y+m >= 0 && y+m < Y ) {
		    hb = refh[y+m][x+n];
		    vb = refv[y+m][x+n];
		  }

		  d += ((ha-hb)*(ha-hb) + (va-vb)*(va-vb))*(1.0/(255*255));
		}

	      /* Store the minimum value (pruning) */
	      if( d < min )
		min = d;
	    }
	}

      /* Store the maximum value */
      if( d > MAX )
	MAX=d;

      /* Set the warping value of pixel (i,j) */
      warp[j][i] = d;      
    }

  if( MAX > 0 ) {

    if( itp == 1 ) {
      //¡* Interpolation of not computed values */
      for(j=0; j<J; j+=2)
	for(i=1; i<I; i+=2) {
	  v1=(i-1>=0) ? warp[j][i-1] : 0;
	  v2=(i+1<I)  ? warp[j][i+1] : 0;
	  warp[j][i] = (v1+v2)/2;
	}
      for(j=1; j<J; j+=2)
	for(i=0; i<I; i++) {
	  v1=(j-1>=0) ? warp[j-1][i] : 0;
	  v2=(j+1<J)  ? warp[j+1][i] : 0;
	  warp[j][i] = (v1+v2)/2;
	}
    }
    else if( itp == 2 ) {
      /* nterpolation of not computed values */
      for(j=0; j<J; j+=3)
	for(i=1; i<I; i+=3) {
	  v1=(i-1>=0) ? warp[j][i-1] : 0;
	  v2=(i+2<I)  ? warp[j][i+2] : 0;
	  aux=(v2-v1)/3;
	  warp[j][i]   = v1+aux;
	  if( i+1 < I )
	    warp[j][i+1] = v2-aux;
	}
      for(j=1; j<J; j+=3)
	for(i=0; i<I; i++) {
	  v1=(j-1>=0) ? warp[j-1][i] : 0;
	  v2=(j+2<J)  ? warp[j+2][i] : 0;
	  aux=(v2-v1)/3;
	  warp[j][i]   = v1+aux;
	  if( j+1 < J )
	    warp[j+1][i] = v2-aux;
	}
    }

    /* Set mapping values to range [0,255] */
    for(y=0; y<J; y++)
      for(x=0; x<I; x++)
	warp[y][x] = 255*(warp[y][x]/MAX);

    if( fm ) {
      /* Save the gray-scale mapping image */
      fprintf(fm, "P2\n%d %d\n255\n", I, J);
      for(j=0; j<J; j++) {
	for(i=0; i<I; i++)
	  fprintf(fm, " %d", 255-(int)warp[j][i]);
	fprintf(fm, "\n");
      }
      fclose(fm);
    }

    /* Mapping values binarization */
    otsu(warp, I, J, 255);
  }

  return MAX;
}



int main(int argc, char *argv[]) {
  int **imgh, **imgv, **refh, **refv;
  int I, J, X, Y, w, c;
  int i, j, cp, fg;
  int option, om=0, ob=0, oe=0, of=0, od=3, proc=5, ip=0;
  FILE *fm=NULL, *fb=NULL, *fe=NULL;
  char aux[512];
  float **warp, MAX;
  gray **img, **ref;

  if( argc < 5 ) {
    printInfo(argv[0]);
    return -1;
  }

  w = atoi(argv[3]);
  if( w <= 0 ) {
    fprintf(stderr, "Error: parameter 'w' must be a positive integer\n");
    return -1;
  }

  c = atoi(argv[4]);
  if( c <= 0 || c%2==0 ) {
    fprintf(stderr, "Error: parameter 'c' must be a odd positive integer\n");
    return -1;
  }

  while( (option=getopt(argc-4,&argv[4],"fm:b:d:e:i:"))!=-1) {
    proc++;
    switch (option) {
    case 'm':
      proc++;
      om=1;
      fm = fopen(optarg, "w");
      if( !fm ) {
	fprintf(stderr, "Error creating file '%s' (option -m)\n\n", optarg);
	printInfo(argv[0]);
	return -1;
      }
      break;
    case 'b':
      proc++;
      ob=1;
      fb = fopen(optarg, "w");
      if( !fb ) {
	fprintf(stderr, "Error creating file '%s' (option -b)\n\n", optarg);
	printInfo(argv[0]);
	return -1;
      }
      break;
    case 'e':
      proc++;
      oe=1;
      fe = fopen(optarg, "w");
      if( !fe ) {
	fprintf(stderr, "Error creating file '%s' (option -e)\n\n", optarg);
	printInfo(argv[0]);
	return -1;
      }
      break;
    case 'd':
      proc++;
      od=atoi(optarg);
      if( od<=0 ) {
	fprintf(stderr, "Error. Size must be a positive integer (option -d)\n\n", optarg);
	printInfo(argv[0]);
	return -1;
      }
      break;
    case 'i':
      proc++;
      ip=atoi(optarg);
      if( ip!=0 && ip!=1 && ip!=2 ) {
	fprintf(stderr, "Error. Interpolation must be 0, 1 or 2 (option -i)\n\n", optarg);
	printInfo(argv[0]);
	return -1;
      }
      break;
    case 'f':
      of=1;
      break;
    default:
      printInfo(argv[0]);
      return -1;
    }
  }

  if( proc != argc ) {
    fprintf(stderr, "Error: Unknown arguments\n\n");
    printInfo(argv[0]);
    return -1;
  }

  img=pgmtextfea(argv[1], &imgh, &imgv, &I, &J, od);
  ref=pgmtextfea(argv[2], &refh, &refv, &X, &Y, od);

  if( of ) {
    /* Save the derivatives image representation */
    strcpy(aux, argv[1]);  strcat(aux, ".dh.pgm");
    write_pgm(imgh, I, J, aux);
    strcpy(aux, argv[1]);  strcat(aux, ".dv.pgm");
    write_pgm(imgv, I, J, aux);
    strcpy(aux, argv[2]);  strcat(aux, ".dh.pgm");
    write_pgm(refh, X, Y, aux);
    strcpy(aux, argv[2]);  strcat(aux, ".dv.pgm");
    write_pgm(refv, X, Y, aux);
  }

  warp  = (float **)malloc(sizeof(float*)*J);
  for(i=0; i<J; i++) 
    warp[i]  = (float*)malloc(sizeof(float)*I);


  /* Compute bidme( test, ref, w, c ) */
  MAX=bidme(imgh, imgv, I, J, refh, refv, X, Y, warp, w, c, fm, ip);


  if( ob ) {
    /* Save the binary mapping image */
    fprintf(fb, "P1\n%d %d\n", I, J);
    for(j=0; j<J; j++) {
      for(i=0; i<I; i++)
	fprintf(fb, "%d", (int)warp[j][i]);
      fprintf(fb, "\n");
    }
    fclose(fb);
  }

  /* Correct pixels */
  cp=0;
  if( MAX > 0 ) {
    /* Foreground pixels */
    fg=0;

    for(i=0; i<I; i++)
      for(j=0; j<J; j++) {
	
	if( img[j][i]<255 ) {
	  fg++;
	  if( warp[j][i] == 0.0 )
	    cp++;
	}
	else
	  if( warp[j][i] > 0.0 )
	    warp[j][i] = 0;

      }
  }
  else {
    cp=fg=1;
    /* Set the quotient cp/fg to 1 */
  }
  
  printf("%8.6f\n", (float)cp/fg);

  if( oe ) {
    /* Save the intersection between the binary mapping and
       the test expression foreground pixels */

    fprintf(fe, "P2\n%d %d\n4\n", I, J);
    for(j=0; j<J; j++) {
      for(i=0; i<I; i++) {
	if( img[j][i]==255 )
	  fprintf(fe, " 4");
	else {
	  if( warp[j][i]==0.0 ) fprintf(fe, " 3");
	  else                  fprintf(fe, " 0");
	}
      }
      fprintf(fe, "\n");
    }
    fclose(fe);
  }


  return 0;
}
