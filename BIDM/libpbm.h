/* pbm.h - header file for libpbm portable bitmap library
*/

#ifndef _PBM_H_
#define _PBM_H_

#include <stdio.h>

#define PBM_WHITE 0
#define PBM_BLACK 1


/* Magic constants. */

#define PBM_MAGIC1 'P'
#define PBM_MAGIC2 '1'
#define RPBM_MAGIC2 '4'
#define PBM_FORMAT (PBM_MAGIC1 * 256 + PBM_MAGIC2)
#define RPBM_FORMAT (PBM_MAGIC1 * 256 + RPBM_MAGIC2)
#define PBM_TYPE PBM_FORMAT


/* Macro for turning a format number into a type number. */

#define PBM_FORMAT_TYPE(f) \
  ((f) == PBM_FORMAT || (f) == RPBM_FORMAT ? PBM_TYPE : -1)


typedef unsigned char bit;

/* Declarations of routines. */

bit ** pbm_allocarray(int cols, int rows);
void pbm_freearray(bit ** img, int rows); 

bit** pbm_readpbm(FILE* file, int* colsP, int* rowsP);

void pbm_writepbm (FILE* file, bit** bits, int cols, int rows, int forceplain);


#endif /*_PBM_H_*/

