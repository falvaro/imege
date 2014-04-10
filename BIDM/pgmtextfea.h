#ifndef _PGMTEXTFEA_
#define _PGMTEXTFEA_

#include <limits.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "libpgm.h"

unsigned char** pgmtextfea(char *path, int ***imgh, int ***imgv, int *X, int *Y, int size);

#endif
