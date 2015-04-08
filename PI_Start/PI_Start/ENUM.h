//
//  ENUM.h
//  PI_Start
//
//  Created by João Magalhães on 12/03/15.
//  Copyright (c) 2015 João Magalhães. All rights reserved.
//

#ifndef __PI_Start__ENUM__
#define __PI_Start__ENUM__

#include <xmmintrin.h>
#include "lll.h"

extern int dim;
extern double **mu;
extern double *B;

void initENUM();
int* ENUM(int ini, int fim);

int* basicENUM();

#endif /* defined(__PI_Start__ENUM__) */