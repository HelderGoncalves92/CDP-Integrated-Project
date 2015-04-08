//
//  BKZ.h
//  PI_Start
//
//  Created by João Magalhães on 12/03/15.
//  Copyright (c) 2015 João Magalhães. All rights reserved.
//

#ifndef __PI_Start__BKZ__
#define __PI_Start__BKZ__

#include "lll.h"
#include "ENUMwPrunning.h"
#include "ENUM.h"

extern double **mu;
extern double *B;

void initBKZ(int dimension);
double** BKZ(double** bases, int beta, double delta);

#endif /* defined(__PI_Start__BKZ__) */