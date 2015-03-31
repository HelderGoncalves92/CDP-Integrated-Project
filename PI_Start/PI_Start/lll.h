//
//  lll.h
//  PI_Start
//
//  Created by Hélder Gonçalves on 06/03/15.
//  Copyright (c) 2015 Hélder Gonçalves. All rights reserved.
//

#ifndef __PI_Start__lll__
#define __PI_Start__lll__

#include "simpleMath.h"

extern int dimVector;
extern int numVector;


void shiftVector(double** base, int k, int kl);
void initStructsLLL(double** base, int rows, int cols);
void lll(double** base, int delta, int indiceMax);



#endif /* defined(__PI_Start__lll__) */
