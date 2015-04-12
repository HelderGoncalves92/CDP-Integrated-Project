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
 #include <xmmintrin.h>


 
extern int dim;
extern double **mu;
extern double *B;
 
void initStructsLLL(int dim);
void computeGSO(long** base);
void shiftVector(long** base, int k, int kl);
void lll(long** base, double delta, int kmax);
 

#endif /* defined(__PI_Start__lll__) */
