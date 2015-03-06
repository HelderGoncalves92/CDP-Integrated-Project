//
//  simpleMath.h
//  PI_Start
//
//  Created by Hélder Gonçalves on 06/03/15.
//  Copyright (c) 2015 Hélder Gonçalves. All rights reserved.
//

#ifndef __PI_Start__simpleMath__
#define __PI_Start__simpleMath__

#include <iostream>
#include <math.h>


double vectorNorm(double *v, int dim);
void normalizeVector(double *v, int dim);
double innerProduct(double *x, double *y, int dim);
void projection(double *u, double *v, double *res, int dim);


#endif /* defined(__PI_Start__simpleMath__) */
