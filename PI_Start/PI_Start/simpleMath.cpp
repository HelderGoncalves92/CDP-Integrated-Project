//
//  simpleMath.cpp
//  PI_Start
//
//  Created by Hélder Gonçalves on 06/03/15.
//  Copyright (c) 2015 Hélder Gonçalves. All rights reserved.
//

#include "simpleMath.h"


double vectorNorm(double *v, int dim){
    int i;
    double res = 0;
    
    for(i=0; i<dim; i++)
        res += pow(v[i], 2);
    
    return ( res/dim );
}

void normalizeVector(double *v, int dim){
    
    int i=0;
    double norm = vectorNorm(v, dim);
    
    for(i=0; i<dim; i++)
        v[i] /= norm;
    
}

double innerProduct(double *x, double *y, int dim){
    
    int i;
    double result=0;
    
    for(i=0; i<dim; i++){
        result += x[i]*y[i];
    }
    return result;
}


void projection(double *u, double *v, double *res, int dim){
    
    int i;
    double scalar = innerProduct(v, u, dim) / innerProduct(u, u, dim);
    
    for(i=0; i<dim; i++){
        res[i] = u[i] * scalar;
    }
}
