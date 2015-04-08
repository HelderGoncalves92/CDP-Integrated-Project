//
//  lll.cpp
//  PI_Start
//
//  Created by Hélder Gonçalves on 07/04/15.
//  Copyright (c) 2015 Hélder Gonçalves. All rights reserved.
//

#include "lll.h"

double **baseORT;
double **u;

double *c;
double *B;

void initStructsLLL(int dimension){
    int i;
    dim=dimension;
    
    //Allocate all necessary memory
    baseORT = (double**)_mm_malloc(dim*sizeof(double*),64);
    u = (double**)_mm_malloc(dim*sizeof(double*),64);
    
    for(i=0; i<dim; i++){
        
        //Base Orthogonal
        baseORT[i] = (double*)_mm_malloc(dim*sizeof(double),64);
        u[i] = (double*)_mm_malloc(dim*sizeof(double),64);
    }
    
    c = (double*)_mm_malloc(dim*sizeof(double),64);
    B = (double*)_mm_malloc(dim*sizeof(double),64);
}


void shiftVector(double** base, int k, int kl){
    
    double *aux, *aux2;
    
    if(k!=kl){
        aux2 = base[kl];
    
        //insert B[kl] right before B[k] and shift the other pointers
        for (; k<kl; k++) {
            aux = base[k];
            base[k] = aux2;
            aux2 = aux;
        }
        base[k]=aux2;
    }
}

void computeGSO(double** base){
    
    int i, j, k;
    
    memcpy(baseORT[0], base[0], dim*sizeof(double));
    B[0] = innerProduct(baseORT[0], baseORT[0], dim);
    
    for(i=1; i<dim; i++){
        memcpy(baseORT[i], base[i], dim*sizeof(double));
        
        for(j=0; j<i; j++){
            u[i][j] = innerProduct(base[i], baseORT[j], dim) / B[j];
            
            for(k=0; k<dim; k++)
                baseORT[i][k] -= u[i][j] * baseORT[j][k];
            
        }
        B[i] = innerProduct(baseORT[i], baseORT[i], dim);
    }
    
    for (i=0; i<dim; i++) {
        c[i] = pow(vectorNorm(baseORT[i], dim), 2);
    }
}


double breakCondition(int k, int kl){
    
    int i;
    double res=c[kl];
    
    for(i=k-1; i<kl; i++){
        res += pow(u[kl][i],2) * c[i];
    }
    
    return res;
}


void sizeReduction(double** base, int k){
    int i, j;
    
    for(i=k-1; i>=0; i--){
        
        for(j=0; j<dim; j++){
            base[k][j] -= round(u[k][i]) * base[i][j];
        }
        
        for(j=0; j<i; j++){
            u[k][j] -= round(u[k][i]) * u[i][j];
        }
    }
    
    //Update the GSO accordingly
    computeGSO(base);
}

void lll(double** base, double delta){
    
    int i, k=1, kl;
    
    computeGSO(base);
    
    while (k<dim) {
        sizeReduction(base, k);
        
        kl = k;
        while(k>=1 && ( (delta*c[k-1]) >= breakCondition(k, kl)) ){
            k--;
        }
        
        for(i=0; i<k-1; i++){
            u[k][i] = u[kl][i];
        }
        
        //Shift vectors
        shiftVector(base, k, kl);
        
        k++;
    }
}