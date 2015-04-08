//
//  lll.cpp
//  PI_Start
//
//  Created by Hélder Gonçalves on 07/04/15.
//  Copyright (c) 2015 Hélder Gonçalves. All rights reserved.
//

#include "lll.h"

double **baseORT;

void initStructsLLL(int dimension){
    int i;
    dim=dimension;
    
    //Allocate all necessary memory
    baseORT = (double**)_mm_malloc((dim+1)*sizeof(double*),64);
    mu = (double**)_mm_malloc((dim+1)*sizeof(double*),64);
    
    for(i=0; i<=dim; i++){
        
        //Base Orthogonal
        baseORT[i] = (double*)_mm_malloc(dim*sizeof(double),64);
        mu[i] = (double*)_mm_malloc(dim*sizeof(double),64);
    }
    
    B = (double*)_mm_malloc((dim+1)*sizeof(double),64);
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
            mu[i][j] = innerProduct(base[i], baseORT[j], dim) / B[j];
            
            for(k=0; k<dim; k++)
                baseORT[i][k] -= mu[i][j] * baseORT[j][k];
            
        }
        B[i] = innerProduct(baseORT[i], baseORT[i], dim);
    }
}


double breakCondition(int k, int kl){
    
    int i;
    double res=B[kl];
    
    for(i=k-1; i<kl; i++){
        res += pow(mu[kl][i],2) * B[i];
    }
    
    return res;
}


void sizeReduction(double** base, int k){
    int i, j;
    double r;
    
    for(i=k-1; i>=0; i--){
        r = round(mu[k][i]);
        for(j=0; j<dim; j++){
            base[k][j] -= r * base[i][j];
        }
        
        for(j=0; j<i; j++){
            mu[k][j] -= r * mu[i][j];
        }
    }
    
    //Update the GSO accordingly
    computeGSO(base);
}

void lll(double** base, double delta, int kmax){
    
    int i, k=1, kl;
    computeGSO(base);
    
    while (k<kmax) {
        sizeReduction(base, k);
        
        kl = k;
        while(k>=1 && ( (delta*B[k-1]) >= breakCondition(k, kl)) ){
            k--;
        }
        
        for(i=0; i<k; i++){
            mu[k][i] = mu[kl][i];
        }
        
        //Shift vectors
        shiftVector(base, k, kl);
        
        k++;
    }
}