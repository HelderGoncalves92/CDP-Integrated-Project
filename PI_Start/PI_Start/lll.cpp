//
//  lll.cpp
//  PI_Start
//
//  Created by Hélder Gonçalves on 07/04/15.
//  Copyright (c) 2015 Hélder Gonçalves. All rights reserved.
//

#include "lll.h"

double **baseORT;

//Alocate memory to (dim+1) vectors because of BKZ insertion
void initStructsLLL(int dimension){
    int i;
    
    //Define global variable
    dim=dimension;
    
    //Allocate all necessary memory (Orthogonal Basis and its coefficients)
    baseORT = (double**)_mm_malloc((dim+1)*sizeof(double*),64);
    mu = (double**)_mm_malloc((dim+1)*sizeof(double*),64);
    
    for(i=0; i<=dim; i++){
        baseORT[i] = (double*)_mm_malloc(dim*sizeof(double),64);
        mu[i] = (double*)_mm_malloc(dim*sizeof(double),64);
    }
    
    //Quadratics Norms
    B = (double*)_mm_malloc((dim+1)*sizeof(double),64);
}

//Insert B[kl] right before B[k] and shift the other pointers
void shiftVector(long** base, int k, int kl){
    
    long *aux, *aux2;
    
    //Just if they are differents
    if(k!=kl){
        aux2 = base[kl];
    
        for (; k<kl; k++) {
            aux = base[k];
            base[k] = aux2;
            aux2 = aux;
        }
        base[k]=aux2;
    }
}

void copyVectorToDouble(double* dst, long* src){
    int i;
    for (i=0; i<dim; i++) {
        dst[i]=(double)src[i];
    }
}


//Compute all Coefficients and Norms accordingly the passed basis
void computeGSO(long** base){
    int i, j, k;
    
    //Prepare first vector
    copyVectorToDouble(baseORT[0], base[0]);
    B[0] = innerProduct(baseORT[0], baseORT[0], dim);   //<bi,bi> equals to ||bi||^2
    
    for(i=1; i<dim; i++){
        copyVectorToDouble(baseORT[i], base[i]);
        
        for(j=0; j<i; j++){
            mu[i][j] = innerProductv2(base[i], baseORT[j], dim) / B[j];
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


void sizeReduction(long** base, int k){
    int i, j;
    double r;
    
    for(i=k-1; i>=0; i--){
        //Round mu[k][i] once
        r = round(mu[k][i]);
        for(j=0; j<dim; j++){
            base[k][j] -= r * base[i][j];
        }
        
        for(j=0; j<i; j++){
            mu[k][j] -= r * mu[i][j];
        }
    }
    
    //Update the GSO accordingly the new basis
    computeGSO(base);
}

void lll(long** base, double delta, int kmax){
    
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

        //Update the GSO accordingly the new basis
        computeGSO(base);
        
        k++;
    }
}