//
//  lll.cpp
//  PI_Start
//
//  Created by Hélder Gonçalves on 06/03/15.
//  Copyright (c) 2015 Hélder Gonçalves. All rights reserved.
//

#include "lll.h"


void allocateMemory(double *base[], double *baseORT[], double* coeffs[], double* norms, int dimVector, int numVector){
    int i;
    
    baseORT = (double**)malloc(numVector*sizeof(double*));
    coeffs = (double**)malloc(numVector*sizeof(double*));
    
    for(i=0; i<numVector; i++){
        
        //Base Orthogonal
        baseORT[i] = (double*)malloc(dimVector*sizeof(double));
        memccpy(baseORT[i], base[i], dimVector, sizeof(double));
        
        coeffs[i] = (double*)malloc(dimVector*sizeof(double));
    }
    
    norms = (double*)malloc(dimVector*sizeof(double));
}



void coefficientsGS(double *base[], double *baseORT[], double* coeffs[],double* norms){
    
    int i, j;
    
    for(i=numVector-1; i>0; i--){
        for (j=i-1; j>=0; j--) {
            coeffs[i][j]= innerProduct(&base[i][0], &baseORT[j][0], dimVector) / innerProduct(&baseORT[j][0], &baseORT[j][0], dimVector);
        }
        
    }
    
    for (i=0; i<dimVector; i++) {
        norms[i] = pow(vectorNorm(&baseORT[i][0], dimVector), 2);
    }
}


void gram_Schmidt_Orthonormalization(double* base[], double* baseORT[]){
    
    int i, k, j;
    
    //Allocate Auxmemory
    double *projAux = (double*)malloc(dimVector*sizeof(double));
    
    for(k=0; k<numVector; k++){
        
        normalizeVector(baseORT[k], dimVector);
        
        for(i=k+1; k<numVector; k++){
            //remove component in direction Vk
            projection(baseORT[k], baseORT[i], projAux, dimVector);
            
            for(j=0; j<dimVector; j++)
                baseORT[k][j] -= projAux[j];
            
        }
    }
    free(projAux);
}

double breakCondition_Alg2(double **coeffs, double* norms, int k, int kl){
    
    int i;
    double res=norms[kl];
    
    for(i=k-1; i<= kl-1; i++ ){
        res += pow(coeffs[kl][i],2) * norms[i];
    }
    
    return res;
}



//Algorithm 1
void sizeReduction(double* base[], double* coeffs[], int k){
    int i, j;
    
    for(i=k-1; i>=0; i--){
        
        for(j=0; j<dimVector; j++){
            base[k][j] -= round(coeffs[k][i]) * base[i][j];
        }
        
        for(j=0; j<i; j++){
            coeffs[k][j] -= round(coeffs[k][i])*coeffs[i][j];
        }
    }
}

//Algorithm 2 - Lenstra–Lenstra–Lovász
void lll(double* base[], int delta){
    
    int i, k=1, kl;
    double **baseORT, **coeffs, *norms;
    double *aux, *aux2;
    
    
    //Initializate data structs
    allocateMemory(base, baseORT, coeffs, norms, dimVector, numVector);
    
    //Compute orthogonal basis, Coefficients, norms
    gram_Schmidt_Orthonormalization(base, baseORT);
    coefficientsGS(base, baseORT, coeffs, norms);
    
    while(k<numVector){
        //sizeResuction and Update the Coeffs and norms
        sizeReduction(baseORT, coeffs, k);
        coefficientsGS(base, baseORT, coeffs, norms);

        kl=k;
        while(k>=1 && (delta*norms[k-1] > breakCondition_Alg2(coeffs, norms, k, kl)) ){
            k--;
        }
        
        for(i=0; i<k-1; i++){
            coeffs[k][i] = coeffs[kl][i];
        }
        
        
        //insert B[kl] right before B[k]
        aux=base[k+1];
        base[k+1]=base[kl];
        
        //Correct the other pointers
        for (i=k+2; i<=kl; i++) {
            aux2 = base[k];
            base[k] = aux;
            aux = aux2;
        }
        
        k++;
    }
    
}









