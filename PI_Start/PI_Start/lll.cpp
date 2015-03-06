//
//  lll.cpp
//  PI_Start
//
//  Created by Hélder Gonçalves on 06/03/15.
//  Copyright (c) 2015 Hélder Gonçalves. All rights reserved.
//

#include "lll.h"

void coefficientsGS(double *base[], double *baseORT[], double* coeffs[], int dimVector, int numVector){
    
    int i, j;
    
    for(i=numVector-1; i>0; i--){
        for (j=i-1; j>=0; j--) {
            coeffs[i][j]= innerProduct(&base[i][0], &baseORT[j][0], dimVector) / innerProduct(&baseORT[j][0], &baseORT[j][0], dimVector);
        }
    }
}


void gram_Schmidt_Orthonormalization(double* base[], int dimVector, int numVector){
    
    int i, k, j;
    double *projAux = (double*)malloc(dimVector*sizeof(double));
    
    for(k=0; k<numVector; k++){
        
        normalizeVector(&base[k][0], dimVector);
        
        for(i=k+1; k<numVector; k++){
            //remove component in direction Vk
            projection(&base[k][0], &base[i][0], projAux, dimVector);
            
            for(j=0; j<dimVector; j++)
                base[k][j] -= projAux[j];
            
        }
    }
}


