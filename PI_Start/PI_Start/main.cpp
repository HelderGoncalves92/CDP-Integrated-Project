//
//  main.cpp
//  PI_Start
//
//  Created by Hélder Gonçalves on 06/03/15.
//  Copyright (c) 2015 Hélder Gonçalves. All rights reserved.
//
#include <stdio.h>
#include <NTL/ZZ.h>
#include <NTL/vec_double.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/LLL.h>

#include <iostream>
#include <fstream>
#include <time.h>
#include <omp.h>

#include "BKZ.h"

#define NUM_EVENTS 2

using namespace std;
int dim;
double **mu;
double *B;

void computeNewVector(long* dst_vec, int* src_vec ,long** base){
    int i,l;
    
    for (i = 0; i < dim; i++){
        for (l = 0; l < dim; l++){
            dst_vec[l] += -src_vec[i] * base[i][l];
        }
    }
}

void MatIntFromMatZZ(long** A, NTL::mat_ZZ B) {
    
    int row,col, rows = (int)B.NumRows(), cols= (int)B.NumCols();
    
    for (row = 0; row < rows; ++row) {
        for (col = 0; col < cols; ++col) {
            A[row][col] = to_long(B[row][col]);
            
        }
    }
}

void check_equals(int* vec, int* ntl, int rows){
    int mult = 1, equal = 1, i;
    if(vec[0] == -ntl[0]){
	mult = -1;
    }
    for(i = 0;i < rows;i++){
	if(vec[i] != mult*ntl[i]){equal=0;}
    }
    if(equal==1){
	printf("Vectors match\n");
    }else{
	printf("WARNING vectors do not match\n");
    }
} 

int main(int argc, const char * argv[]) {

    NTL::mat_ZZ BB;
    std::ifstream input_file(argv[1]);
    
    //Try to read Base file
    if (input_file.is_open()) {
        input_file >> BB;
        input_file.close();
    }
    else{
        printf("File was not read\n");
        exit(-1);
    }
    
    //To call LLL or BKZ of NTL
    NTL::G_BKZ_FP(BB, 0.99, 20 ); //BKZ janela 20

    
    //Basis Matrix - Memory Allocation
    int rows = (int)BB.NumRows(), cols= (int)BB.NumCols();
    long** BB_ = (long**)_mm_malloc((rows+1) *sizeof(long*),64);
    
    for(int row = 0; row < rows+1; row++)
        BB_[row] = (long*)_mm_malloc(cols*sizeof(long),64);
    
    
    //Convert ZZ data to double
    MatIntFromMatZZ(BB_, BB);
    
    //Init all Structs (Vectors an Matrix)
    long* fvec = (long*)calloc(cols ,sizeof(long));
    initStructsLLL(cols);
    initENUM();
    
    //Compute all Coefficients and Norms accordingly the basis
    computeGSO(BB_);

    double time = omp_get_wtime();
    
    int* vec = ENUM(0, dim-1);
    
    time = omp_get_wtime() -time;
    
    //Compute SVP with ENUM vector
    computeNewVector(fvec, vec, BB_);
    
    double norm = sqrt(innerProductv3(fvec, fvec, dim));
    
    //Final Outputs
    cout << "****** ENUM VEC ******" << endl;
    for (int i=0; i<rows; i++)
        cout << vec[i]<< " ";
    cout << "\n******** SVP ********" << endl;
    for (int i=0; i<rows; i++)
        cout << fvec[i]<< " ";

    
    cout << "\n******* NORM *******" << endl;
    cout << norm;
    cout  << endl;
    
    cout << "Time: " << time << endl;
    return 0;
}
