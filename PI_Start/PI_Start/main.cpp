//
//  main.cpp
//  PI_Start
//
//  Created by Hélder Gonçalves on 06/03/15.
//  Copyright (c) 2015 Hélder Gonçalves. All rights reserved.
//

#include <NTL/ZZ.h>
#include <NTL/vec_double.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/LLL.h>

#include <iostream>
#include <fstream>
#include "BKZ.h"

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
    NTL::G_LLL_FP(BB,0.99);
    
    //Basis Matrix - Memory Allocation
    int rows = (int)BB.NumRows(), cols= (int)BB.NumCols();
    long** BB_ = (long**)_mm_malloc((rows+1) *sizeof(long*),64);
    
    for(int row = 0; row < rows+1; row++)
        BB_[row] = (long*)_mm_malloc(cols*sizeof(long),64);
    
    
    //Convert ZZ data to double
    MatIntFromMatZZ(BB_, BB);
    
    //Init all Structs (Vectors an Matrix)
    long* fvec = (long*)calloc(cols ,sizeof(long));
    initBKZ(cols);
    
    computeGSO(BB_);
    int* vec = ENUM(0, dim-1);
    
    computeNewVector(fvec, vec, BB_);
    
    //FINAL OUTPUT
   
  /*  for (int i=0; i<rows; i++) {
        for (int j=0; j<cols; j++) {
            cout << BB_[i][j]<< " ";
        }
        cout << endl;
    }
    */
    
    cout << "\n********************" << endl;
    for (int i=0; i<rows; i++)
        cout << vec[i]<< " ";
    cout << "\n********************" << endl;
    for (int i=0; i<rows; i++)
        cout << fvec[i]<< " ";
    cout << "Finish" << endl;
    return 0;
}
