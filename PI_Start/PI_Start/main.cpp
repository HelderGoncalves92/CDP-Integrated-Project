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




void MatIntFromMatZZ(double** A, NTL::mat_ZZ B) {
    
    int row,col, rows = B.NumRows(), cols= B.NumCols();
    
    for (row = 0; row < rows; ++row) {
        for (col = 0; col < cols; ++col) {
            A[row][col] = to_long(B[row][col]);
            
        }
    }
}

int main(int argc, const char * argv[]) {
    
    
    NTL::mat_ZZ B;
    std::ifstream input_file(argv[1]);
    
    //Try to read Base file
    if (input_file.is_open()) {
        input_file >> B;
        input_file.close();
    }
    else{
        printf("File was not read\n");
        exit(-1);
    }
    
    
    //NTL::G_BKZ_FP(B, 0.99, 20 ); //BKZ janela 20
    //NTL::G_LLL_FP(B,0.99);
    
    //Basis Matrix - Memory Allocation
    int rows = (int)B.NumRows(), cols= (int)B.NumCols();
    double** B_ = (double**)_mm_malloc((rows+1) *sizeof(double*),64);
    
    for(int row = 0; row < rows+1; row++)
        B_[row] = (double*)_mm_malloc(cols*sizeof(double),64);
    
    
    //Convert ZZ data to double
    MatIntFromMatZZ(B_, B);
    
    //Init all Structs
    initBKZ(cols);
    
    BKZ(B_, 20, 0.99);
    //lll(B_, 0.99, dim);
    
    
    cout << "********************" << endl;
    for (int i=0; i<rows; i++) {
        for (int j=0; j<cols; j++) {
            cout << B_[i][j]<< " ";
        }
        cout << endl;
    }
    
    
    cout << "Finish" << endl;
    return 0;
}
