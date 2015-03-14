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
#include "lll.h"

extern int dimVector;
extern int numVector;




void MatIntFromMatZZ(long** A, NTL::mat_ZZ B) {
    
    int row,col, rows = B.NumRows(),cols= B.NumCols();
    
    for (row = 0; row < rows; ++row) {
        for (col = 0; col < cols; ++col) {
            A[row][col] = to_long(B[row][col]);
            
        }
    }
}

int main(int argc, const char * argv[]) {
    
    
    NTL::mat_ZZ B;
    std::ifstream input_file(argv[1]);
    if (input_file.is_open()) {
        input_file >> B;
        input_file.close();
    }
    else{
        printf("File was not read\n");
        exit(-1);
    }
    
    
    NTL::G_BKZ_FP(B, 0.99, 10 );
    NTL::G_LLL_FP(B,0.99);
    
    
    int rows = B.NumRows(), cols= B.NumCols(), i;
    long** B_ = (long**)malloc(B.NumCols()*sizeof(long*));
    
    for(int row = 0; row < rows; row++)
        B_[row] = (long*)malloc(B.NumRows()*sizeof(long));
    
    
    MatIntFromMatZZ(B_, B);
    
    
    return 0;
}
