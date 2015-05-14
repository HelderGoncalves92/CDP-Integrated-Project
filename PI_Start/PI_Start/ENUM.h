//
//  ENUM.h
//  PI_Start
//
//  Created by João Magalhães on 12/03/15.
//  Copyright (c) 2015 João Magalhães. All rights reserved.
//

#ifndef __PI_Start__ENUM__
#define __PI_Start__ENUM__

#include <xmmintrin.h>
#include "lll.h"

extern int dim;
extern double **mu;
extern double *B;

/******* STRUCTS *********/
typedef struct sEnum{
    double *cT, *y;
    int *uT, *d, *delta, *v, bound;
    struct sEnum *next;
}*Enum, NEnum;


typedef struct slist{
    Enum head;
    Enum tail;
    int count;
}*LEnum,NLEnum;


LEnum initEnum(int bound);
int* EnumSET(Enum set);


#endif /* defined(__PI_Start__ENUM__) */