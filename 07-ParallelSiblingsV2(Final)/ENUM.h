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


extern short dim;
extern double **mu;
extern double *B;

/******* STRUCTS *********/
typedef struct sEnum{
    short bound, type, sibling, level, *vec;
    struct sEnum *next;
}*Enum, NEnum;


typedef struct slist{
    Enum head;
    Enum tail;
    short count;
}*LEnum,NLEnum;


void initEnum(short nthreads);
short* ENUM();


#endif /* defined(__PI_Start__ENUM__) */