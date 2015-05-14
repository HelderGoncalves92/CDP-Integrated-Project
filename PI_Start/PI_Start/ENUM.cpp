//Algoritmo em http://www.csie.nuk.edu.tw/~cychen/Lattices/Lattice%20Basis%20Reduction_%20Improved%20Practical%20Algorithms%20and%20Solving%20Subset%20Sum%20Problems.pdf - pagina 16

#include "ENUM.h"
#include "math.h"
#include <iostream>
#include <stdio.h>

int *u;
double cL;
LEnum list = NULL;

Enum newEnumElem(int bound){
    int dimesion = dim+1;
    Enum st = (Enum)_mm_malloc(sizeof(NEnum),64);
    st->next = NULL;
    st->bound = bound;
    
    st->cT = (double*)_mm_malloc(dimesion * sizeof(double), 64);
    st->y =  (double*)_mm_malloc(dimesion * sizeof(double), 64);
    st->v =     (int*)_mm_malloc(dimesion * sizeof(int), 64);
    st->delta = (int*)_mm_malloc(dimesion * sizeof(int), 64);
    st->d =     (int*)_mm_malloc(dimesion * sizeof(int), 64);
    st->uT =    (int*)_mm_malloc(dimesion * sizeof(int), 64);

    return st;
}

void initEnumElem(Enum st){
    int i;
    
    //Init all vectors
    st->d[0] = st->uT[0] = 1;
    st->y[0] = st->delta[0] = st->v[0] = 0;
    
    for(i = 1; i <= st->bound; i++){
        st->cT[i] = st->uT[i] = st->y[i] = st->delta[i] = st->v[i] = 0;
        st->d[i] = 1;
    }

}

LEnum initEnum(int bound){
    u = (int*)_mm_malloc(dim * sizeof(int), 64);
    cL = B[0];
    
    Enum first = newEnumElem(bound);
    initEnumElem(first);
    
    list = (LEnum)_mm_malloc(sizeof(NLEnum), 64);
    
    list->count = 1;
    list->head = list->tail = first;
    
    return list;
    
}


//ENUM accordingly C. P. Schnorr && M. Euchner
int* EnumSET(Enum set){
    
    double aux;
    int s = 0, t = 0, i;
    int bound = set->bound;
    
    
   set->cT[t] = set->cT[t + 1] + (set->y[t]*set->y[t] + 2*set->uT[t]*set->y[t] + set->uT[t]*set->uT[t]) * B[t];
    
	while(t < bound){
        printf("%d\n",t);
        if (set->cT[t] < cL){
			if (t > 0){
                //moveDown
				t--;
                aux = 0;
                
				for (i = t + 1; i <= s; i++){
					aux += set->uT[i] * mu[i][t];
				}
                
                set->y[t] = aux;
				set->uT[t] = set->v[t] = int(round(-aux));
                set->delta[t] = 0;
                
				if (set->uT[t] > -aux){
					set->d[t] = -1;
				}
				else{
					set->d[t] = 1;
				}
                
                //Prepare cT[t] to next iteration
                set->cT[t] = set->cT[t + 1] + (set->y[t]*set->y[t] + 2*set->uT[t]*set->y[t] + set->uT[t]*set->uT[t]) * B[t];
                
			}
			else{
                //UpdateVector
				cL = set->cT[0];
                memcpy(&u[0], set->uT, bound*sizeof(int));
			}
		}
		else{
            //moveUp
			t++;
            
			s = (s<t)?t:s; //Get max value
            
			if(t < s){
				set->delta[t] = -set->delta[t];
			}
			if(set->delta[t]*set->d[t] >= 0){
				set->delta[t] += set->d[t];
			}
			set->uT[t] = set->v[t] + set->delta[t];
            
            //Prepare cT[t] to next iteration
            set->cT[t] = set->cT[t + 1] + (set->y[t]*set->y[t] + 2*set->uT[t]*set->y[t] + set->uT[t]*set->uT[t]) * B[t];
            
		}
	}
	return u;
}


void ENUM(int nthreads, int min){
    
    
    
    
}

