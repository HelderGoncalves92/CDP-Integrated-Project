//Algoritmo em http://www.csie.nuk.edu.tw/~cychen/Lattices/Lattice%20Basis%20Reduction_%20Improved%20Practical%20Algorithms%20and%20Solving%20Subset%20Sum%20Problems.pdf - pagina 16

#include "ENUM.h"
#include "math.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

int *u, *toInit_0,*toInit_1,*toInit_d0;
double cL;
LEnum list = NULL;

Enum newEnumElem(int bound, int simbling){
    int dimesion = bound+1;
    Enum st = (Enum)_mm_malloc(sizeof(NEnum),64);
    
    st->cT = (double*)_mm_malloc(dimesion * sizeof(double), 64);
    st->y =  (double*)_mm_malloc(dimesion * sizeof(double), 64);
    st->v =     (int*)_mm_malloc(dimesion * sizeof(int), 64);
    st->delta = (int*)_mm_malloc(dimesion * sizeof(int), 64);
    st->d =     (int*)_mm_malloc(dimesion * sizeof(int), 64);
    st->uT =    (int*)_mm_malloc(dimesion * sizeof(int), 64);
    
    //Fill vectors
    memcpy(st->v, toInit_0, dimesion*sizeof(int));
    memcpy(st->delta, toInit_0, dimesion*sizeof(int));
    memcpy(st->uT, toInit_0, dimesion*sizeof(int));
    memcpy(st->d, toInit_1, dimesion*sizeof(int));
    
    memcpy(st->y, toInit_d0, dimesion*sizeof(double));
    memcpy(st->cT, B, dimesion*sizeof(double));

    //Prepare to start
    st->cT[bound]=0.0;
    st->next = NULL;
    st->bound = bound-1;
    st->d[bound-1] = st->uT[bound-1] = 1;
    
    return st;
}

void printVec(Enum st){
    int i;
    
    int dimension = st->bound+1;

    printf("uT:");
    for(i=0; i<dimension; i++)
        printf("%d ",st->uT[i]);
    printf("\n");
    
    printf("d:");
    for(i=0; i<dimension; i++)
        printf("%d ",st->d[i]);
    printf("\n");
    
    printf("delta:");
    for(i=0; i<dimension; i++)
        printf("%d ",st->delta[i]);
    printf("\n");

    printf("v:");
    for(i=0; i<dimension; i++)
        printf("%d ",st->v[i]);
    printf("\n");
    
    printf("cT:");
    for(i=0; i<dimension; i++)
        printf("%.2f ",st->cT[i]);
    printf("\n");
    printf("Y:");
    for(i=0; i<dimension; i++)
        printf("%.2f ",st->y[i]);
    printf("\n");

        printf("\n");
}


void initEnum(){
    
    //Allocate memory
    u = (int*)_mm_malloc(dim*sizeof(NLEnum), 64);
    toInit_0 = (int*)_mm_malloc((dim+1)*sizeof(int), 64);
    toInit_1 = (int*)_mm_malloc((dim+1)*sizeof(int), 64);
    toInit_d0 = (int*)_mm_malloc((dim+1)*sizeof(double), 64);
    
    list = (LEnum)_mm_malloc(sizeof(NLEnum), 64);
    list->count = 0;
    list->head = list->tail = NULL;;
    
    //Fill toInit with '0' or '1'
    for(int i=0; i<=dim; i++){
        toInit_0[i] = 0;
        toInit_1[i] = 0;
        toInit_d0[i] = 0;
    }
    
    memcpy(&u[0], toInit_0, dim*sizeof(int));
    u[0]=1;
    cL = B[0];
}


//ENUM accordingly C. P. Schnorr && M. Euchner
void EnumSET(Enum set){
    
    double aux;
    int bound = set->bound;
    int s = bound, t = bound, i;
    bound++;
    
    set->cT[t] = set->cT[t + 1] + (set->y[t]*set->y[t] + 2*set->uT[t]*set->y[t] + set->uT[t]*set->uT[t]) * B[t];
    int cic=0;
	while(t < bound){
        cic++;
        //printVec(set);
        
        if (set->cT[t] < cL){
			if (t > 0){
                //moveDown
         //       printf("%d:DOWN\n",t);
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
                printf("%d:UPDATE\n",s);
                //printVec(set);
                
				cL = set->cT[0];
                memcpy(&u[0], set->uT, bound*sizeof(int));
                memcpy(&u[bound], toInit_0, (dim-bound)*sizeof(int));
			}
		}
		else{
            //moveUp
        //    printf("%d:UP\n",t);
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
    printf("cic:%d\n",cic);
}


int* ENUM(int nthreads){
    
    Enum set = newEnumElem(49, 0);
    
    EnumSET(set);
    
    return u;
}

