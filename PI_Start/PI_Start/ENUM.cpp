//Algoritmo em http://www.csie.nuk.edu.tw/~cychen/Lattices/Lattice%20Basis%20Reduction_%20Improved%20Practical%20Algorithms%20and%20Solving%20Subset%20Sum%20Problems.pdf - pagina 16

#include "ENUM.h"
#include "math.h"
#include <iostream>
#include <stdio.h>

double cL, *cT, *y;
int *u, *uT, *d, *v, *delta;



void initENUM(){
    //Alocate memory for every vector
    cT = (double*)_mm_malloc((dim+1) * sizeof(double), 64);
    y = (double*)_mm_malloc((dim+1) * sizeof(double), 64);
    v = (int*)_mm_malloc((dim+1) * sizeof(int), 64);
    delta = (int*)_mm_malloc((dim+1) * sizeof(int), 64);
    d = (int*)_mm_malloc((dim+1) * sizeof(int), 64);
    u = (int*)_mm_malloc((dim+1) * sizeof(int), 64);
    uT = (int*)_mm_malloc((dim+1) * sizeof(int), 64);
}


//ENUM accordingly C. P. Schnorr && M. Euchner
int* ENUM(int ini, int fim){
    int s = ini, t = ini, i;
    
    for(i=1; i<dim; i++){
        for(int j=0; j<i; j++){
            printf("MU:%f\n", mu[i][j]);
        }
        printf("B:%f\n", B[i]);
    }
    
    //Init all vectors
	cL = B[ini];
	u[ini] = uT[ini] = 1;
	y[ini] = delta[ini] = v[ini]= 0;
	d[ini] = 1;

	for(i = ini + 1; i <= fim+1; i++){
		cT[i] = u[i] = uT[i] = y[i] = delta[i] = v[i] = 0;
        v[i] = 0;
		d[i] = 1;
    }
    
    printf("%d %d\n",t,fim);
	while(t <= fim){
        
        printf("***STAR_IT***\nT:: %d\ncT[t+1]: %f\ny[t]: %f\nuT[t]: %d\nB[t]: %f\n",t,cT[t + 1],y[t],uT[t],B[t]);
        
		cT[t] = cT[t + 1] + pow(y[t] + double(uT[t]),2) * B[t];
        
        printf("Res:cT[t]: %f\n cL: %f\n",cT[t],cL);
		if (cT[t] < cL){
            printf("True1:t: %d\n",t);
			if (t > ini){
				t--;
                printf("True2:t: %d\n",t);
				y[t] = 0;
				for (i = t + 1; i <= s; i++){
					y[t] += uT[i] * mu[i][t];
				}
                printf("True2:Res:y[t] %f\n",y[t]);
                
				uT[t] = v[t] = int(round(-y[t]));
                printf("round:y[t]: %d\n",uT[t]);
				delta[t] = 0;

				if (uT[t] > -y[t]){
                    printf("TRUE\n");
					d[t] = -1;
				}
				else{
                    printf("FALSE\n");
					d[t] = 1;
				}
			}
			else{
                //printf("%f | %f -||- %d | %d\n",cT[t], cL,t,ini);
                printf("Else2\n");
				cL = cT[ini];
				for (i = ini; i <= fim; i++){
					u[i] = uT[i];
				}
			}
		}
		else{
			t++;
            
			s = maxn(s,t);
            printf("Else1:S:%d\n",s);
            
			if(t < s){
                printf("Else1:DELTA_A:%d\n",delta[t]);
				delta[t] = -delta[t];
                printf("Else1:DELTA_D:%d\n",delta[t]);
			}
			if(delta[t]*d[t] >= 0){
                printf("Else1:D[T]:%d\n",d[t]);
				delta[t] += d[t];
                
			}
            
            printf("Else1:V[T]:%d\n",v[t]);
            
			uT[t] = v[t] + delta[t];
            printf("Else1:uT[t]:%d\n",uT[t]);
		}
	}
    
    // printf("%d | %d | %d\n",a,b,c);
	return u;
}














