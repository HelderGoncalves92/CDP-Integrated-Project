#include "ENUMwPrunning.h"
#include "xmmintrin.h"

int* EnumWPrun (double* b, int ini, int fim){
	double C = b[ini];
	//O i começa a ini pois e o primeiro indice
	int i = ini, j, k, tam = fim - ini;
	double *dist, *c;
	double **E;
	int *res, *delta, *d, *u, *uL;
	int last_nonzero = 0;

	//Alocate vectors
	dist = (double*)_mm_malloc((dim + 1) * sizeof(double), 64);
	c = (double*)_mm_malloc(dim*sizeof(double), 64);
	delta = (int*)_mm_malloc(dim*sizeof(int), 64);
	d = (int*)_mm_malloc(dim*sizeof(int), 64);
	u = (int*)_mm_malloc(dim*sizeof(int), 64);
	uL = (int*)_mm_malloc(dim*sizeof(int), 64);
	res = (int*)_mm_malloc(dim*sizeof(int), 64);
	E = (double**)_mm_malloc((dim + 1)*sizeof(double*), 64);

	for (j = 0; j <= dim; j++){
		E[j] = (double*)_mm_malloc(dim*sizeof(double), 64);
		for (k = 0; k < dim; k++){
			E[j][k] = 0;
		}
	}
	for (j = 0; j < dim; j++){
		res[j] = 0;
	}

	u[ini] = uL[ini] = 1;

	for(j = ini; j <= fim+1; j++){
		u[j] = 0;
		uL[j] = 0;
		dist[j] = 0;
	}
	for(j = ini; j <= fim; j++){
		c[j] = 0;
		delta[j] = 0;
		d[j] = j + 1 - ini;
	}

	while(true){
		dist[i] = dist[i + 1] + pow((u[i] - c[i]),2) * b[i];
		if(dist[i] < C){
			if(i != 0){
				//move down
				i--;
				d[i - 1] = fmax(d[i - 1], d[i]);
				for(j = d[i]; j <= i+1; j--){
					E[j][i] = E[j+1][i] + u[j] * mu[j][i];
				}
				c[i] = -E[i + 1][i];
				u[i] = round(c[i]);
				delta[i] = 1;
			}else{
				//update best vector
				C = dist[i];
				for(j = 0; j < dim; j++){
					uL[j] = u[j];
				}
			}
		}else{

			if(i == tam){
				for (j = ini; j <= fim; j++){
					res[j] = uL[j];
				}
				return res;
			}
			//move up
			i++;
			d[i - 1] = i - ini + 1;
			if(i - ini >= last_nonzero){
				last_nonzero = i - ini;
				u[i]++;
			}else{
				if(u[i] > c[i]){
					u[i] = u[i] - delta[i];
				}else{
					u[i] = u[i] + delta[i];
				}
				delta[i]++;
			}
		}
	}

}

