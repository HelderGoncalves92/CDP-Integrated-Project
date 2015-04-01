#include "ENUMwPrunning.h"

int* EnumWPrun (double* mu[], double b[], int ini, int fim){
	double C = b[ini];
	//O i começa a ini pois e o primeiro indice
	int i = ini, j, k;
	int tam = fim - ini;
	double dist[tam+1], c[tam];
	double E[tam+1][tam];
	int delta[tam], d[tam], u[tam], uL[tam];
	int last_nonzero = 0;

	int *ret = (int*)calloc(dimVector,sizeof(int));

	u[0] = uL[0] = 1;

	for(j = 1; j < tam; j++){
		u[j] = 0;
		uL[j] = 0;
		dist[j] = 0;
	}
	for(j = 0; j < tam; j++){
		c[j] = 0;
		delta[j] = 0;
		d[j] = j + 1;
	}
	for(j = 0; j < tam + 1; j++){
		for(k = 0; k < tam; k++){
			E[j][k] = 0;
		}
	}

	while(true){
		dist[i - ini] = dist[i + 1 - ini] + pow((u[i - ini] - c[i]),2) * b[i];
		if(dist[i - ini] < C){
			if(i != 0){
				//move down
				i--;
				d[i - 1 - ini] = fmax(d[i - 1 - ini], d[i - ini]);
				for(j = d[i - ini]; j <= i+1; j--){
					E[j][i - ini] = E[j+1][i - ini] + u[j - ini] * mu[j][i];
				}
				c[i] = -E[i + 1 - ini][i - ini];
				u[i - ini] = round(c[i]);
				delta[i - ini] = 1;
			}else{
				//update best vector
				C = dist[i - ini];
				for(j = 0; j < tam; j++){
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
			d[i - 1 - ini] = i - ini + 1;
			if(i - ini >= last_nonzero){
				last_nonzero = i - ini;
				u[i - ini]++;
			}else{
				if(u[i - ini] > c[i]){
					u[i - ini] = u[i - ini] - delta[i - ini];
				}else{
					u[i - ini] = u[i - ini] + delta[i - ini];
				}
				delta[i - ini]++;
			}
		}
	}

}

