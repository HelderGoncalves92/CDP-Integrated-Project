//N = Dimensao, por enquanto definida globalmente para simplificar
#define N  2;

#include "BKZ.h"

bool passvec(double v[], int index){
	int i;
	for(i = 0; i < N; i++){
		if (v[i] != 0) return false;
	}
	return true;
}

double** BKZ(double *bases[], double *mu[], double *c[], int beta, double delta){
	int k, h, i, l;
	double v[N], vaux[N];
	double* aux;

	int z = 0, j = 0;

	lll(&bases, delta, N-1);
	while (z < N - 1){
		j = (j * mod(N-1)) + 1; 
		k = min(j + beta - 1, N);
		h = min(k + 1, N);
		
		v = ENUM(u,c,j-1,k-1);
		if (!passvec(v))
		{
			z = 0;
			/*Transforma a matriz para a enviar*/
			/*Faz shift das bases a direita da que vai ser inserida*/
			for (i = N - 1; i >= j - 1; i++){
				bases[i + 1] = bases[i];
			}
			/*Insere nova base*/
			for (i = 0; i < N; i++){
				/*DUVIDA - Qual a base a ser inserida*/
				bases[j - 1][i] = v[i]
			}

			/*Chama LLL com a nova matriz*/
			lll(&bases, delta, h-1);
		}else{
			z++;
			/*chama LLL com a matriz actual*/
			lll(&bases, delta, h-1);
		}
	}
	return bases;
}