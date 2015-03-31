//N = Dimensao, por enquanto definida globalmente para simplificar
#define N  2;

#include "BKZ.h"

bool passvec(double v[]){
	int i;
	if (v[0]!=1){return false;}
	for(i = 1; i < N; i++){
		if (v[i] != 0) return false;
	}
	return true;
}

double** BKZ(double *bases[], double *u[], double *c[], int beta, double delta){
	int k, h, i, l;
	double v[N], vaux[N];
	double* aux;

	int z = 0, j = 1;

	lll(&bases, delta);
	while (z < N - 1){
		j = (j * mod(N-1)) + 1; 
		k = min(j + beta - 1, N);
		h = min(k + 1, N);
		////cria nova matriz de ortogonalizacao para enviar para o ENUM
		//for (i = 0; i < N; i++){
		//	vaux[i] = vectorNorm(&bases[i],2);
		//}
		v = ENUM();
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
			lll(&bases, delta);
		}else{
			z++;
			/*chama LLL com a matriz actual*/
			lll(&bases, delta);
		}
	}
	return bases;
}