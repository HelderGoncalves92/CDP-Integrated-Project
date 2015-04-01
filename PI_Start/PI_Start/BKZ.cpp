
#include "BKZ.h"

bool passvec(int* v, int index){
	int i;
	if (v[index] != 1){ return false; }
	for(i = 0; i < dimVector; i++){
		if (i != index && v[i] != 0) return false;
	}
	return true;
}

double** BKZ(double** bases, double *mu[], double *c, int beta, double delta){
	int k, h, i;
    int *v;

	int z = 0, j = 0;

	lll(bases, delta, dimVector-1);
	while (z < dimVector - 1){
		j = (j%(dimVector-1)) + 1;
		k = fmin(j + beta - 1, dimVector);
		h = fmin(k + 1, dimVector);
		
		v = EnumWPrun(mu, c, j - 1, k - 1);
		if (!passvec(v, j-1))
		{
			z = 0;
			//Transforma a matriz para a enviar
			//Insere nova base
			for (i = 0; i < dimVector; i++){
				//alterar vetor
                bases[dimVector][i] = v[i] ;
			}
			//Faz shift da base
			shiftVector(bases, j - 1, dimVector);


			//Chama LLL com a nova matriz
			lll(bases, delta, h-1);
		}else{
			z++;
			//chama LLL com a matriz actual
			lll(bases, delta, h-1);
		}
	}
	return bases;
}