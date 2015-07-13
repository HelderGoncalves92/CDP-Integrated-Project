
#include "BKZ.h"

void initBKZ(int dimension){
    
    dim=dimension;
    
    initENUM();
    initStructsLLL(dimension);
}


bool passvec(int* v, int index){
	int i;

	if (v[index] != 1){ return false; }
	for(i = 0; i < dim; i++){
        if (i != index && v[i] != 0){
            return false;
        }
	}
	return true;
}

void BKZ(long** bases, int beta, double delta){
	int k, h, i, l;
    int *v;
	int z = 0, j = 0;

	lll(bases, delta, dim);
	while (z < dim - 1){
		j = (j % (dim-1)) + 1;
		k = fmin(j + beta - 1, dim);
		h = fmin(k + 1, dim);

		v = ENUM(j - 1, k - 1);
		if (!passvec(v, j-1))
		{
			z = 0;
			//Transforma a matriz para a enviar
			//Insere nova base
            for(i=0;i<dim; i++)
                 bases[dim][i] = 0;
            
			for (i = 0; i < dim; i++){
				//alterar vetor
				for (l = 0; l < dim; l++){
					bases[dim][l] += v[i] * bases[i][l];
				}
			}
			//Faz shift da base
			shiftVector(bases, j - 1, dim);

			//Chama LLL com a nova matriz, acaba em h pois vai ter mais 1 vetor a computar pela LLL
			lll(bases, delta, h+1);
		}else{
			z++;
			//chama LLL com a matriz actual
			lll(bases, delta, h);
		}
        
    }
}