//Algoritmo em http://www.csie.nuk.edu.tw/~cychen/Lattices/Lattice%20Basis%20Reduction_%20Improved%20Practical%20Algorithms%20and%20Solving%20Subset%20Sum%20Problems.pdf - pagina 16

#include "ENUM.h"
#include "math.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

int nthreads;
int *u, **uT, **d, **delta, **v;
double **cT, **y, cL;
LEnum list = NULL;

pthread_mutex_t u_Mutex, toPop_Mutex;


void printVec(int bound, int id);

void initEnum(int n_threads){
    int i, auxDim = dim+1;
    nthreads = n_threads;
    
    //Allocate memory
    u = (int*)_mm_malloc(dim*sizeof(NLEnum), 64);
    list = (LEnum)_mm_malloc(sizeof(NLEnum), 64);
    
    list->count = 0;
    list->head = list->tail = NULL;;
    
    //Fill u with the shortest vector
    for(int i=1; i<=dim; i++)
        u[i] = 0;
    
    cL = B[0];
    u[0]=1;
    
    //Prepare memory to each thread
    delta = (int**)_mm_malloc(n_threads*sizeof(int*), 64);
    uT = (int**)_mm_malloc(n_threads*sizeof(int*), 64);
    d = (int**)_mm_malloc(n_threads*sizeof(int*), 64);
    v = (int**)_mm_malloc(n_threads*sizeof(int*), 64);
    cT = (double**)_mm_malloc(n_threads*sizeof(double*), 64);
    y = (double**)_mm_malloc(n_threads*sizeof(double*), 64);
    
    for(i=0; i<n_threads; i++){
        delta[i] = (int*)_mm_malloc(auxDim*sizeof(int), 64);
        uT[i] = (int*)_mm_malloc(auxDim*sizeof(int), 64);
        d[i] = (int*)_mm_malloc(auxDim*sizeof(int), 64);
        v[i] = (int*)_mm_malloc(auxDim*sizeof(int), 64);
        cT[i] = (double*)_mm_malloc(auxDim*sizeof(double), 64);
        y[i] = (double*)_mm_malloc(auxDim*sizeof(double), 64);
    }
    
    //Prepare to Pthreads
    pthread_mutex_init(&u_Mutex, NULL);
    pthread_mutex_init(&toPop_Mutex, NULL);
}

void startSet(int id, int bound, int sibling, int type){
    int i;
    
    //Clean vectors from preciously executions
    for(i=0; i<=bound; i++)
        y[id][i] = 0.0;
    
    for(i=0; i<=bound; i++)
        delta[id][i] = 0;
    for(i=0; i<=bound; i++)
        d[id][i] = 1;
    for(i=0; i<=bound; i++)
        v[id][i] = 0;
    for(i=0; i<=bound; i++)
        uT[id][i] = 0;
    
    
    //Prepare to start by type
    if(type == 1){
        cT[id][bound] = 0.0;
        cT[id][bound-1] = B[bound-1];
        uT[id][bound-1] = 1;
        uT[id][0] = 1;
        
        double aux;
        for(i=bound-1; i>=0; i--){
            aux = y[id][i] + uT[id][i];
            cT[id][i] = cT[id][i + 1] + (aux * aux) * B[i];
        
    }else{
        uT[id][0] = 1;
        for(i=0; i<=bound; i++)
            cT[id][i] = 0.0;
    }
}


Enum newEnumElem(int bound, int sibling, int type){
    Enum st = (Enum)_mm_malloc(sizeof(NEnum),64);
    st->next = NULL;
    st->bound = bound-1;
    st->type = type;
    st->sibling = sibling;

    return st;
}

void addTail(Enum newSet){
    if(list->count != 0){
        list->tail->next=newSet;
        list->tail = newSet;
        list->count++;
    }else{
        list->head = list->tail = newSet;
        list->count=1;
    }
}

Enum pop(){
    Enum aux;
    
    //Lock critical zone
    pthread_mutex_lock(&toPop_Mutex);
    
    aux = list->head;
    list->head = list->head->next;
    list->count--;
    
    //Unlock critical zone
    pthread_mutex_unlock(&toPop_Mutex);
 
    return aux;
}



//ENUM accordingly C. P. Schnorr && M. Euchner
void EnumSET(Enum set, int id){
    
    double aux;
    int s, t, i;
    int bound = set->bound;
    
    startSet(id, bound+1, set->sibling, set->type);
    s = t = 0;

    bound++;
    aux = y[id][t] + uT[id][t];
    cT[id][t] = cT[id][t + 1] + (aux * aux) * B[t];
    
    while(t < bound){
        
        if (cT[id][t] < cL){
            if (t > 0){
                //moveDown
                t--;
                aux = 0;
                
                for (i = t + 1; i <= s; i++){
                    aux += uT[id][i] * mu[i][t];
                }
                
                y[id][t] = aux;
                uT[id][t] = v[id][t] = int(round(-aux));
                delta[id][t] = 0;
                
                if (uT[id][t] > -aux){
                    d[id][t] = -1;
                }
                else{
                    d[id][t] = 1;
                }
                
                //Prepare cT[t] to next iteration
                aux = y[id][t] + uT[id][t];
                cT[id][t] = cT[id][t + 1] + (aux * aux) * B[t];
                
            }
            else{
                //UpdateVector
                printf("%d:UPDATE\n",s);
                //printVec(set);
                
                //Lock critical zone
                pthread_mutex_lock(&u_Mutex);
                    cL = cT[id][0];
                    memcpy(&u[0], uT[id], bound*sizeof(int));
                    for(i=bound; i<dim; i++)
                        u[i] = 0;
                
                //Unlock critical zone
                pthread_mutex_unlock(&u_Mutex);
            }
        }
        else{
            //moveUp
            t++;
            s = (s<t)?t:s; //Get max value
            
            if(t < s){
                delta[id][t] = -delta[id][t];
            }
            if(delta[id][t]*d[id][t] >= 0){
                delta[id][t] += d[id][t];
            }
            uT[id][t] = v[id][t] + delta[id][t];
            
            //Prepare cT[t] to next iteration
            aux = y[id][t] + uT[id][t];
            cT[id][t] = cT[id][t + 1] + (aux * aux) * B[t];
            
        }
    }
    printf("Fim:%d\n",bound);
}

void* threadHander(void* vID){
    
    int id = *((int *) vID);
    Enum set = NULL;
    
    while (list->count>0) {
        set = pop();
        printf("%d - %d\n",id, set->bound);
        EnumSET(set, id);
    }

    return NULL;
}


int* ENUM(){
    
    int i, n=1, MAX_DEPTH = dim*0.7;
    Enum set = NULL;
    
    //Prepare list with tasks
    for(i=dim; i>MAX_DEPTH; i--){
        set = newEnumElem(i, 0, 1);
        addTail(set);
        n++;
    }
    set = newEnumElem(i, 0, 0);
    addTail(set);
    
    
    pthread_t tHandles[nthreads];
    
    //Threads Start
    for (i = 0; i < nthreads; i++) {
        int *threadNum = (int*)malloc(sizeof (int));
        *threadNum = i;
        pthread_create(&tHandles[i], NULL, threadHander, (void *)threadNum);

    }
    
    
    //Threads Join
    for (i = 0; i < nthreads; i++)
        pthread_join(tHandles[i], NULL);
    
    return u;
}


void printVec(int bound, int id){
    int i;
    
    int dimension = bound+1;
    
    printf("uT:");
    for(i=0; i<dimension; i++)
        printf("%d ",uT[id][i]);
    printf("\n");
    
    printf("d:");
    for(i=0; i<dimension; i++)
        printf("%d ",d[id][i]);
    printf("\n");
    
    printf("delta:");
    for(i=0; i<dimension; i++)
        printf("%d ",delta[id][i]);
    printf("\n");
    
    printf("v:");
    for(i=0; i<dimension; i++)
        printf("%d ",v[id][i]);
    printf("\n");
    
    printf("cT:");
    for(i=0; i<dimension; i++)
        printf("%.2f ",cT[id][i]);
    printf("\n");
    printf("Y:");
    for(i=0; i<dimension; i++)
        printf("%.2f ",y[id][i]);
    printf("\n");
    
    printf("\n");
}
