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

void startSet(int id, Enum set){
    int i, bound = set->bound+1;
    
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
    if(set->type == 0){
        uT[id][0] = 1;
        for(i=1; i<=bound; i++)
            cT[id][i] = 0.0;
        
    } else if(set->type == 1){
        //Base Start
        cT[id][bound] = 0.0;
        cT[id][bound-1] = B[bound-1];
        uT[id][bound-1] = 1;
        
    } else
        if(set->type == 2 || set->type == 3){
            cT[id][bound] = 0.0;
            cT[id][bound-1] = B[bound-1];
            uT[id][bound-1] = 1;
            
            uT[id][bound-2] = set->sibling;
            delta[id][bound-2] = set->sibling;
            
            
            for(i=bound-2; i>= bound - set->level-2; i--)
                cT[id][i] = cT[id][i + 1] + (y[id][i]*y[id][i] + 2*uT[id][i]*y[id][i] + uT[id][i]*uT[id][i]) * B[i];
           // set->bound--;
            
        }
}


Enum newEnumElem(int bound, int sibling, int type, int level){
    Enum st = (Enum)_mm_malloc(sizeof(NEnum),64);
    st->next = NULL;
    st->bound = bound-1;
    st->type = type;
    st->sibling = sibling;
    st->level = level;

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
    
    if(!list->head){
        pthread_mutex_unlock(&toPop_Mutex);
        return NULL;
    }
        
    aux = list->head;
    list->head = list->head->next;
    list->count--;
    
    //Unlock critical zone
    pthread_mutex_unlock(&toPop_Mutex);
 
    return aux;
}

void moveDown(int id, int bound, int s, int times){
    
    int i, t = bound;
    
    double aux;
  //  printVec(dim, id);
    
    while(times>0){
        times--;
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
        cT[id][t] = cT[id][t + 1] + (y[id][t]*y[id][t] + 2*uT[id][t]*y[id][t] + uT[id][t]*uT[id][t]) * B[t];
    //    printVec(dim, id);
    }

}


void moveDownUP(int id, int bound, int s, int times){
    
    int i, t = bound;
    
    double aux;
    printVec(dim, id);
    
    while(times>0){
        times--;
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
        cT[id][t] = cT[id][t + 1] + (y[id][t]*y[id][t] + 2*uT[id][t]*y[id][t] + uT[id][t]*uT[id][t]) * B[t];
        printVec(dim, id);
        
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
        cT[id][t] = cT[id][t + 1] + (y[id][t]*y[id][t] + 2*uT[id][t]*y[id][t] + uT[id][t]*uT[id][t]) * B[t];
        printVec(dim, id);
    }
    
}



//ENUM accordingly C. P. Schnorr && M. Euchner
void EnumSET(Enum set, int id){
    
    double aux;
    int s, t, i, bound = set->bound;
    int toCopy = bound + 1;
    
    startSet(id, set);

    //Start on leaf (like Schnorr)
    if(set->type==0){
        s = t = 0;
        cT[id][t] = cT[id][t + 1] + (y[id][t]*y[id][t] + 2*uT[id][t]*y[id][t] + uT[id][t]*uT[id][t]) * B[t];

    //One Gamma
    } else if(set->type == 1){
        s = t = bound;

    //Starts on given sibling and go on
    } else if(set->type == 2){
        s = bound;
        bound = bound - set->level;
        t=bound;
        
    //Search just in one sibling
    } else {
        s = bound;
        bound = bound - set->level;
        moveDown(id, bound, s, 1);
        bound--;
        t=bound;

    }

    
   /*
    moveDown(id, bound, s, 1);
    
    bound--;
    t = bound;

    printVec(dim, id);
    moveDownUP(id, bound, s, 1);
    
    printVec(dim, id);
    */
/*

    cT[id][t] = cT[id][t + 1] + (y[id][t]*y[id][t] + 2*uT[id][t]*y[id][t] + uT[id][t]*uT[id][t]) * B[t];
    printVec(dim, id);
  */
    while(t <= bound){
       //printVec(dim, id);
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
                cT[id][t] = cT[id][t + 1] + (y[id][t]*y[id][t] + 2*uT[id][t]*y[id][t] + uT[id][t]*uT[id][t]) * B[t];
                
            }
            else{
                //UpdateVector
                printf("%d:UPDATE-%d\n",s, set->type);
                //printVec(dim, id);
                
                //Lock critical zone
                pthread_mutex_lock(&u_Mutex);
                    cL = cT[id][0];
                    memcpy(&u[0], uT[id], toCopy*sizeof(int));
                    for(i=toCopy; i<dim; i++)
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
            cT[id][t] = cT[id][t + 1] + (y[id][t]*y[id][t] + 2*uT[id][t]*y[id][t] + uT[id][t]*uT[id][t]) * B[t];
            
        }
    }
    printf("Fim:%d\n",bound);
}

void* threadHander(void* vID){
    
    int id = *((int *) vID);
    Enum set = NULL;
    
    while (list->count>0) {
        set = pop();
        if(set){
            printf("%d - %d\n",id, set->bound);
            EnumSET(set, id);
        }
    }

    return NULL;
}


int* ENUM(){
    
    int i, n=1, MAX_DEPTH = dim*0.8;
    Enum set = NULL;
    
    set = newEnumElem(55, 0, 3, 1);
    addTail(set);
    n++;
    set = newEnumElem(55, 1, 3, 1);
    addTail(set);
    n++;
    set = newEnumElem(55, -1, 3, 1);
    addTail(set);
    n++;
    set = newEnumElem(55, 2, 2, 1);
    addTail(set);
    n++;
    
    set = newEnumElem(54, 0, 3, 1);
    addTail(set);
    n++;
    set = newEnumElem(54, 1, 3, 1);
    addTail(set);
    n++;
    set = newEnumElem(54, -1, 3, 1);
    addTail(set);
    n++;
    set = newEnumElem(54, 2, 2, 1);
    addTail(set);
    n++;
    
    
    
    
    //Prepare list with tasks
    for(i=dim-2; i>MAX_DEPTH; i--){
        set = newEnumElem(i, 1, 1, 1);
        addTail(set);
        n++;
    }
    set = newEnumElem(i, 1, 0, 1);
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
    
  /*  printf("v:");
    for(i=0; i<dimension; i++)
        printf("%d ",v[id][i]);
    printf("\n");
    */
    printf("cT:");
    for(i=0; i<dimension; i++)
        printf("%.2f ",cT[id][i]);
    printf("\n");
  /*  printf("Y:");
    for(i=0; i<dimension; i++)
        printf("%.2f ",y[id][i]);
    printf("\n");
    */
    printf("\n");
}
