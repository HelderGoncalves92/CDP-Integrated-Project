//Algoritmo em http://www.csie.nuk.edu.tw/~cychen/Lattices/Lattice%20Basis%20Reduction_%20Improved%20Practical%20Algorithms%20and%20Solving%20Subset%20Sum%20Problems.pdf - pagina 16

#include "ENUM.h"
#include "math.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

short nthreads;
short *u, **uT, **d, **delta, **v;
double **cT, **y, cL;
LEnum list = NULL;

pthread_mutex_t u_Mutex, toPop_Mutex;


void printVec(short bound, short id);

void initEnum(short n_threads){
    short i, auxDim = dim+1;
    nthreads = n_threads;
    
    //Allocate memory
    u = (short*)_mm_malloc(dim*sizeof(NLEnum), 64);
    list = (LEnum)_mm_malloc(sizeof(NLEnum), 64);
    
    list->count = 0;
    list->head = list->tail = NULL;;
    
    //Fill u with the shortest vector
    for(short i=1; i<=dim; i++)
        u[i] = 0;
        
    
    cL = B[0];
    u[0]=1;
    
    //Prepare memory to each thread
    delta = (short**)_mm_malloc(n_threads*sizeof(short*), 64);
    uT = (short**)_mm_malloc(n_threads*sizeof(short*), 64);
    d = (short**)_mm_malloc(n_threads*sizeof(short*), 64);
    v = (short**)_mm_malloc(n_threads*sizeof(short*), 64);
    cT = (double**)_mm_malloc(n_threads*sizeof(double*), 64);
    y = (double**)_mm_malloc(n_threads*sizeof(double*), 64);
    
    for(i=0; i<n_threads; i++){
        delta[i] = (short*)_mm_malloc(auxDim*sizeof(short), 64);
        uT[i] = (short*)_mm_malloc(auxDim*sizeof(short), 64);
        d[i] = (short*)_mm_malloc(auxDim*sizeof(short), 64);
        v[i] = (short*)_mm_malloc(auxDim*sizeof(short), 64);
        cT[i] = (double*)_mm_malloc(auxDim*sizeof(double), 64);
        y[i] = (double*)_mm_malloc(auxDim*sizeof(double), 64);
    }
    
    //Prepare to Pthreads
    pthread_mutex_init(&u_Mutex, NULL);
    pthread_mutex_init(&toPop_Mutex, NULL);
}


void moveDown(short id, short t, short s){
    
    double aux = 0.0;
    short i;
    
    for (i = t + 1; i <= s; i++){
        aux += uT[id][i] * mu[i][t];
    }
    
    y[id][t] = aux;
    uT[id][t] = v[id][t] = short(round(-aux));
    delta[id][t] = 0;
    
    if (uT[id][t] > -aux)
        d[id][t] = -1;
    else
        d[id][t] = 1;
    
    //Prepare cT[t] to next iteration
    cT[id][t] = cT[id][t + 1] + (y[id][t]*y[id][t] + 2*uT[id][t]*y[id][t] + uT[id][t]*uT[id][t]) * B[t];
}





void moveUP(short id, short t, short s){

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

int startSet(short id, Enum set){
	short i;
	short bound = set->bound+1;
    
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
            /*
            uT[id][bound-(set->level+1)] = set->sibling;
            delta[id][bound-(set->level+1)] = set->sibling;
            
            
            for(i=bound-2; i>= bound - set->level-2; i--)
                cT[id][i] = cT[id][i + 1] + (y[id][i]*y[id][i] + 2*uT[id][i]*y[id][i] + uT[id][i]*uT[id][i]) * B[i];
            */
            short t=bound-2;
            //printVec(dim, id);

            if(set->vec != NULL){
                
                short j=0;
                while(t>=(bound-set->level-1)){
                    moveDown(id, t, bound-1);
                    
                    if(cT[id][t] > cL){
                        printf("%d:Lower cL |Type:%d, Sib:%d, Level:%d\n",bound, set->type,set->sibling, set->level);
                        return 1;
                    }
                        
                    while(uT[id][t] != set->vec[j]){
                        
                        moveUP(id, t, bound-1);
                        if(abs(uT[id][t]) > abs(set->vec[j])){
                            printf("%d:ERRO (Jump Sib) Type:%d, Sib:%d, Level:%d\n",bound, set->type,set->sibling, set->level);
                            return 1;
                        }
                    }
                    
                    
                    t--;
                    j++;
                }
            
            }else{
            
                while(t>=(bound-set->level-1)){
                    moveDown(id, t, bound-1);
                    
                    if(cT[id][t] > cL){
                        printf("%d:Lower cL |Type:%d, Sib:%d, Level:%d\n",bound, set->type,set->sibling, set->level);
                        return 1;
                    }
                    
                    if(t == (bound-set->level-1)){
                        
                        while(uT[id][t] != set->sibling){
                            
                            moveUP(id, t, bound-1);
                            if(abs(uT[id][t]) > abs(set->sibling)){
                                printf("%d:ERRO (Jump Sib) Type:%d, Sib:%d, Level:%d\n",bound, set->type,set->sibling, set->level);
                                return 1;
                            }
                        }
                        
                    } else if(uT[id][t] != 0){
                        printf("%d:ERRO (Not '0') Type:%d, Sib:%d, Level:%d\n",bound, set->type,set->sibling, set->level);
                        printVec(dim, id);

                        return 1;
                    }
                    
                    t--;
                }
            }
        }
    return 0;
}


Enum newEnumElem(short bound, short sibling, short type, short level, short *vec){
    Enum st = (Enum)_mm_malloc(sizeof(NEnum),64);
    st->next = NULL;
    st->bound = bound-1;
    st->type = type;
    st->sibling = sibling;
    st->level = level;
    st->vec = vec;

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


//ENUM accordingly C. P. Schnorr && M. Euchner
void EnumSET(Enum set, short id){
    
    double aux;
	short s, t;
	short i;
	short bound = set->bound;
    short toCopy = bound + 1;
    
    if(startSet(id, set))
        return;
  //  printVec(dim, id);

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
        bound = bound - set->level-1;
        moveDown(id, bound, s);
        t=bound;

    }

    while(t <= bound){
      // printVec(dim, id);
        if (cT[id][t] < cL){
            if (t > 0){
                //moveDown
                t--;
                aux = 0.0;
                
                for (i = t + 1; i <= s; i++){
                    aux += uT[id][i] * mu[i][t];
                }
                
                y[id][t] = aux;
                uT[id][t] = v[id][t] = short(round(-aux));
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
                printf("%d:UPDATE Type:%d, Sib:%d, Level:%d -->%f / %f\n",s+1, set->type,set->sibling, set->level ,cT[id][0],cL);
                //printVec(dim, id);
                
                //Lock critical zone
                pthread_mutex_lock(&u_Mutex);
                    cL = cT[id][0];
                    memcpy(&u[0], uT[id], toCopy*sizeof(short));
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
 //   printf("Fim:%d || %d\n",toCopy, set->sibling);
}

void* threadHander(void* vID){
    
    short id = *((short *) vID);
    Enum set = NULL;

    short* vec = (short*)_mm_malloc(4*sizeof(short), 64);
    vec[0]=-1;
    vec[1]= 1;
    vec[2]= 0;

    set = newEnumElem(dim-id, 1, 3, 3, vec);
    
    //if(id==0){
     //EnumSET(set, id);
        //printVec(dim, id);
    //}

    while (list->count>0) {
        set = pop();
        if(set){
          //  printf("%d - %d\n",id, set->bound);
            EnumSET(set, id);
        }
    }

    return NULL;
}


//From vector '0 0 0 0 0'
void creatTasks(short bound, short level){
	short i, j;
	//short depth = bound - level;
    Enum set = NULL;

    short* vec = (short*)_mm_malloc(4*sizeof(short), 64);
    vec[0]=-1;
    vec[1]=1;
    vec[2]=0;

    	
    for(i=1; i<=4; i++){
	set = newEnumElem(bound, 1, 3, 3, vec);
	addTail(set);
        set = newEnumElem(bound, -1, 3, 3, vec);
        addTail(set);
        set = newEnumElem(bound, 2, 2, 3, vec);
        addTail(set);
    }
    set = newEnumElem(bound, 0, 2, 3, vec);
    addTail(set);
    

    for(i=1; i<=level; i++){

        set = newEnumElem(bound, 1, 3, i, NULL);
        addTail(set);
        set = newEnumElem(bound, -1, 3, i, NULL);
        addTail(set);
        set = newEnumElem(bound, 2, 2, i, NULL);
        addTail(set);
    }
    set = newEnumElem(bound, 0, 2, i, NULL);
    addTail(set);
}


short* ENUM(){
    
	short i, n = 1, creatRange = dim - nthreads, MAX_DEPTH = 0.4*dim, divRange = 0.6*dim;
    Enum set = NULL;
	for (i = dim; i > creatRange; i--){
		creatTasks(i, 3);
	}
    
    for(i=creatRange; i>divRange; i--){
        set = newEnumElem(i, 0, 3, 1, NULL);
        addTail(set);
		n++;
        set = newEnumElem(i, 1, 3, 1, NULL);
        addTail(set);
        n++;
        set = newEnumElem(i, -1, 3, 1, NULL);
        addTail(set);
        n++;
        set = newEnumElem(i, 2, 2, 1, NULL);
        addTail(set);
        n++;
    }
    
    //Prepare list with tasks
    for(; i>MAX_DEPTH; i--){
        set = newEnumElem(i, 1, 1, 1, NULL);
        addTail(set);
        n++;
    }
    set = newEnumElem(i, 1, 0, 1, NULL);
    addTail(set);
    
    
    pthread_t tHandles[nthreads];
    
    //Threads Start
    for (i = 0; i < nthreads; i++) {
        short *threadNum = (short*)malloc(sizeof (short));
        *threadNum = i;
        pthread_create(&tHandles[i], NULL, threadHander, (void *)threadNum);

    }
    
    //Threads Join
    for (i = 0; i < nthreads; i++)
        pthread_join(tHandles[i], NULL);
    
    return u;
}


void printVec(short bound, short id){
    short i;
    
    short dimension = bound;
    
    printf("uT:");
    for(i=0; i<dimension; i++)
        printf("%d ",uT[id][i]);
    printf("\n");
  /*
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
    */
   /* printf("cT:");
    for(i=0; i<dimension; i++)
        printf("%.2f ",cT[id][i]);
    printf("\n");
    printf("Y:");
    for(i=0; i<dimension; i++)
        printf("%.2f ",y[id][i]);
    printf("\n");
    printf("\n");
   */
}
