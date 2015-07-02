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
LEnum list_Urgent = NULL;

pthread_mutex_t u_Mutex, toHead_Mutex;


void printVec(short bound, short id);

void initEnum(short n_threads){
    short i, auxDim = dim+1;
    nthreads = n_threads;
    
    //Allocate memory
    u = (short*)_mm_malloc(dim*sizeof(NLEnum), 64);
    list = (LEnum)_mm_malloc(sizeof(NLEnum), 64);
    list_Urgent = (LEnum)_mm_malloc(sizeof(NLEnum), 64);
    
    list->count = 0;
    list->head = list->tail = NULL;
    list_Urgent->count = 0;
    list_Urgent->head = list_Urgent->tail = NULL;;

    
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
    pthread_mutex_init(&toHead_Mutex, NULL);
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
        delta[id][bound-1]=1;
        
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
                printf("My bound: %d, sibling: %d\n",t, set->sibling);
                set->vec[3] = set->sibling;
                printf("VECTOR: 0: %d, 1: %d, 2: %d, 3: %d\n",set->vec[0], set->vec[1], set->vec[2], set->vec[3]);
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
        
        } else if (set->type == 4){

            cT[id][set->level] = set->cT;
            cT[id][set->level+1] = set->cT2;
            memcpy(&uT[id][set->level], set->vec, (bound-set->level+1)*sizeof(short));
           /* printf("uT:");
            for(int i=0; i<(bound- set->level); i++)
                printf("%d ",set->vec[i]);
            printf("\n");
            */free(set->vec);
            
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
    if(vec!=NULL){
        st->vec = (short*)_mm_malloc(4*sizeof(short),64);
        memcpy(&st->vec[0], vec, 4*sizeof(short));
    }else{
        st->vec = NULL;
    }
    return st;
}

Enum newEnumElemV2(short bound, short level, short type, short *vec, double cT, double cT2){
    Enum st = (Enum)_mm_malloc(sizeof(NEnum),64);
    st->next = NULL;
    
    st->bound = bound;
    st->type = type;
    st->level = level;
    st->cT=cT;
    st->cT2=cT2;
    
    st->vec = (short*)_mm_malloc((bound-level+1)*sizeof(short),64);
    memcpy(st->vec, vec, (bound-level+1)*sizeof(short));
    
   /* printf("uT:");
    for(int i=0; i<(bound-level+1); i++)
        printf("%d ",vec[i]);
    printf("\n");
    */
    return st;
}

void addTail(LEnum l, Enum newSet){
    //Lock critical zone
    pthread_mutex_lock(&toHead_Mutex);
    
    if(l->count != 0){
        l->tail->next=newSet;
        l->tail = newSet;
        l->count++;
    }else{
        l->head = l->tail = newSet;
        l->count=1;
    }
    
    //Lock critical zone
    pthread_mutex_unlock(&toHead_Mutex);
}

void addHead(LEnum l, Enum newSet){
    
    //Lock critical zone
    pthread_mutex_lock(&toHead_Mutex);
    
    if(l->count != 0){
        newSet->next=l->head;
        l->head = newSet;
        l->count++;
        
    }else{
        l->head = l->tail = newSet;
        l->count=1;
    }
    
    //Unlock critical zone
    pthread_mutex_unlock(&toHead_Mutex);
}

Enum pop(LEnum l){
    Enum aux;
    
    //Lock critical zone
    pthread_mutex_lock(&toHead_Mutex);
    
    if(!l->head){
        pthread_mutex_unlock(&toHead_Mutex);
        return NULL;
    }
        
    aux = l->head;
    l->head = l->head->next;
    l->count--;
    
    //Unlock critical zone
    pthread_mutex_unlock(&toHead_Mutex);
 
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

    //printVec(dim,id);
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
    } else if(set->type == 3){
        s = bound;
        bound = bound - set->level-1;
        moveDown(id, bound, s);
        t=bound;
        
    } else if(set->type == 4){
       
        s = set->bound;
        t = bound = set->level;
        delta[id][t] = uT[id][t];
        //cT[id][t] = cT[id][t + 1] + (y[id][t]*y[id][t] + 2*uT[id][t]*y[id][t] + uT[id][t]*uT[id][t]) * B[t];
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
   // printf("%dEND Type:%d, Sib:%d, Level:%d\n",s+1, set->type,set->sibling, set->level);
}

void EnumCreatTasks(Enum set, short id, short depth){
    
    double aux;
    short i, s, t, bound;
    
    //Always in type 1 -- One Gamma
    s = t = bound = set->bound;
    startSet(id, set);

    while(t <= bound){
        
        if (cT[id][t] < cL){
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
            
            if (t == depth-1){
                

                if(uT[id][bound]==2) return;

                Enum st = newEnumElemV2(bound, t, 4, &uT[id][t], cT[id][t], cT[id][t+1]);

                if(uT[id][bound-1]==0  ||  s>dim-3){
                    addHead(list_Urgent, st);}
                else addHead(list, st);
                
/*                if(s>dim-3){
                    if(uT[id][bound-1]==0)
                        addHead(list_Urgent, st);
                    else
                        addTail(list_Urgent, st);
                    
                } else if(uT[id][bound-1]==0)
                    addTail(list_Urgent, st);
                else
                    addHead(list, st);
  */
                //moveUp
                t++;
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
}

void* threadHander(void* vID){
    
    short id = *((short *) vID);
    Enum set = NULL;

    while (list->count>0) {
        if(!(set = pop(list_Urgent)))
            set = pop(list);
        
        if(set)
            EnumSET(set, id);
    }

    return NULL;
}


void creatTasks(short bound, short level){
	short i, j;
	//short depth = bound - level;
    Enum set = NULL;
    bool finished = false;
    short* vec = (short*)_mm_malloc(4*sizeof(short), 64);
    
    
//    for(i=1; i<=level; i++){
//
//        set = newEnumElem(bound, 1, 3, i, NULL);
//        addTail(set);
//        set = newEnumElem(bound, -1, 3, i, NULL);
//        addTail(set);
//        set = newEnumElem(bound, 2, 2, i, NULL);
//        addTail(set);
//    }
//    set = newEnumElem(bound, 0, 2, i, NULL);
//    addTail(set);
}


short* ENUM(){
    
    short i, n = 1, creatRange = dim - 5, MAX_DEPTH = 0.7*dim, divRange = 0.8*creatRange;
    Enum set = NULL;
    
    
    for(i=creatRange; i>divRange; i--){
        set = newEnumElem(i, 0, 3, 1, NULL);
        addTail(list_Urgent, set);
        n++;
        set = newEnumElem(i, 1, 3, 1, NULL);
        addTail(list, set);
        n++;
        set = newEnumElem(i, -1, 3, 1, NULL);
        addTail(list, set);
        n++;
        set = newEnumElem(i, 2, 2, 1, NULL);
        addTail(list, set);
        n++;
    }
    
    for(; i>MAX_DEPTH; i--){
        set = newEnumElem(i, 1, 1, 1, NULL);
        addTail(list, set);
        n++;
    }
    
    //Prepare list with tasks
    set = newEnumElem(MAX_DEPTH, 1, 0, 1, NULL);
    addTail(list, set);

    
    //Threads Start
    pthread_t tHandles[nthreads];

    for (i=1; i<nthreads ; i++) {
        short *threadNum = (short*)malloc(sizeof (short));
        *threadNum = i;
        pthread_create(&tHandles[i], NULL, threadHander, (void *)threadNum);
    }
    
    //more tasks
    for (i = dim; i>creatRange ; i--){
        set = newEnumElem(i, 1, 1, 1, NULL);
        EnumCreatTasks(set, 0, i-5);
    }
    
    short *threadNum = (short*)malloc(sizeof (short));
    *threadNum = 0;
    pthread_create(&tHandles[0], NULL, threadHander, (void *)threadNum);
    
    
    
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
   
    printf("cT:");
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
