#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include<time.h>
#include "ResearchProjectCodeSubFunctionHeader.h"


//Creates the barrier for parallel function
pthread_barrier_t barrier;

void *CHEBYJACOBIThreadKernel(void *arguments){
    register int i,j,k;
    int ThreadVal;
    //Thread Value is the thread number
    ThreadVal = pthread_self();
    printf("In Kernel Thread %d\n",ThreadVal);
    //The arguments from void* so that it can be accessed
    struct arg_struct *args = (struct arg_struct *)arguments;
    int StartValue,EndValue,ExceptionStartValue,ExceptionEndValue;
    StartValue = args->ArrayThread[ThreadVal-1][0];
    EndValue = args->ArrayThread[ThreadVal-1][1];
    //printf("%d,%d\n",StartValue,EndValue);

    /*
    * Below if statements handle the case that the thread starts or ends at the maximum limits to avoid
    * overwriting boundary conditions
    */
    if(StartValue == 0){
        ExceptionStartValue = 1;
    }
    else{
        ExceptionStartValue = StartValue;
    }
    if(EndValue == (args->ROWS-1)){
        ExceptionEndValue = EndValue-1;
    }
    else{
        ExceptionEndValue = EndValue;
    }
    /*
    * Begin Initialisation
    */
    //Set all Array Values to 0
    for (i = StartValue; i <= EndValue; i++) {
        for (j = 0; j < args->COLS; j++) {
            args->ArrayNew[i][j] = 0.0;
            args->ArrayOld[i][j] = 0.0;
        }
    }
    pthread_barrier_wait(&barrier);
    double j0(double x);
    double BesselK;
    BesselK = 2.4048/args->RADIUS;
    for (i = StartValue; i <= EndValue; i++) {
        args->ArrayOld[args->ROWS-1][i] = j0(BesselK*(args->RStep*i));
        args->ArrayNew[args->ROWS-1][i] = j0(BesselK*(args->RStep*i));
    }


    pthread_barrier_wait(&barrier);
    //printf("Barrier 1 Reached, Thread = %d\n",ThreadVal);
    /*
    *End Initialisation
    */



    /*
    * Begin Relaxation
    */
    double Omega,Convergence,Factor1,Factor2,epsilon;
    Omega = 1;
    Convergence = args->Tolerance +1;
    while(Convergence > args->Tolerance){
        Convergence = 0;
        Omega = OmegaFunc(args->Iteration,Omega, args->Rho);
        args->ConvList[ThreadVal-1] = 0;
        pthread_barrier_wait(&barrier);
        //printf("Barrier 2 Reached, Thread = %d\n",ThreadVal);
        /* First Half Pass */
        for(k = 0;k <2;k++){
            for(i = ExceptionStartValue+k;i<=ExceptionEndValue;i=i+2){
                for(j = k;j<args->COLS-1;j = j+2){
                    if(j == 0){
                        Factor1 = 1;
                        epsilon = ((2*Factor1*args->ArrayOld[i][j+1] + args->ArrayOld[i+1][j] + args->ArrayOld[i-1][j])) - 4*args->ArrayOld[i][j];
                        args->ArrayNew[i][j] = args->ArrayOld[i][j] + (Omega*epsilon)/4;
                    }
                    //Relaxation
                    else{
                        Factor1 = 1 + 1/(double)(2*j);
                        Factor2 = 1 - 1/(double)(2*j);
                        epsilon = ((Factor1*args->ArrayOld[i][j+1] + Factor2*args->ArrayOld[i][j-1] + args->ArrayOld[i+1][j] + args->ArrayOld[i-1][j])) - 4*args->ArrayOld[i][j];
                        args->ArrayNew[i][j] = args->ArrayOld[i][j] + (Omega*epsilon)/4;
                    }

                    //Updates ConvList to allow for Convergence calculation at end
                    if(args->ArrayNew[i][j] != 0.0){
                        args->ConvList[ThreadVal-1] += fabs(args->ArrayOld[i][j]-args->ArrayNew[i][j])/fabs(args->ArrayNew[i][j]);
                    }
                    //Updates OldArray value with NewArray value
                    args->ArrayOld[i][j] = args->ArrayNew[i][j];
                }
            }
        }
        if(ThreadVal == 1){
            args->Iteration = args->Iteration + 0.5;
        }

        Omega = OmegaFunc(args->Iteration,Omega, args->Rho);
        pthread_barrier_wait(&barrier);
        /* Second Half Pass */
        for(k = 0;k <2;k++){
            for(i = ExceptionStartValue+k;i<=ExceptionEndValue;i=i+2){
                for(j = 1-k;j<args->COLS-1;j = j+2){
                    //Relaxation
                    if(j == 0){
                        Factor1 = 1 ;
                        epsilon = ((2*Factor1*args->ArrayOld[i][j+1] + args->ArrayOld[i+1][j] + args->ArrayOld[i-1][j])) - 4*args->ArrayOld[i][j];
                        args->ArrayNew[i][j] = args->ArrayOld[i][j] + (Omega*epsilon)/4;
                    }
                    //Relaxation
                    else{
                        Factor1 = 1 + 1/(double)(2*j);
                        Factor2 = 1 - 1/(double)(2*j);
                        epsilon = ((Factor1*args->ArrayOld[i][j+1] + Factor2*args->ArrayOld[i][j-1] + args->ArrayOld[i+1][j] + args->ArrayOld[i-1][j])) - 4*args->ArrayOld[i][j];
                        args->ArrayNew[i][j] = args->ArrayOld[i][j] + (Omega*epsilon)/4;
                    }

                    //Updates ConvList to allow for Convergence calculation at end
                    if(args->ArrayNew[i][j] != 0.0){
                        args->ConvList[ThreadVal-1] += fabs(args->ArrayOld[i][j]-args->ArrayNew[i][j])/fabs(args->ArrayNew[i][j]);
                    }
                    //Updates OldArray value with NewArray value
                    args->ArrayOld[i][j] = args->ArrayNew[i][j];
                }
            }
        }
        pthread_barrier_wait(&barrier);
        if(ThreadVal == args->ThreadNum){
            args->Iteration = args->Iteration + 0.5;
        }
        //printf("ConvList[%d] = %f\n",ThreadVal-1,args->ConvList[ThreadVal-1]);
        //printf("Barrier 5 Reached, Thread = %d\n",ThreadVal);
        for(i = 0;i<args->ThreadNum;i++){
            Convergence += args->ConvList[i];
        }
        Convergence = Convergence/(args->ROWS*args->COLS);
        pthread_barrier_wait(&barrier);
        //printf("Barrier 6 Reached, Thread = %d\n",ThreadVal);
        /*
        if(ThreadVal == 1){
            printf("Iteration = %f, Convergence = %.16f\n",args->Iteration,Convergence);
        }
        */
    }
    printf("Iteration = %f, Convergence = %.16f\n",args->Iteration,Convergence);
    return 0;
}




int CHEBYJACOBI(int ROWS,int COLS,double HEIGHT,double RADIUS,double TOLERANCE,int THREADNUM){
/*
* The Intent of this code is to merge the Chebyshev SOR with parallelisation via Pthreads
*/  register int i,j;
    double Rho,HStep,RStep;

    // Step Size in Z axis
    HStep = HEIGHT/(ROWS-1);
    //Step Size in R axis, R*2 as representing from -r to r due to being a cylinder
    RStep = (RADIUS)/(COLS-1);

    //Rho value for Jacobi method, used to calculate overrelaxation parameter
    Rho = RhoJacobi(ROWS,COLS,HEIGHT,RADIUS,HStep,RStep);

    //Memory Allocation for Old and New Array respectively
    double **ArrayNew = (double **) malloc(ROWS * sizeof(double*));
    for (i = 0; i < ROWS; i++) {
        ArrayNew[i] = (double *) malloc(COLS * sizeof(double));
    }
    double **ArrayOld = (double **) malloc(ROWS * sizeof(double*));
    for (i = 0; i < ROWS; i++) {
        ArrayOld[i] = (double *) malloc(COLS * sizeof(double));
    }

    /*
    * BEGIN THREAD INITIALISATION
    */

    //List of thread ids for each thread
    //56 because max number of cores for SPECTRE
    pthread_t thread[2000];
    pthread_attr_t attr;
    struct sched_param fifo_param;

    //Assigns a start and stop value for each thread,not yet
    int **ThreadInfo = (int**) malloc(THREADNUM*sizeof(int*));
    for(i = 0;i<THREADNUM;i++){
        ThreadInfo[i] = (int *) malloc(2 * sizeof(int));
        ThreadInfo[i][0] = 0;
        ThreadInfo[i][1] = 0;
    }
    //ThreadWorkDistributor gives start and end value for each thread
    ThreadWorkDistributor(ROWS,THREADNUM,ThreadInfo); //Tested, this bit works as intended

    //Barrier Initialisation
    pthread_barrier_init(&barrier, NULL, THREADNUM);
    //Create arguments struct for pthreads
    double *ConvList = (double*) malloc(THREADNUM*sizeof(double));

    struct arg_struct args;
    args.ArrayThread = ThreadInfo;
    args.ArrayNew = ArrayNew;
    args.ArrayOld = ArrayOld;
    args.ConvList = ConvList;
    args.ROWS = ROWS;
    args.COLS = COLS;
    args.HEIGHT = HEIGHT;
    args.RADIUS = RADIUS;
    args.RStep = RStep;
    args.HStep = HStep;
    args.Rho = Rho;
    args.Tolerance = TOLERANCE;
    args.ThreadNum = THREADNUM;
    printf("Args created\n");

    clock_t begin = clock();



    pthread_attr_init(&attr);
    pthread_attr_setinheritsched(&attr, PTHREAD_EXPLICIT_SCHED);
    pthread_attr_setschedpolicy(&attr, SCHED_FIFO);


    fifo_param.sched_priority = sched_get_priority_min(SCHED_FIFO);
    pthread_attr_setschedparam(&attr, &fifo_param);

    int policy;
    if(pthread_attr_getschedpolicy(&attr, &policy) != 0)
        fprintf(stderr, "Unable to get policy.\n");
    else{
        if(policy == SCHED_OTHER)
            printf("SCHED_OTHER\n");
        else if(policy == SCHED_RR)
            printf("SCHED_RR\n");
        else if(policy == SCHED_FIFO)
            printf("SCHED_FIFO\n");
    }


    for (i = 0; i < THREADNUM; i++) {
        //Creates each thread, args(thread id, ?,Function to run on, Value as pointer)
        pthread_create(&thread[i], &attr, CHEBYJACOBIThreadKernel, &args);
        }



    for (i = 0; i < THREADNUM; i++) {
        //Terminates all the threads
        pthread_join(thread[i], NULL);
    }
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Iterations Completed, In time %f Secs\n",time_spent);
    /*
    * END THREADS
    */


    /*
    * Writes to File
    */


    printf("Final Iteration = %.18d\n",(int)args.Iteration);
    printf("Started writing to file\n");
    FILE *FinalArr;
    FinalArr = fopen("ResearchProjectResultBessel100.csv", "w");
    for (i = 0; i < ROWS; i++) {
        for (j = 0; j < COLS; j++) {
            fprintf(FinalArr, "%.18f\n", ArrayNew[i][j]);
            //printf( "%.18f\n", ArrayOld[i][j]);
        }
    }
    fclose(FinalArr);
    printf("Finished writing to file\n");


    /*
    * Ends Write to File
    */
/*
    * Free Memory
    */
    for(i=0;i<ROWS;i++){
        free(ArrayNew[i]);
        free(ArrayOld[i]);
    }
    free(ArrayNew);
    free(ArrayOld);

    free(ConvList);

    for(i = 0;i<THREADNUM;i++){
        free(ThreadInfo[i]);
    }
    free(ThreadInfo);
    /*
    * Memory Freed
    */
return 0;
}




int main(){
    printf("Simulation Started\n");
    int ROWS, COLS,THREADNUM;
    double HEIGHT, RADIUS, TOLERANCE;
    ROWS = 100;
    COLS = ROWS;
    HEIGHT = 100;
    RADIUS = 100;
    TOLERANCE = 1*pow(10,-8);
    THREADNUM = 2;
    CHEBYJACOBI(ROWS,COLS,HEIGHT,RADIUS,TOLERANCE,THREADNUM);
    printf("Simulation Finished\n");
    return 0;
}




