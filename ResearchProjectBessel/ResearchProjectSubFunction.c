#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double RhoJacobi(int ROWS,int COLS,double HEIGHT,double RADIUS,double HStep,double RStep){
    /*
    * Value needed for the optimal omega value as given in ~ pg865 of Numerical Reciples
    * Is used in the Omega Func but needs to only be calculated once during the simulation run.
    */
    double RhoJacobi;
    RhoJacobi = (cos(M_PI/ROWS) + pow(HStep/RStep,2.)*cos(M_PI/COLS))/(1+pow(HStep/RStep,2.));
    printf("RhoJacobi = %f\n",RhoJacobi);
    RhoJacobi = RhoJacobi;
    return RhoJacobi;

}

double OmegaFunc(double Iteration,double Omega,double Rho){
    /*
    * This code recursively calculates the new Omega for the Chebyshev overrelaxation
    * As in omega changes with iterations and so this function is called every half step
    * due to red-black ordering
    */
    double om;
    if(Iteration == 0){
        om = 1;
    }
    else if(Iteration == 0.5){
        om = 1/(1-(Rho*Rho)/2);
    }
    else{
        om = 1/(1.0-(Rho*Rho*Omega)/4);
    }

    // Implements a check that Omega is between 0 and 2 since otherwise result will not converge
    // Can be removed for optimisation purposes when simulation works
    if(Omega > 2 || Omega < 0){
        printf("Error Occured in OmegaFunc \n Omega is not in the valid range [0,2]");
        exit(1);
    }
    return om;
}

void ThreadWorkDistributor(int ROWS, int THREADNUM,int **ThreadInfo){
    int Remnant;
    int k;
    int RowsPerThread;
    if(THREADNUM > ROWS){
        printf("Error: ThreadWorkDistributor failed due to more threads than rows\n");
        printf("Soln: Reduce the number of threads to less than the number of rows\n");
        exit(1);
    }
    Remnant = ROWS%THREADNUM;
    RowsPerThread = ROWS/THREADNUM;
    ThreadInfo[0][0] = 0;
    ThreadInfo[0][1] = (RowsPerThread)-1;
    if(Remnant>0){
        ThreadInfo[0][1] += 1;
        Remnant--;
    }
    for(k = 1;k<THREADNUM;k++){
        ThreadInfo[k][0] = ThreadInfo[k-1][1]+1;
        ThreadInfo[k][1] = ThreadInfo[k-1][1]+RowsPerThread;
        if(Remnant>0){
            ThreadInfo[k][1] +=1;
            Remnant--;
        }
    }

    for(k = 0;k<THREADNUM-1;k++){
        if(ThreadInfo[k][1] != ThreadInfo[k+1][0] - 1){
            printf("Error: ThreadWorkDistributor Failed due to non-contiguous work distribution\n");
            printf("Soln: Reduce number of threads");
            exit(2);
        }
    }
    if(ThreadInfo[THREADNUM-1][1] != ROWS-1){
            printf("Error: ThreadWorkDistributor Failed due to non-contiguous work distribution\n");
            printf("Soln: Reduce number of threads");
            exit(3);
    }

    return;
}







