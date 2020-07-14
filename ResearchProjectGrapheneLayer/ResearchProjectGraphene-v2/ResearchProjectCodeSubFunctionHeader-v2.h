#ifndef RESEARCHPROJECTCODESUBFUNCTIONHEADER_H_INCLUDED
#define RESEARCHPROJECTCODESUBFUNCTIONHEADER_H_INCLUDED
double RhoJacobi(int ROWS,int COLS,double HEIGHT,double RADIUS,double HStep,double RStep);
double OmegaFunc(double Iteration,double Omega,double Rho);
void ThreadWorkDistributor(int ROWS, int THREADNUM,int **ThreadInfo);
void *CHEBYJACOBIThreadKernel(void *arguments);
int CHEBYJACOBI(int ROWS,int COLS,double HEIGHT,double RADIUS,double TOLERANCE,int THREADNUM);
struct arg_struct{
    int **ArrayThread;
    double **ArrayNew;
    double **ArrayOld;
    double *ConvList;
    int ROWS;
    int COLS;
    double RStep;
    double HStep;
    double RADIUS;
    double HEIGHT;
    double Rho;
    double Tolerance;
    double Iteration;
    int ThreadNum;
    void *barrier;
    };

#endif // RESEARCHPROJECTCODESUBFUNCTIONHEADER_H_INCLUDED
