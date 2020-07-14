#include <stdio.h>
#include <stdlib.h>
#include <math.h>



double Dim1GrapheneLayerUpdateFunc(double VAbove,double VBelow,double HStep,double epsilon2,double epsilon1,double epsilon0,double gamma,double e){
    double UpdateValue,a,b,c;
    a = (HStep*pow(e,3))/(M_PI*pow(gamma,2)*(epsilon1+epsilon2)*epsilon0);
    b = -1;
    c = (epsilon1*VAbove+epsilon2*VBelow)/(epsilon1+epsilon2);
    UpdateValue = (-b-sqrt(fabs(pow(b,2) - 4*a*c)))/(2*a);
    return UpdateValue;
}

double Dim1GeneralUpdateFunc(double VAbove, double VBelow,double epsilon2,double epsilon1){
    double UpdateValue;
    UpdateValue = (epsilon1*VAbove+epsilon2*VBelow)/(epsilon1+epsilon2);
    return UpdateValue;
}


double Dim2GrapheneLayerUpdateFunc(double VAbove,double VBelow,double VLeft,double VRight,int j,double HStep,double epsilon2,double epsilon1,double epsilon0,double gamma,double e){
    double UpdateValue,a,b,c;
    a = (HStep*pow(e,3))/(M_PI*pow(gamma,2)*(epsilon1+epsilon2)*epsilon0);
    b = -1;
    c = (epsilon1*VAbove+epsilon2*VBelow)/(epsilon1+epsilon2);
    UpdateValue = (-b-sqrt(fabs(pow(b,2) - 4*a*c)))/(2*a);
    return UpdateValue;
}

double Dim2GeneralUpdateFunc(double VAbove,double VBelow,double VLeft,double VRight,int j){
    double UpdateValue;
    double Factor1,Factor2;
    Factor1 = 1 + 1/((double)(2*j));
    Factor2 = 1 - 1/((double)(2*j));
    UpdateValue = 0.25*(Factor1*VRight+Factor2*VLeft+VAbove+VBelow);
    return UpdateValue;
}


