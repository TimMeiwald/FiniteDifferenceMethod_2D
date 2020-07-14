#ifndef UPDATEFUNCTIONS_H_INCLUDED
#define UPDATEFUNCTIONS_H_INCLUDED


// 1 Dimensional update functions
double Dim1GrapheneLayerUpdateFunc(double VAbove,double VBelow,double HStep,double epsilon2,double epsilon1,double epsilon0,double gamma,double e);
double Dim1GeneralUpdateFunc(double VAbove, double VBelow,double epsilon2,double epsilon1);


// 2 Dimensional update functions
//double Dim2GrapheneLayerUpdateFunc(double VAbove,double VBelow,double HStep,double epsilon2,double epsilon1,double epsilon0,double gamma,double e);
double Dim2GeneralUpdateFunc(double VAbove,double VBelow,double VLeft,double VRight,int j);
#endif // UPDATEFUNCTIONS_H_INCLUDED
