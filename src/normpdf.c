#include <math.h>

double normpdf(double x,double u,double sgm)
{
    //printf("%f",M_PI);
    return 1.0/(sqrt(2*M_PI)*sgm)*exp(-pow((x-u)/sgm,2)/2);
}
