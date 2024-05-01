#-
#include- DECO.h

#define MASSDIM "2"
#define EOM "1"
#define IBP "1"
#define numFermGen "1"

Auto Symbol a,b;

Local Input = Scalar(a,SU3(10))
    + Scalar(b,SU3([10B]));

#call HilbertSeries(0)
Print +s;
.end