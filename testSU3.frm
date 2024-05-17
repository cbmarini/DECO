#-
#include- DECO.h

#define MASSDIM "3"
#define EOM "1"
#define IBP "1"
#define numFermGen "1"

Auto Symbol a,b;

Local Input = Scalar(a1,SU3(B3))
    + Scalar(a2,SU3(B6))
    + Scalar(b,SU3(B15));

#call HilbertSeries(0)
Print +s;
.end