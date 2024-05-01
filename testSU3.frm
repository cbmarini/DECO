#-
#include- DECO.h

#define MASSDIM "3"
#define EOM "1"
#define IBP "1"
#define numFermGen "1"

Auto Symbol a,b;

Local Input = Scalar(a1,SU3(B))
    + Scalar(a2,SU3(6))
    + Scalar(b,SU3(15));

#call HilbertSeries(0)
Print +s;
.end