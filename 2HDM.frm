#-
Off Statistics;
#include- DECO.h

#define MASSDIM "6"
#define EOM "1"
#define IBP "1"
#define numFermGen "1"


*
*	User input
*	
Symbol h1,h1d,h2,h2d,q,qd,u,ud,d,dd,l,ld,e,ed,B,W,G,p;

Local Input = Scalar(h1,SU2(2),U1(3))
	+ Scalar(h1d,SU2(2),U1(-3))
	+ Scalar(h2d,SU2(2),U1(3))
	+ Scalar(h2d,SU2(2),U1(-3))
	+ LHFermion(q,SU3(3),SU2(2),U1(1))
	+ RHFermion(qd,SU3(B),SU2(2),U1(-1))
	+ RHFermion(u,SU3(3),U1(4))
	+ LHFermion(ud,SU3(B),U1(-4))
	+ RHFermion(d,SU3(3),U1(-2))
	+ LHFermion(dd,SU3(B),U1(2))
	+ LHFermion(l,SU2(2),U1(-3))
	+ RHFermion(ld,SU2(2),U1(3))
	+ RHFermion(e,U1(-6))
	+ LHFermion(ed,U1(6))
	+ FieldStrength(B)
	+ FieldStrength(W,SU2(3))
	+ FieldStrength(G,SU3(8));
.sort

#call HilbertSeries(p)
Print +s;
.sort

#call counting(HS)
.end



