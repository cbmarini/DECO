#-
Off Statistics;
#include- HilbertSeries.h

#define MASSDIM "6"
#define EOM "1"
#define IBP "1"
#define numFermGen "1"


*
*	User input
*	
Symbol Phi,Phid,DeL,DeLd,DeR,DeRd,QL,QLd,QR,QRd,LL,LLd,LR,LRd,B,WL,WR,G,p;

Local Input = Scalar(Phi,SU2(2,2),U1(0))
	+ Scalar(Phid,SU2(2,2),U1(0))
	+ Scalar(DeL,SU2(3,1),U1(6))
	+ Scalar(DeLd,SU2(3,1),U1(-6))
	+ Scalar(DeR,SU2(1,3),U1(6))
	+ Scalar(DeRd,SU2(1,3),U1(-6))
	+ LHFermion(QL,SU3(3),SU2(2,1),U1(1))
	+ RHFermion(QLd,SU3(B),SU2(2,1),U1(-1))
	+ RHFermion(QR,SU3(3),SU2(1,2),U1(1))
	+ LHFermion(QRd,SU3(B),SU2(1,2),U1(-1))
	+ LHFermion(LL,SU2(2,1),U1(-3))
	+ RHFermion(LLd,SU2(2,1),U1(3))
	+ RHFermion(LR,SU2(1,2),U1(-3))
	+ LHFermion(LRd,SU2(1,2),U1(3))
	+ FieldStrength(B)
	+ FieldStrength(WL,SU2(3,1))
	+ FieldStrength(WR,SU2(1,3))
	+ FieldStrength(G,SU3(8));
.sort

#call HilbertSeries(p)
Print +s;
.sort

#call counting(HS)
.end



