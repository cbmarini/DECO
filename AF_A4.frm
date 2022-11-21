#-
Off Statistics;
#include- HilbertSeries.h

#define MASSDIM "7"
#define EOM "1"
#define IBP "1"
#define numFermGen "1"


*
*	User input
*	
Symbol hu,hud,hd,hdd,phiT,phiS,xi,xitilde,G,W,B,l,lb,e,eb,mu,mub,tau,taub;

Local Input = 
	+ Scalar(hu,SU2(2),U1(3),A4(1),Zn(R(3,0)))
	+ Scalar(hud,SU2(2),U1(-3),A4(1),Zn(R(3,0)))
	+ Scalar(hd,SU2(2),U1(3),A4(1),Zn(R(3,0)))
	+ Scalar(hdd,SU2(2),U1(-3),A4(1),Zn(R(3,0)))
	+ Scalar(phiT,A4(3),Zn(R(3,0)))
	+ Scalar(phiS,A4(3),Zn(R(3,1)))
	+ Scalar(xi,A4(1),Zn(R(3,1)))	
	+ Scalar(xitilde,A4(1),Zn(R(3,1)))
	+ FieldStrength(G,SU3(8))	
	+ FieldStrength(W,SU2(3))	
	+ FieldStrength(B)
	+ LHFermion(l,SU2(2),U1(-3),A4(3),Zn(R(3,1)),U1R(R(2,1)))	
	+ RHFermion(lb,SU2(2),U1(3),A4(3),Zn(R(3,2)),U1R(R(2,-1)))	
	+ RHFermion(e,U1(-6),A4(1),Zn(R(3,1)),U1R(R(2,-1)))	
	+ LHFermion(eb,U1(6),A4(1),Zn(R(3,2)),U1R(R(2,1)))	
	+ RHFermion(mu,U1(-6),A4(p1),Zn(R(3,1)),U1R(R(2,-1)))	
	+ LHFermion(mub,U1(6),A4(pp1),Zn(R(3,2)),U1R(R(2,1)))	
	+ RHFermion(tau,U1(-6),A4(pp1),Zn(R(3,1)),U1R(R(2,-1)))	
	+ LHFermion(taub,U1(6),A4(p1),Zn(R(3,2)),U1R(R(2,1)))	
	;
.sort

#call HilbertSeries(p)
Print +s;
.sort

#call counting(HS)
.sort

Brackets+ l;
.sort

Local HS= HS[l^2];
Print +s;
.sort

#call counting(HS)
.end



