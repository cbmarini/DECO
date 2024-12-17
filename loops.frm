#-
Off Statistics;
#include- DECO.h

*
*	User input
*	
#define MASSDIM "5"
#define EOM "1"
#define IBP "1"
#define numFermGen "1"

Symbol h,hd,h2,h2d,q,qd,u,ud,d,dd,l,ld,e,ed,B,W,G,p;

#do SU2rep = 1,6
	Local Input = 
		#include SM.inc
		+ Scalar(h2,SU2(`SU2rep'),U1(3))
		+ Scalar(h2d,SU2(`SU2rep'),U1(-3))
		;
	.sort

	#call HilbertSeries(p)
	#write "SU(2) rep of second Higgs field is `SU2rep':\n	HS = %E",HS
	.sort
#enddo
.end



