*
*	Declarations
*
#define DEFScalar "0"
#define DEFFermion "0"
#define DEFLHFermion "0"
#define DEFRHFermion "0"
#define DEFDiracFermion "0"
#define DEFFieldStrength "0"
#define DEFGravity "0"

CFunction SU3, SU2, U1, U1R, S4, A4, Zn;
CFunction Scal,Scalar,Grav,Gravity,FIELD,fill,f,FieldStrength,FS,R,DF;
Autodeclare CFunction char,spin,LH,RH,Dirac;
Autodeclare Symbol x,y,z,rep,w,im,var,r,p;
Symbol i,ii,jj,field,rep,mass,n,B,maxn;




#procedure HilbertSeries(p)
*******************************************************
*
* Main procedure that runs the show
* 	A lot of the computations are outsourced 
*	to separate procedures 
*
*******************************************************
.sort

*
*	Variable that keeps track of mass dimension
*
Symbol mass(:{2*'MASSDIM'});
.sort

*
*	Process user input
*
#do field={Scalar,LHFermion,RHFermion,DiracFermion,FieldStrength,Gravity}
	if (match(`field'(?x))) redefine DEF`field' "1";
#enddo
if (match(LHFermion(?x)) || match(RHFermion(?x)) || match(DiracFermion(?x))) redefine DEFFermion "1";
.sort

#do field={Scalar,LHFermion,RHFermion,DiracFermion,FieldStrength,Gravity}
        Hide;
        Local `field'Input = Input;
        if ( count(`field',1) == 0 ) Discard;
        .sort
#enddo
Unhide;
Drop Input;
.sort

#do field={Scalar,LHFermion,RHFermion,DiracFermion,FieldStrength,Gravity}
	id `field'(x?,?x1,SU2(?rep),?x2) = `field'(x,?x1,?x2)*charSU2(?rep,n);
	id `field'(x?,?x1,SU3(?rep),?x2) = `field'(x,?x1,?x2)*charSU3(?rep,n);
	id `field'(x?,?x1,U1(?rep),?x2) = `field'(x,?x1,?x2)*charU1(?rep,n);
	id `field'(x?,?x1,U1R(?rep),?x2) = `field'(x,?x1,?x2)*charU1R(?rep,n);
	id `field'(x?,?x1,A4(?rep),?x2) = `field'(x,?x1,?x2)*charA4(?rep,n);
	id `field'(x?,?x1,S4(?rep),?x2) = `field'(x,?x1,?x2)*charS4(?rep,n);
	id `field'(x?,?x1,Zn(?rep),?x2) = `field'(x,?x1,?x2)*charZn(?rep,n);
	id `field'(x?) = `field'(x,n);
#enddo
.sort

*
*	Internal notation for the groups and representations
*
#do group={U1,SU2,SU3,A4,S4}
	id char`group'(?rep,n) = char`group'tmp(R(?rep),n);
	repeat id char`group'tmp(?x1,R(x?,?rep),n) = char`group'tmp(?x1,R(x),R(?rep),n);
	id char`group'tmp(?x1,R,n) = char`group'tmp(?x1,n);
#enddo
id charU1R(?x) = charU1Rtmp(?x);
id charZn(?x) = charZntmp(?x);
.sort
*
*	Construct the plethystic exponentials
*	This is done in a separate procedure
*
#call momentumP

#do field={Scalar,FieldStrength,Fermion,Gravity}
	#if (`DEF`field'')
		#call expandPE(`field')
	#else
		Local `field'PE = 1;
	#endif
	Brackets+ mass;
	.sort
	Hide;
#enddo

#if (`IBP')
        #call expandPE(Ibp)
#else
        Local IbpPE = 1;
#endif

Brackets+ mass;
.sort

*
*	Multiply the PE's
*
Hide;
#call multiplyPE

*
*	Taking residues for the Lorentz group
*
#call SO31symmetry

*
*	Fill in field content and their representations
*
#do ii=1,1
	id, once Scal(jj?) = ScalarInput*replace_(n,jj);
	if (match(Scal(jj?))) redefine ii "0";
	.sort
#enddo
Drop ScalarInput;
id Scalar(field?,n?) = field^n;
.sort: Scalar input inserted;

#do ii=1,1
	id, once LHF(jj?) = LHFermionInput*replace_(n,jj);
	if (match(LHF(jj?))) redefine ii "0";
	.sort
#enddo
Drop LHFermionInput;
id LHFermion(field?,n?) = field^n;
.sort: LHFermion input inserted;

#do ii=1,1
	id, once RHF(jj?) = RHFermionInput*replace_(n,jj);
	if (match(RHF(jj?))) redefine ii "0";
	.sort
#enddo
Drop RHFermionInput;
id RHFermion(field?,n?) = field^n;
.sort: RHFermion input inserted;

#do ii=1,1
	id, once DF(jj?) = DiracFermionInput*replace_(n,jj);
	if (match(DF(jj?))) redefine ii "0";
	.sort
#enddo
Drop DiracFermionInput;
id DiracFermion(field?,n?) = field^n;
.sort: DiracFermion input inserted;

#do ii=1,1
	id, once FS(jj?) = FieldStrengthInput*replace_(n,jj);
	if (match(FS(jj?))) redefine ii "0";
	.sort
#enddo
Drop FieldStrengthInput;
id FieldStrength(field?,n?) = field^n;
.sort: FieldStrength input inserted;

#do ii=1,1
	id, once Grav(jj?) = GravityInput*replace_(n,jj);
	if (match(Grav(jj?))) redefine ii "0";
	.sort
#enddo
Drop GravityInput;
id Gravity(field?,n?) = field^n;
.sort: Gravity input inserted;

*
*	Remaining symmetries
*
#$terminate = 0;

#do group={U1,U1R,SU2,SU3,Zn,A4,S4}
    #do ii=1,1
	id char`group'tmp(R(?rep),?x,n?) = char`group'(?rep,n)*char`group'tmp(?x,n);
	id char`group'tmp(n?) = 1;
	#call `group'symmetry
	
	if (match(char`group'tmp(?x))) redefine ii "0";
	.sort: `group' symmetries done;
    #enddo
#enddo

#endprocedure






#procedure momentumP
*******************************************************
*
* Construction of the momentum generating function
*
*******************************************************
Table Momentum(1:'MASSDIM');

Local P = 0 +
#do l=1,'MASSDIM'
        + fill('l')*sum_(ii,0,'MASSDIM',f('l')^ii/fac_(ii))
#enddo
;
.sort

#do kk=1,1
        id, once mass^ii?*f(n?) 
		= mass^ii*sum_(jj,1,'MASSDIM','p'^(n*jj)*mass^(2*n*jj)*charSO31(R(1/2,1/2),n*jj)/jj);
        if ( count(f,1) ) redefine kk "0";
        .sort
#enddo

Brackets+ fill;
.sort

Fillexpression Momentum = P(fill);
Drop P;
.sort: momentum generating function done;
#endprocedure







#procedure expandPE(field)
*******************************************************
*
* Procedure that expands the Plethystic exponentials
*
*******************************************************
#switch `field'
	#case Scalar
	Local ScalarArgumentPE 
		= sum_(ii,1,'MASSDIM',mass^(2*ii)*Scal(ii)*(1-`EOM'*`p'^(2*ii)*mass^(2*2*ii))*Momentum(ii)/ii);
	#break
        #case FieldStrength
        Local FieldStrengthArgumentPE 
		= sum_(ii,1,'MASSDIM',mass^(2*2*ii)*FS(ii)*(charLorentzSpin1(ii)-'EOM'*2*'p'^ii*mass^(2*ii)*charSO31(R(1/2,1/2),ii)+'EOM'*2*'p'^(2*ii)*mass^(2*2*ii))*Momentum(ii)/ii);
        #break
        #case Fermion
        Local FermionArgumentPE 
		= sum_(ii,1,'MASSDIM','numFermGen'*mass^(3*ii)*(-1)^(ii+1)*(charFermion(ii) - 'EOM'*'p'^ii*mass^(2*ii)*charFermionEOM(ii))*Momentum(ii)/ii);
        #break
        #case Gravity
        Local GravityArgumentPE 
		= sum_(i,1,'MASSDIM',mass^(2*2*i)*Grav(i)*(spin2(i)-'EOM'*'p'^i*mass^(2*i)*spin32(i)+'EOM'*'p'^(2*i)*mass^(2*2*i)*charLorentzSpin1(i))*Momentum(i)/i);
        #break
        #case Ibp
        Local IbpArgumentPE 
		= sum_(ii,1,'MASSDIM', -'p'^ii*mass^(2*ii)*charSO31(R(1/2,1/2),ii)/ii);
        #break
#endswitch
Brackets+ mass;
.sort

Hide;
Local `field'PE = sum_(ii,0,'MASSDIM',f^ii/fac_(ii));
#do kk=1,1
	id, once f*mass^n? = sum_(ii,1,{2*'MASSDIM'}-n, mass^ii*FIELD(ii))*mass^n;
	if ( count(f,1,mass,1) > {2*'MASSDIM'} ) Discard;
	if( match(f) ) redefine kk "0";
	id FIELD(n?) = `field'ArgumentPE[mass^n];
	.sort
#enddo

Unhide;
Drop `field'ArgumentPE;
.sort: `field' PE expanded; 
#endprocedure






#procedure multiplyPE
*******************************************************
*
* Constructing the HS at mass dimension
*    `MASSDIM' by multiplying the PE's
*
*******************************************************
Local HS = sum_(i,0,2*'MASSDIM',mass^i*FieldStrengthPE[mass^i]*f(i));
.sort

id f(n?) = sum_(i,0,2*'MASSDIM'-n,FermionPE[mass^i]*mass^i);
Brackets+ mass;
.sort

*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Local HS = sum_(i,0,2*'MASSDIM',mass^i*HS[mass^i]*f(i));
.sort

id f(n?) = sum_(i,0,2*'MASSDIM'-n,GravityPE[mass^i]*mass^i);
Brackets+ mass;
.sort

*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Local HS = sum_(i,0,2*'MASSDIM',mass^i*HS[mass^i]*f(i));
.sort

id f(n?) = sum_(i,0,2*'MASSDIM'-n,ScalarPE[mass^i]*mass^i);
Brackets+ mass;
.sort

*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Local HS = sum_(i,0,2*'MASSDIM',HS[mass^i]*f(i));
.sort

id f(n?) = IbpPE[mass^(2*'MASSDIM'-n)];
.sort: generating function expanded;

Unhide;
Drop IbpPE;
Drop ScalarPE;
Drop GravityPE;
Drop FieldStrengthPE;
Drop FermionPE;
.sort
#endprocedure






#procedure SO31symmetry
**************************************************************************
*       Lorentz symmetry
*
*       All Lorentz characters are inserted.
*       - Note: charcter for scalars is just equal to 1, so no need to
*       insert this seperately.
*       - Use as basic building blocks the (1/2,0) and (0,1/2) reps via
*       their tensor products as much as possible to obtain speed up.
*       Only for spin (2,0)+(0,2) this is not a speed up.
*       - Insert the character once at a time with id, once s.t. Form can
*       sort the expression as much as possible in intermediate steps.
*       - First insert the left handed characters (with variable y1) and
*       than take the integral over y1. Some terms will be set to zero
*       before we have to insert the right handed characters (variable y2)
*       and integrate/take residues.
**************************************************************************

*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*       Characters for derivative
*       Use: (1/2,0) * (0,1/2) = (1/2,1/2)
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#do i=1,1
        id, once charSO31(R(1/2,1/2),n?) = charSO31(R(1/2,0),n)*charSO31(R(0,1/2),n);
        if ( match(charSO31(R(1/2,1/2),n?)) ) redefine i "0";
        .sort
#enddo

*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*       Characters for left handed, right handed and Dirac fermions
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#do i=1,1
        id, once charFermion(n?) 
		= 'DEFLHFermion'*charSO31(R(1/2,0),n)*LHF(n) 
		+ 'DEFRHFermion'*charSO31(R(0,1/2),n)*RHF(n) 
		+ 'DEFDiracFermion'*(charSO31(R(1/2,0),n)+ charSO31(R(0,1/2),n))*DF(n);
        if ( count(charFermion,1) ) redefine i "0";
        .sort
#enddo

*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*       Characters to subtract the EOM of the left handed, right handed
*       and Dirac fermions
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#do i=1,1
        id, once charFermionEOM(n?) 
		= 'DEFRHFermion'*charSO31(R(1/2,0),n)*RHF(n)
		+ 'DEFLHFermion'*charSO31(R(0,1/2),n)*LHF(n) 
		+ 'DEFDiracFermion'*(charSO31(R(0,1/2),n) + charSO31(R(1/2,0),n))*DF(n);
        if ( count(charFermionEOM,1) ) redefine i "0";
        .sort
#enddo


*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*       FieldStrength Lorentz characters
*       Use that we can express spin 1 in terms of spin 1/2
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#do i=1,1
        id, once charLorentzSpin1(n?) = charSO31(R(1/2,0),2*n) + charSO31(R(0,1/2),2*n) + 2;
        if ( count(charLorentzSpin1,1) ) redefine i "0";
        .sort
#enddo

*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*       Characters to subtract the EOM of gravity/Weyl tensor
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#do i=1,1
        id, once spin32(n?) 
		= (charSO31(R(1/2,0),3*n) + charSO31(R(1/2,0),n))*charSO31(R(0,1/2),n) 
		+ charSO31(R(1/2,0),n)*(charSO31(R(0,1/2),3*n) + charSO31(R(0,1/2),n));
        if ( count(spin32,1) ) redefine i "0";
        .sort
#enddo

*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*       Filling in the explicit form of the (Left) characters
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#do i=1,1
        id, once charSO31(R(1/2,0),n?) = y1^n+1/y1^n;
        if ( match(charSO31(R(1/2,0),n?)) ) redefine i "0";
        .sort
#enddo

#do i=1,1
        id, once spin2(n?) = y1^(4*n)+y1^(2*n)+1+1/y1^(2*n)+1/y1^(4*n) + spin2Right(n);
        if ( count(spin2,1) ) redefine i "0";
        .sort
#enddo


*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*       Residues of Lorentz group
*       See Tab. V for the Haar measure.
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Brackets+ y1;
.sort

Local HS = HS[1] - HS[y1^(-2)];
.sort

*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*       Filling in the explicit form of the right handed characters
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#do i=1,1
        id, once charSO31(R(0,1/2),n?) = y2^n+1/y2^n;
        if ( match(charSO31(R(0,1/2),n?)) ) redefine i "0";
        .sort
#enddo

#do i=1,1
        id, once spin2Right(n?) = y2^(4*n)+y2^(2*n)+1+1/y2^(2*n)+1/y2^(4*n);
        if ( count(spin2Right,1) ) redefine i "0";
        .sort
#enddo

Brackets+ y2;
.sort

Local HS = HS[1] - HS[y2^(-2)];
.sort

#endprocedure










#procedure SU2symmetry
*******************************************************
*
*       SU(2) symmetry
*
*******************************************************
id charSU2(1,n?) = 1;

#do i=1,1
        id, once charSU2(3,n?) = charSU2(2,2*n) + 1;
        if ( match(charSU2(3,n?)) ) redefine i "0";
        .sort
#enddo

#do i=1,1
        id, once charSU2(2,n?) = y^n + 1/(y^n);
        if ( match(charSU2(2,n?)) ) redefine i "0";
        .sort
#enddo

#call terminate(SU2)

*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*       Residue of SU(2)
*       See Tab. V for the Haar measure.
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Brackets+ y;
.sort

Local HS = HS[1] - HS[y^(-2)];
.sort
#endprocedure









#procedure U1symmetry
*******************************************************
*
*	U(1) symmetry
*
*******************************************************
id charU1(rep?int_,n?) = x^(rep*n);
.sort

#call terminate(U1)

Brackets+ x;
.sort

Local HS = HS[1];
.sort
#endprocedure








#procedure SU3symmetry
*******************************************************
*
*       SU(3) symmetry
*
*******************************************************
id charSU3(1,n?) = 1;
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*       Character of fundamental representation (3)
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#do i=1,1
        id, once charSU3(3,n?) = 1/z2^n + z2^n/z1^n + z1^n;
        if ( match(charSU3(3,n?)) ) redefine i "0";
        .sort
#enddo

*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*       Character of anti fundamental representation (B)
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#do i=1,1
        id, once charSU3(B,n?) = 1/z1^n + z1^n/z2^n + z2^n;
        if ( match(charSU3(B,n?)) ) redefine i "0";
        .sort
#enddo

*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*      Character of adjoint representation (8)
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#do i=1,1
       id, once charSU3(8,n?) = z1^n*z2^n+ z2^(2*n)/z1^n + z1^(2*n)/z2^n + 2 + z1^n/z2^(2*n) + z2^n/z1^(2*n) + 1/(z1^n*z2^n);
       if ( match(charSU3(8,n?)) ) redefine i "0";
       .sort
#enddo

#call terminate(SU3)


*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*       Residue of SU(3)
*       See Tab. V for the Haar measure.
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Brackets+ z1;
.sort

Local HS = HS[1]/z2 - HS[z1^1]*z2 - HS[z1^(-2)]/z2^2 + HS[1]*z2^2 + HS[z1^(-3)]/z2 - HS[z1^(-2)]*z2;
Brackets+ z2;
.sort

Local HS = HS[z2^(-1)];
.sort

#endprocedure








#procedure A4symmetry
*******************************************************
*
*       A4 symmetry
*
*******************************************************
id charA4(1,n?) = var1^n;
id charA4(p1,n?) = varp^n;
id charA4(pp1,n?) = varpp^n;
id charA4(3,n?) = var31^n + var32^n + var33^n;
.sort

#call terminate(A4)

Local HS = 1/12*(
		+ 1*HS*replace_(var1,1,varp,1  ,varpp,1  ,var31,1,var32, 1,var33, 1)
		+ 3*HS*replace_(var1,1,varp,1  ,varpp,1  ,var31,1,var32,-1,var33,-1)
		+ 4*HS*replace_(var1,1,varp,w  ,varpp,w^2,var31,1,var32,w^2,var33,w)
		+ 4*HS*replace_(var1,1,varp,w^2,varpp,w  ,var31,1,var32,w,var33,w^2)
	);
*
* Use w = Exp(2*pi*I/3), so
* 	w^3 = 1 and 1 + w + w^2 = 0.
*
id w^3 = 1;
id w^2 = -1-w;
.sort
#endprocedure









#procedure S4symmetry
*******************************************************
*
*       S4 symmetry
*
*******************************************************
id charS4(1,n?) = var1^n;
id charS4(p1,n?) = varp^n;
id charS4(2,n?) = var21^n+var22^n;
id charS4(3,n?) = var31^n + var32^n + var33^n;
id charS4(p3,n?) = varp31^n + varp32^n + varp33^n;
.sort

#call terminate(S4)

*
* Sum over conjugace classes
*   (in the order C1, C3,C6, C6' and C8).
*
Local HS = 1/24*(
		+ 1*HS*replace_(var1,1,varp,1 ,var21,1 ,var22,1  ,var31,1 ,var32, 1,var33,1  ,varp31,1 ,varp32,1  ,varp33,1)
		+ 3*HS*replace_(var1,1,varp,1 ,var21,1 ,var22,1  ,var31,1 ,var32,-1,var33,-1 ,varp31,1 ,varp32,-1 ,varp33,-1)
		+ 6*HS*replace_(var1,1,varp,-1,var21,-1,var22,1  ,var31,-1,var32,1 ,var33,1  ,varp31,-1,varp32,-1 ,varp33,1)
		+ 6*HS*replace_(var1,1,varp,-1,var21,-1,var22,1  ,var31,-1,var32,im,var33,-im,varp31,im,varp32,-im,varp33,1)
		+ 8*HS*replace_(var1,1,varp,1 ,var21,w ,var22,w^2,var31,1 ,var32,w^2,var33,w ,varp31,1,varp32,w   ,varp33,w^2)
	);
*
* Use: 
*   im^2 = -1.
*   w = Exp(2*pi*I/3), so w^3 = 1 and 1 + w + w^2 = 0.
*
id im^2 = -1;
id w^3 = 1;
id w^2 = -1-w;
.sort
#endprocedure






#procedure Znsymmetry
*******************************************************
*
*       Zn symmetry
*
*******************************************************
#$maxn = 0;
id,once charZn(maxn?$maxn,rep?int_,n?) = w^(rep*n);
id charZn($maxn,rep?int_,n?) = w^(rep*n);
.sort

#call terminate(Zn)

id w^$maxn = 1;
id w = 0;
.sort
#endprocedure







#procedure U1Rsymmetry
*******************************************************
*
*       U(1)_R symmetry
*
*******************************************************
#$rem = 0;
id,once charU1R(r?$rem,rep?int_,n?) = xr^(rep*n);
id charU1R($rem,rep?int_,n?) = xr^(rep*n);
.sort

#call terminate(U1R)

Brackets+ xr;
.sort

Local HS = HS[xr^$rem];
.sort
#endprocedure








#procedure terminate(group)
*******************************************************
*
*       This procecure checks if the charges are 
*	correctly defined by the user.
*
*******************************************************
if (match(char`group'(?x))) $terminate = 1;
.sort
	
#if (`$terminate' == 1)
	#write "Program terminated during symmetries of `group'.\nCheck user input, e.g. charges etc."
	#terminate
#endif
.sort
	
#endprocedure










#procedure counting(HS)
*****************************************************************
* This procedure counts the number of operators in the Local
*       expression Hilbert
*****************************************************************
        Local HilbertCounting = `HS';
        .sort

        Hide `HS';
        dropsymbols;
        .sort

        #$number = HilbertCounting;
        #write "Number of operators at mass dimension `MASSDIM' is `$number'."
        .sort

        Drop HilbertCounting;
        .sort

        Unhide `HS';
        .sort
#endprocedure




#procedure saveto(expr, file)
*****************************************************************
* Save a given Local expression `expr' to a Mathematica 
*	file `file'.
*****************************************************************
    format Mathematica;
    format nospaces;
    format 78;
    #write <`file'> "(\n      %E\n)" `expr';
    #close <`file'>
#endprocedure
