C
      subroutine strufacnewm (c0r,c0i,cHr,cHi,cHrb,cHib) 
C     --------------------------------------------------     
C	written by Stefan Joksch ( modified by Michael Krisch )
C	modified and extended by Michael Krisch
C	Calculation of the components of the dielectric susceptibilities
C
C	Last interaction: 11.2.92
C	Modified by Carlos Giles for use with INPRO program (4/12/92).
C       Modified:  Nov  2013  RJD, ASD/APS. Changed order of variables in "parameter" and
C                       "structure" common blocks so that the integer variables now appear last to avoid
C                       alignment issues and compilation warnings.
C                       Fixed issue with the Debye-Waller factor variable dwf and the variable zmass not being set to
C                       dimension to 2. The variable zmass is not being used, but it caused values to be
C                       incorrect for variables follwing it in the in the common block. Hence, the variables par,
C                       pcr, tdeb, str, all had incorrect values. In particular, the str variable was incorrectly set
C                       (typically to 0 and the code continued at label 200 even for NaCl structures, which should
C                       have continued at label 300). Also the value of the first element of the dwf for a crystal
C                       compound was used. See proposed changes "CRJD ..." below to apply dwf for each element.

	Implicit none
C
	Integer*4 	hm,km,lm,ne,nue,moz,j,k,jsta,jend,str
C
	Real*4	f1,f2
	Real*8	lambda,zmass,thetab,vez,f0,dwf
	Real*8	x,y,z,fak,re,pi,u!,rho
	Real*8	par,pcr,tdeb
	Real*8	c0r,c0i
C
	Complex*16	cHr,cHi,cHrb,cHib,i
C
	Dimension	x(8),y(8),z(8)
	Dimension	f0(2),f1(2),f2(2),nue(2),moz(2)
	Dimension	zmass(2),tdeb(2),dwf(2)
C
	Parameter 	(re=2.81777d-15,pi=3.141592741012573)
	Parameter 	(u=1.6604d-24)
C
	Common/general/lambda,thetab
	Common/parameter/zmass,par,pcr,tdeb,ne,nue,moz,str
	Common/structure/vez,f0,f1,f2,dwf,hm,km,lm
C ________________________________________________________________________
C|                                                                        |
C|*** initialization and parameter setup ***                              |
C|________________________________________________________________________|
C
	i = cmplx(0.d0,1.d0)
100	continue
      c0r=0.d0
      c0i=0.d0
      cHr=0.d0
      cHrb=0.d0
      cHi=0.d0
      cHib=0.d0
      do 105 j=1,8
         x(j)=0.d0
         y(j)=0.d0
         z(j)=0.d0
105   continue
C
C	splitting up into the different crystal structures
	if (str.eq.1) goto 200     ! ZnS structure
	if (str.eq.2) goto 300     ! NaCl structure
	if (str.eq.3) goto 400     ! graphite structure
	if (str.eq.4) goto 500     ! hcp structure
C ________________________________________________________________________
C|                                                                        |
C|*** crystal structures **                                               |
C|________________________________________________________________________|
C
C	ZnS - stucture 
C	Material: Si, Ge, diamond, GaAs, GaP, InAs,InP, InSb, SiC
C
200   x(1)=0.
      x(2)=0.0
      x(3)=0.500
      x(4)=0.500
      x(5)=0.250
      x(6)=0.250
      x(7)=0.750
      x(8)=0.750
      y(1)=0.0
      y(2)=0.500
      y(3)=0.0
      y(4)=0.500
      y(5)=0.250
      y(6)=0.750
      y(7)=0.250
      y(8)=0.750
      z(1)=0.
      z(2)=0.500
      z(3)=0.500
      z(4)=0.000
      z(5)=0.250
      z(6)=0.750
      z(7)=0.750
      z(8)=0.250               
	goto 950
C	
C	sodium cloride structure
C	material: NaCl, CsF, KCl, LiF, NaCl
C
300	x(1) = 0.0
	x(2) = 0.5
	x(3) = 0.5
	x(4) = 0.0
	x(5) = 0.5
	x(6) = 0.5
	x(7) = 0.0
	x(8) = 0.0
	y(1) = 0.0
	y(2) = 0.5
	y(3) = 0.0
	y(4) = 0.5
	y(5) = 0.5
	y(6) = 0.0
	y(7) = 0.5
	y(8) = 0.0
	z(1) = 0.0
	z(2) = 0.0
	z(3) = 0.5
	z(4) = 0.5
	z(5) = 0.5
	z(6) = 0.0
	z(7) = 0.0
	z(8) = 0.5
	goto 950
C
C	graphite structure   ( hexagonal )
C	
400	x(1) = 0.0
	x(2) = 0.0
	x(3) = 0.5
	x(4) = 2.0/3.0
	y(1) = 0.0
	y(2) = 0.0
	y(3) = 2.0/3.0
	y(4) = 1.0/3.0
	z(1) = 0.0
	z(2) = 0.5
	z(3) = 0.0
	z(4) = 0.5
	goto 950
C
C	hexagonal closed packed structure
C	material: Be
C
500	x(1) = 0.0
	x(2) = 1.0/3.0
	y(1) = 0.0
	y(2) = 2.0/3.0
	z(1) = 0.0
	z(2) = 0.5
	goto 950
C
950	continue
C ________________________________________________________________________
C|                                                                        |
C|*** calculation of the susceptibilities ***                             |
C|________________________________________________________________________|
C
C	do loop for the number of different atoms in the unit cell.
C	For the moment the program is limited to the case where there are
C	equal numbers of atoms in the unit cell.
C
	do 970 k = 1,ne
 		c0r = c0r + nue(k)*(moz(k) + 1.d0*f1(k))
		c0i = c0i + nue(k)*f2(k)
C		start and end values for the do loop thru atoms in unit 
C		cell.
C
		jsta = 1 + (k - 1)*nue(k)
		jend = k*nue(k)
C		do loop thru atoms in unit cell
C
		do 955 j = jsta,jend
C
C		 	  "real" part of chih
			   cHr = cHr + ((f0(k) + 1.d0*f1(k))*
     1			   exp(2.*pi*i*(hm*x(j) + km*y(j) + lm*z(j))))
C
C		  	"imaginary" part of chih
		   	cHi = cHi + f2(k)*
     1		        exp(2.*pi*i*(hm*x(j) + km*y(j) + lm*z(j)))
C
C		   	"real part" of chih(bar)
		   	cHrb = cHrb + ((f0(k) + 1.d0*f1(k))*
     1		        exp(-2.*pi*i*(hm*x(j) + km*y(j) + lm*z(j))))
C
C		   	"imaginary part of chih(bar)
		   	cHib = cHib + f2(k)*
     1		        exp(-2.*pi*i*(hm*x(j) + km*y(j) + lm*z(j)))
C
955		continue
CRJD, dwf applied for each element
        cHr=cHr*dwf(k)
        cHrb=cHrb*dwf(k)
        cHi=cHi*dwf(k)
        cHib=cHib*dwf(k)
970	continue
C
C	Conversion from structure factor to suceptibilities
C
      fak=(-1.)*re*lambda**2/(pi*vez)
      c0r=c0r*fak
      c0i=c0i*fak
      cHr=cHr*fak
      cHrb=cHrb*fak
      cHi=cHi*fak
      cHib=cHib*fak

	return
	end              
