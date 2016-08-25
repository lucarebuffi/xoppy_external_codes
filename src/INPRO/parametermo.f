	subroutine parametermo(z,ele,temp)
c	
c	This subroutine assigns the atomic symbol of the elements to 
c	the choice of the crystal made in the main program and provides
C	the necessary crystal parameters for further calculations.
c	
c	Subroutine written by Michael Krisch
c	last interaction:  12.2.92
c	**Modified for the INPRO program by Carlos Giles (21/1/93)
c	Thermal expansion coefficient corrections for some crystals
c	have been implemented.
C       Modified:  Nov  2013  RJD, ASD/APS. Changed order of variables in "parameter" common block
C                       so that the integer variables now appear last to avoid
C                       alignment issues and compilation warnings.
c
	Integer*4     z
	Integer*4     ne,nue,moz,str
C
	Real*8	      zmass,par,pcr,tdeb,temp
C
C
	Character*2   ele
C
	Dimension     ele(2),moz(2),nue(2),zmass(2),tdeb(2)
c
	Common/parameter/zmass,par,pcr,tdeb,ne,nue,moz,str
c
	if (z.eq.1)  then
		ele(1) = 'SI'         ! silicon
		moz(1) = 14
		ne = 1	
		nue(1) = 8
		zmass(1) = 28.0855d0
		par = 5.430941d-10*(1.d0+2.57d-06*(temp-298.d0))  !at 25 Celsius
		tdeb(1) = 532.d0 
		str = 1
	endif
C
	if (z.eq.2) then
		ele(1) = 'GE'         ! germanium
		moz(1) = 32
		ne = 1
		nue(1) = 8
		zmass(1) = 72.59d0
		par = 5.65735d-10*(1.d0+6.d-06*(temp-298.d0))  !at 25 Celsius
		tdeb(1) = 293.d0
		str = 1
	endif
C
	if (z.eq.3)  then
		ele(1) = 'C'          ! carbon   (diamond)
		moz(1) = 6
		ne = 1
		nue(1) = 8
		zmass(1) = 12.011d0
		par = 3.56679d-10*(1.d0+1.18d-06*(temp-293.d0))  ! at 20 celsius
		tdeb(1) = 2021.d0
		str = 1
	endif
C
	if (z.eq.4)  then
		ele(1) = 'GA'         ! gallium arsenite
		ele(2) = 'AS'    
		moz(1) = 31
		moz(2) = 33
		ne = 2
		nue(1) = 4
		nue(2) = 4
		zmass(1) = 69.72d0
		zmass(2) = 74.9216d0
		par = 5.6537d-10
		tdeb(1) = 182.5d0
		tdeb(2) = 254.5d0
		str = 1
	endif
C
	if (z.eq.5)  then
		ele(1) = 'GA'         ! gallium phosphite
		ele(2) = 'P'  
		moz(1) = 31
		moz(2) = 15
		ne = 2
		nue(1) = 4
		nue(2) = 4
		zmass(1) = 69.72d0
		zmass(2) = 30.97376d0
		par = 5.4505d-10
		tdeb(1) = 182.5d0
		tdeb(2) = 200.d0          ! estimated value
		str = 1
	endif
C
	if (z.eq.6)  then
		ele(1) = 'IN'         ! indium arsenite
		ele(2) = 'AS'    
		moz(1) = 49
		moz(2) = 33
		ne = 2
		nue(1) = 4
		nue(2) = 4
		zmass(1) = 114.82d0
		zmass(2) = 74.9216d0
		par = 6.0360d-10
		tdeb(1) = 103.5d0
		tdeb(2) = 254.5d0
		str = 1
	endif
C
	if (z.eq.7)  then 
		ele(1) = 'IN'         ! indium phosphite
		ele(2) = 'P'     
		moz(1) = 49
		moz(2) = 15
		ne = 2
		nue(1) = 4
		nue(2) = 4
		zmass(1) = 114.82d0
		zmass(2) = 30.97376d0
		par = 5.8687d-10
		tdeb(1) = 103.5d0
		tdeb(2) = 200.d0          ! estimated value
		str = 1
	endif
C
	if (z.eq.8)  then
		ele(1) = 'IN'         ! indium antimonite
		ele(2) = 'SB'   
		moz(1) = 49
		moz(2) = 51
		ne = 2
		nue(1) = 4
		nue(2) = 4
		zmass(1) = 114.82d0
		zmass(2) = 121.75d0
		par = 6.4782d-10
		tdeb(1) = 103.5d0
		tdeb(2) = 170.d0
		str = 1
	endif
C
	if (z.eq.9)  then
		ele(1) = 'SI'         ! silicon carbide
		ele(2) = 'C'    
		moz(1) = 14
		moz(2) = 6
		ne = 2
		nue(1) = 4
		nue(2) = 4
		zmass(1) = 28.0855d0
		zmass(2) = 12.011d0
		par = 4.348d-10
		tdeb(1) = 532.d0
		tdeb(2) = 2021.d0
		str = 1
	endif
C
	if (z.eq.10) then
		ele(1) = 'CS'         ! cesium fluoride
		ele(2) = 'F'    
		moz(1) = 55
		moz(2) = 9
		ne = 2
		nue(1) = 4
		nue(2) = 4
		zmass(1) = 132.9054d0
		zmass(2) = 18.998403d0
		par = 6.008d-10
		tdeb(1) = 184.d0          ! Debye temperature for CsF
		tdeb(2) = 184.d0          ! Debye temperature for CsF
		str = 2
	endif
C
	if (z.eq.11) then
		ele(1) = 'K'          ! potassium chloride
		ele(2) = 'CL'    
		moz(1) = 19
		moz(2) = 17
		ne = 2
		nue(1) = 4
		nue(2) = 4
		zmass(1) = 39.0983d0
		zmass(2) = 35.453d0
		par = 6.29294d-10*(1.d0+36.9d-06*(temp-283.d0))
		tdeb(1) = 218.d0          ! Debye temperature of KCl
		tdeb(2) = 218.d0          ! Debye temperature of KCl
		str = 2
	endif
C
	if (z.eq.12) then
		ele(1) = 'LI'         ! lithium fluoride
		ele(2) = 'F'     
		moz(1) = 3
		moz(2) = 9
		ne = 2
		nue(1) = 4
		nue(2) = 4
		zmass(1) = 6.941d0
		zmass(2) = 18.998403d0
		par = 4.0271d-10*(1.d0+32.9d-06*(temp-283.d0))
		tdeb(1) = 650.d0          ! Debye temperature of LiF
		tdeb(2) = 650.d0          ! Debye temperature of LiF
		str = 2
	endif
C
	if (z.eq.13) then
		ele(1) = 'NA'         ! sodium chloride
		ele(2) = 'CL'    
		moz(1) = 11
		moz(2) = 17
		ne = 2
		nue(1) = 4
		nue(2) = 4
		zmass(1) = 22.98977d0
		zmass(2) = 35.453d0
		par = 5.63978d-10*(1.d0+39.5d-06*(temp-283.d0))  ! at 18 celsius
		tdeb(1) = 270.d0          ! Debye temperature of NaCl
		tdeb(2) = 270.d0          ! Debye temperature of NaCl
		str = 2
	endif
C
	if (z.eq.14) then 
		ele(1) = 'C'          ! carbon   (graphite)
		moz(1) = 6
		ne = 1
		nue(1) = 4
		zmass(1) = 12.011d0
		par = 2.456d-10
		pcr = 6.696d-10
		tdeb(1) = 2021.d0
		str = 3
	endif
C
	if (z.eq.15) then
		ele(1) = 'BE'         ! beryllium
		moz(1) = 4
		ne = 1
		nue(1) = 2
		zmass(1) = 9.01218d0
		par = 2.2866d-10
		pcr = 3.5833d-10      ! at about 22 Celsius
		tdeb(1) = 1188.d0
		str = 4
	endif
C
	end	
