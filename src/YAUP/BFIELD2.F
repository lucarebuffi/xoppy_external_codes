C
C This program is the same as bfield.f provided with
C yaup with the only exceptions of the field coefficients:
C bfield uses ampl=0.95*3.44*exp(-r*(5.08-1.54*r)) which stands
C        for a Nd-Fe-B magnet and 
C bfield uses ampl=0.95*3.33*exp(-r*(5.47-1.80*r)) which stands
C        for a Nd-Fe-B magnet.
C
C Modification made by M. Sanchez del Rio, ESRF, 94/11/23
C
c
c       ---------------------------------                                      
        double precision function ampl(z)                                      
c       ---------------------------------                                      
c                                                                              
c  this is the undulator field-strength function.  the B field                
c  is assumed be an amplitude-modulated sinusoid, e.g.             
c       Btot(z) = B(z)*sin(2*pi/period*z).                              
c                                                                              
        implicit double precision (a-h,o-z)                                    
        common/param/ ulen,per,rk0,gzmin,dg                                     
c
c   the field strength is for the APS.  A linearly modulated gap
c   is assumed.  DG is the degree of gap taper.
c                                                                             
        gap=gzmin*(1.0+z*dg/ulen)                                         
        r=gap/per                                                              
        ampl=0.95*3.33*exp(-r*(5.47-1.80*r))
        return                                                                 
        end                                                                    
c
c       ----------------------------------                                      
        double precision function phase(z)                                      
c       ----------------------------------                                      
c                                                                              
c   phase errors (if any)
c                                                                              
        implicit double precision (a-h,o-z)                                    
        common/param/ ulen,per,rk0,gzmin,dg                                     
c
	phase=0.d0
	return
	end
c
c	--------------
	program bfield
c	--------------
c
c   generate a b-field file in format acceptable to YAUP
c   iit, dept of physics, bib, 7/92
c                                                                              
        implicit double precision (a-h,o-z)                                    
        parameter (maxper=150,maxpp=100,maxpts=maxper*maxpp)           
        dimension b(maxpts),phi(maxpts)                                         
        character*30 outfn                                                     
        common/param/ ulen,per,rk0,gzmin,dg                                     
c                                                                              
c   undulator papramers                                                   
c                                                                              
        write(*,10)                                     
10      format(/1x,'Generate a B-field file in YAUP format.'/)                                                         
        write (*,20)
20	format(1x,'Undulator period (cm) ? ',$)                                   
        read(*,*) per                                                            
        write(*,30) maxper
30	format(1x,'Number of periods (',i3,' max) ? ',$)                             
        read(*,*) nper                                                           
        write(*,40) maxpp
40	format(1x,'Points per period (40 sugg, ',i3,' max) ? ',$)                    
        read(*,*) npts
c
42	write(*,44)
44	format(/1x,'Enter : '/3x,'[1] for planar undulator'
     +	       /3x,'[2] for tapered undulator'
     +	       /1x,'Choose [1-2] : ',$)
	read(*,*) imode
c
	if (imode.eq.1) then
	   write(*,46)
46	   format(/1x,'Planar undulator case :'
     +	         /1x,'K-factor ? ',$)
	   read(*,*) rk
	   b0=rk/(0.934d0*per)
	else if (imode.eq.2) then                                                          
           write(*,50)
50	   format(/1x,'Tapered undulator case :'
     +	          /1x,'Gap at zmin (cm) ? ',$)                                        
           read(*,*) gzmin                                                          
           write(*,60)
60	   format(1x,'Gap taper (%) ? ',$)
           read(*,*) dg
	   dg=dg/100.
	else
	   write(*,*) char(7)
	   goto 42
	end if
c                                                         
        write(*,70)
70	format(1x,'Output filename ? ',$)                                         
        read (*,'(a)') outfn
c                                                      
        write(*,80) 
80	format(/1x,'Working...')                                              
c                                                                              
c   some constants                                                             
c                                                                              
        twopi=8.d0*datan(1.d0)                                                 
        rk0=twopi/per                                                          
        nptst=nper*npts+1                                                      
        ulen=dble(nper)*per                                                    
        dz=ulen/dble(nptst-1)                                                  
c
c   generate the values for the modulating function and the phase errors                                                                              
c   bmin is the minimum value of the modulating function                               
c                                                                              
        bmin=1.d30                                                            
        do i=1,nptst                                                           
           z=dble(i-1)*dz
	   b(i)=b0                                                   
           if (imode.eq.2) b(i)=ampl(z)
	   phi(i)=phase(z)                                                    
           if (b(i).lt.bmin) bmin=b(i)                                        
        end do                                                                 
c
c   open unformatted file and save the results.  the ascii file is for 
c   reference purposes only.
c                                                                              
        open(1,file=outfn,form='unformatted',status='unknown')                 
        write(1) per,nper,npts,bmin                                            
Csrio        write(1) (b(i),phase(dreal(i)), i=1,nptst)
        write(1) (b(i),phase(dfloat(i)), i=1,nptst)
        close(1)                                                               
c                                             
c	open(2,file='bfield.ascii',status='unknown')
c	do i=1,nptst
c	   z=dble(i-1)*dz                              
c	   write(2,90) z,b(i),phase(i),b(i)*dsin(rk0*z+phase(i))                
c90         format(1x,4(e13.4)) 
c	end do                                                
c	close(2)                                                               
c
	write(*,130) char(7)
130	format(1x,'Done.',a1)
        stop                                                                   
        end                                                                    
