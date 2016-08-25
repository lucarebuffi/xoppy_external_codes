	program u2txt
c
c   read an unformatted trajectroy or b-field file and write it to
c   an ASCII file
c
	implicit double precision (a-h,o-z)
	parameter (maxpts=5000) 
	parameter (twopi=6.2831853d0)
	character*80 fname
	dimension ct(maxpts),x(maxpts),betax(maxpts),betaz(maxpts)
	dimension bfield(maxpts),phaserr(maxpts)
	equivalence (betax,phaserrr), (betaz,bfield)
c
	write(*,10)
10	format(/1x,'u2txt - convert an unformatted YAUP file ',
     +	           'to text format.'/)
c
	ncode=1
15	write(*,20) 
20	format(1x,'Please choose file type.  Enter :'
     +	      /5x,'[1] for b-field file;'
     +	      /5x,'[2] for trajectory file;'
     +	      /1x,'Then ? '$)
	read(*,*) mode
	if (mode.lt.1.or.mode.gt.2) goto 900
c
	ncode=2
35	write(*,40)
40	format(/1x,'Input file name ? '$)
	read(*,'(a)') fname
        open(50,file=fname,form='unformatted',status='old',err=900)           
        read(50) per,nper,npts                                    
        nptst=npts*nper+1
	dz=per/dble(npts)
	rewind(50)                                                   
c
        if (nptst.gt.maxpts) then
	   write(*,41) nptst,maxpts,char(7)
41	   format(/1x,'The input file contains too many data points (',
     +	   i5,' ).  Only ',i5,' points will be transfered.',a1)
	end if 	                                        
c
	if (mode.eq.1) then
           read(50) per,nper,npts,bmin                                    
           read(50) (bfield(i),phaserr(i),i=1,nptst)                     
 	else if (mode.eq.2) then
           read(50) per,nper,npts,ener,efund                                 
           read(50) (ct(i),x(i),betax(i),betaz(i),i=1,nptst)                
	end if
        close(50)
c
	ncode=3	   
45	write(*,50) 
50	format(1x,'Output filename ? '$)
	read(*,'(a)') fname
	open(60,file=fname,status='unknown',err=900)
c
	write(*,'(/1x,a)') 'Working...'
	if (mode.eq.1) then
	   xk0=twopi/per
	   write(60,*) 
     +		'# Columns: z(cm), ampl(tesla), phserr, total(tesla)'
	   write(60,*) 
     +		'# total = ampl * sin ( twopi/period*z + phserr ) '
     	   write(60,*) 
     +		'# period=',real(per),'  nper=',nper,'  npts=',npts  
	   do i=1,nptst
	       zc=dble(i-1)*dz
	       total=bfield(i)*sin(xk0*zc+phaserr(i))
	       write(60,60) zc,bfield(i),phaserr(i),total
60	       format(1x,4(1pe15.5))
	   end do
	else if (mode.eq.2) then
	   write(60,*) 
     +		'# Columns: z(cm), ct(cm), x(cm), betax(cm), betaz(cm)'
     	   write(60,*) 
     +		'# period=',real(per),'  nper=',nper,'  npts=',npts  
	   write(60,70) ((i-1)*dz,ct(i),x(i),betax(i),betaz(i),
     +	                  i=1,nptst)
70	   format(5(1x,1pe15.5))
	end if
	close(60)
c
	write(*,80) char(7)
80	format(1x,'Done.',a1)
	stop                                 
c
900	write(*,910) char(7)
910	format(/1x,'Fumble fingers!  Try again!',a1)
	goto (15,35,45) ncode
	stop
	end
