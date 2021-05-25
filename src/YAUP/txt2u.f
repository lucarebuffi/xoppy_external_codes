	program txt2u
c
c   take an ASCII file and convert it to YAUP format.  The text
c   file should be column-formatted, and contain three colums : z, B(z),
c   and phi(z), where the z's are equidistant with step PERIOD/NPTS.
c   See YAUP.DOC for definitions of PERIOD and NPTS.  There should
c   be NPTS*NPER+1 lines in the ASCII file.	
c
	implicit double precision (a-h,o-z)
	parameter (maxpts=5000) 
	character fname*80,line*100
	dimension B(maxpts),phi(maxpts)
c
	write(*,10)
10	format(/1x,'txt2u - convert text b-field files to ',
     +	           'YAUP format.')
c
	ncode=1
20	write(*,30)
30	format(/1x,'Input file name ? '$)
	read(*,'(a)') fname
        open(50,file=fname,status='old',err=900)           
c
	ncode=2
40	write(*,50)
50	format(1x,'Output file name ? '$)
	read(*,'(a)') fname
        open(60,file=fname,form='unformatted',status='unknown',
     +	         err=900)           
c
	ncode=3
60	write(*,70)
70	format(1x,'Magnet period (cm) ? '$)
	read(*,*,err=900,end=900) per
c
	ncode=4
80	write(*,90)
90	format(1x,'Number of periods ? '$)
	read(*,*,err=900,end=900) nper
c
	ncode=5
100	write(*,110)
110	format(1x,'Number of points per period ? '$)
	read(*,*,err=900,end=900) npts
c
	write(*,115)
115	format(/1x,'Working...')
c
	line='#'
	do while (index(line(1:5),'#').ne.0)
	   read(50,'(a)',err=920,end=920) line
	end do
	backspace(50)
c
	nptst=npts*nper+1
	do i=1,nptst
	   read(50,*,err=920,end=920) z,B(i),phi(i)
	end do
c
	bmin=1.d30
	do i=1,nptst
	   if (B(i).lt.bmin) bmin=B(i)
	end do
c
	write(60) per,nper,npts,bmin                                    
	write(60) (B(i),phi(i),i=1,nptst)
	close(50)
	close(60)
c
	write(*,130) char(7)
130	format(1x,'Done.',a1)
	stop	
c
900	write(*,910) char(7)
910	format(/1x,'Fumble fingers!  Try again!',a1)
	goto (20,40,60,80,100) ncode
	stop
c
920	write(*,930) nptst,per/dble(npts),char(7)
930	format(/1x,
     +	'Wrong format.  The ASCII file should contain three columns:'
     +	/1x,'z, B(z) [tesla], phi(z) [radians], and ',i4,' lines.'
     +  /1x,'The z-s should be equidistant with step ',f6.3,' cm.',a1)
	close(50)
	close(60)
	stop
	end
