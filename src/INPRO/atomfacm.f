C
      subroutine atomfacm (noz,f)
C     -------------------------- 
C     written by Stefan Joksch  ( modified by Michael Krisch )
C
C     Calculation of the atomic form factors from Ferguson et al., Computer
C     Physics Commun. no.5,328(1973) 
C
C	Last interaction: 13.2.92
C	Modified by Carlos Giles (4/12/92) for use with INPRO program.
C
      implicit none
C
      Integer*4	noz!,nf		! nf is never used
      Real*8	lambda,thetab,f,a,sil,pi!,zmass	! zmass is never used
      Dimension	a(8)
      Parameter	(pi=3.141592741012573)
C
	Common/general/lambda,thetab
C
	if (noz.eq.3) goto 3
	if (noz.eq.4) goto 4
	if (noz.eq.6) goto 6
	if (noz.eq.8) goto 8
	if (noz.eq.9) goto 9
	if (noz.eq.11) goto 11
	if (noz.eq.13) goto 13
	if (noz.eq.14) goto 14
	if (noz.eq.15) goto 15
	if (noz.eq.17) goto 17
	if (noz.eq.19) goto 19
	if (noz.eq.31) goto 31
	if (noz.eq.32) goto 32
	if (noz.eq.33) goto 33
	if (noz.eq.49) goto 49
	if (noz.eq.51) goto 51
	if (noz.eq.79) goto 79
	if (noz.eq.82) goto 82
	goto 110
3	continue                          ! lithium 1+
      a(1)=+1.9967634d+00
      a(2)=+2.5303268d-01
      a(3)=-1.0037103d+01
      a(4)=+1.5826518d+01
      a(5)=-1.0373381d+01
      a(6)=+2.6859613d+00
      a(7)=+7.3697077d-02
      a(8)=-1.0669516d-01
	goto 100
4	continue                          ! beryllium
      a(1)=+4.0670396d+00
      a(2)=-5.4888561d+00
      a(3)=-8.2941237d+01
      a(4)=+4.8108942d+02
      a(5)=-1.0860257d+03
      a(6)=+1.2108034d+03
      a(7)=-6.6365824d+02
      a(8)=+1.4283028d+02
	goto 100
6	continue                          ! carbon
      a(1)=+6.0924194d+00
      a(2)=-5.3597449d+00
      a(3)=-7.5934085d+01
      a(4)=+2.8441856d+02
      a(5)=-4.2670810d+02
      a(6)=+3.2007368d+02
      a(7)=-1.1890875d+02
      a(8)=+1.7425608d+01
	goto 100
8       continue                          ! oxygen 
      a(1)=+7.9844418d+00
      a(2)=+2.2363734d+00
      a(3)=-1.3357587d+02
      a(4)=+4.2344496d+02
      a(5)=-6.1400534d+02
      a(6)=+4.7106806d+02
      a(7)=-1.8499087d+02
      a(8)=+2.9196688d+01
        goto 100
9	continue                          ! fluor 1-
      a(1)=+1.0063515d+01
      a(2)=-8.4127753d-01
      a(3)=-1.2557989d+02
      a(4)=+3.6464358d+02
      a(5)=-4.6767656d+02
      a(6)=+3.1259040d+02
      a(7)=-1.0614832d+02
      a(8)=+1.4468730d+01
	goto 100
11	continue                          ! sodium+
      a(1)=+9.9790919d+00
      a(2)=+1.9325372d+00
      a(3)=-7.5737270d+01
      a(4)=+1.5529564d+02
      a(5)=-1.4417644d+02
      a(6)=+6.9901699d+01
      a(7)=-1.7009753d+01
      a(8)=+1.6106931d+00
	goto 100
13	continue                          ! aluminium
      a(1)=+1.3156815d+01
      a(2)=-1.4342490d+01
      a(3)=-8.1559757d+01
      a(4)=+4.1947532d+02
      a(5)=-3.4821404d+02
      a(6)=+8.5331688d+02
      a(7)=-4.2050136d+02
      a(8)=+8.0870938d+01
	goto 100
14	continue                          ! silicon
      a(1)=+1.4160424d+01
      a(2)=-1.2171781d+01
      a(3)=-1.3369605d+02
      a(4)=+6.3456326d+02
      a(5)=-1.2328922d+03
      a(6)=+1.1990638d+03
      a(7)=-5.7440945d+02
      a(8)=+1.0792687d+02
	goto 100
15	continue                          ! phosphorus
      a(1)=+1.5143130d+01
      a(2)=-9.0363137d+00
      a(3)=-1.8233776d+02
      a(4)=+8.0233940d+02
      a(5)=-1.4912678d+03
      a(6)=+1.4010514d+03
      a(7)=-6.5290629d+02
      a(8)=+1.1998103d+02
	goto 100
17	continue                          ! chlorine-
      a(1)=+1.8303162d+01
      a(2)=-1.3848861d+01
      a(3)=-1.7747297d+02
      a(4)=+6.6688997d+02
      a(5)=-1.0093129d+03
      a(6)=+7.6116324d+02
      a(7)=-2.8303522d+02
      a(8)=+4.1373777d+01
	goto 100
19	continue                          ! potassium +1
      a(1)=+1.8101457d+01
      a(2)=-7.0293778d-01
      a(3)=-1.9873972d+02
      a(4)=+6.1264439d+02
      a(5)=-8.3537663d+02
      a(6)=+5.8755644d+02
      a(7)=-2.0790667d+02
      a(8)=+2.9309386d+01
	goto 100
31	continue                          ! gallium
      a(1)=+3.1129586d+01
      a(2)=-1.1428730d+01
      a(3)=-1.8245992d+02
      a(4)=+5.9271312d+02
      a(5)=-9.4919163d+02
      a(6)=+8.4400564d+02
      a(7)=-3.8943384d+02
      a(8)=+7.2214989d+01
	goto 100
32	continue                          ! germanium
      a(1)=+3.2140100d+01
      a(2)=-1.5317000d+01
      a(3)=-1.6431000d+02
      a(4)=+5.4906000d+02
      a(5)=-8.8343000d+02
      a(6)=+7.8016000d+02
      a(7)=-3.5529000d+02
      a(8)=+6.4858000d+01
	goto 100
33	continue                           ! arsenic
      a(1)=+3.3142400d+01
      a(2)=-1.3531000d+01
      a(3)=-2.0634000d+02
      a(4)=+7.3262000d+02
      a(5)=-1.2292000d+03
      a(6)=+1.1041000d+03
      a(7)=-5.0444000d+02
      a(8)=+9.1829000d+01
	goto 100
49	continue                          ! indium
      a(1)=+4.9163313d+01
      a(2)=-9.4541647d+00
      a(3)=-3.9204341d+02
      a(4)=+1.3545237d+03
      a(5)=-2.1924577d+03
      a(6)=+1.8936612d+03
      a(7)=-8.3688468d+02
      a(8)=+1.4853842d+02
	goto 100
51	continue                          ! antimonite
      a(1)=+5.1165712d+01
      a(2)=-9.4387877d+00
      a(3)=-4.0391082d+02
      a(4)=+1.3898317d+03
      a(5)=-2.2444440d+03
      a(6)=+1.9356132d+03
      a(7)=-8.5451791d+02
      a(8)=+1.5155529d+02
	goto 100
79	continue                           ! GOLD
      a(1)=+7.9189972d+01
      a(2)=-1.0962515d+01
      a(3)=-5.3401084d+02
      a(4)=+1.7585190d+03
      a(5)=-2.7721977d+03
      a(6)=+2.3559558d+03
      a(7)=-1.0301140d+03
      a(8)=+1.8147029d+02
	goto 100
82	continue                           ! LEAD
      a(1)=+8.2192030d+01
      a(2)=-1.1078163d+01
      a(3)=-5.4710693d+02
      a(4)=+1.7964753d+03
      a(5)=-2.8287392d+03
      a(6)=+2.4029996d+03
      a(7)=-1.0506179d+03
      a(8)=+1.8510322d+02
	goto 100
100	continue
      sil=dsin(thetab)/lambda*1.0d-10
	f=a(1)+a(2)*sil+a(3)*sil**2+a(4)*sil**3+a(5)*sil**4+
     1a(6)*sil**5+a(7)*sil**6+a(8)*sil**7
	goto 1000
110	continue
c      write (*,*) ('  I can not calculate f0 for this material !')
1000	return
      end
