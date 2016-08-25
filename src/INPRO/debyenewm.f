C
      subroutine debyenewm (zmass,t,tdeb,deb)
C     --------------------------                
C
C       Calculation of the Debye Waller factor from R.W.James;
C       The optical principles of x-rays, 1965
C
C	Program written by M. Krisch, based on some ideas of St. Joksch
C
C	Last interaction: 13.2.92
C	Modified for use with matparmo.for and INPRO program,
C	Carlos Giles 4/12/92.
C
	implicit none
C
        Real*8         dm,fi
        Real*8         lambda,thetab,sid,deb,t,tdeb,tmt,hp,kb,u,zmass
C
        Parameter      (hp=6.6262d-34,kb=1.38062d-23,u=1.66053d-27)
C
	Common/general/lambda,thetab
C
C
      sid=dsin(thetab)/lambda
      tmt=tdeb/t
      dm=6.*(hp/kb*sid)*hp*sid/zmass/u/tdeb/tmt
      fi=1.+(tmt**2)/36.-(tmt**4)/3600.
      deb=dexp(-1.*dm*fi)
      return
      end      
