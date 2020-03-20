      SUBROUTINE trapzd8(func,a,b,s,n)
C  Modified by Roger J. Dejus XFD/APS, April, 1995.
C  Changed to double precision. Routine is from subroutine trapzd.f
C  from "Numerical Recipes" Section 4.2.

      INTEGER n
      REAL*8 a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL*8 del,sum,tnm,x
      if (n.eq.1) then
        s=0.5d0*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5d0*del
        sum=0.0d0
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5d0*(s+(b-a)*sum/tnm)
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software |a.
