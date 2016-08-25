/* ------------------------------------------------------------
   File Name:   mamspec.c
   Author:      John M. Boone, Ph.D.
   Date:        Dec 15, 1997
   Description: Mammography Spectrum Program
   
Note:  This code is supplied as is and you should always check your
       calculated results.  This program was created for ease of use, 
       and was not optimized for memory efficiency or operational speed.
       Good luck.
==================================================================
   Parameters passed to subroutine:  MAMSPEC( itube,xkv,en,spec )
INPUTS:
   itube: (integer)
              1 = molybdenum anode
              2 = rhodium anode
              3 = tungsten anode
   xkv:  (floating point)
              needs to be in the range 18.0 to 42.0
OUTPUTS:
   en:   (floating point array, at least 100 elements long)
              lists the energies of the output spectrum              
   spec: (floating point array, al least 100 elements long)
              lists the fluence values of the output spectrum
------------------------------------------------------------ */

/* standard include files (DOS) */
#include "stdio.h"
#include "math.h"
#include "stdlib.h"
/* #include "dos.h" */

/* include files which contain polynomial coefficients */
#include "moly.h"
#include "rhodium.h"
#include "tungsten.h"

int mamspec( itube,xkv,en,aspec )
int itube;
float xkv,en[],aspec[];
{
	double arg,dsum,sum1,sum2;
	float aa[150][4],mr,xval,e1,e2,energy,w1,w2;
	float aenergy[150];
	float x,y,z,xmax,xmin,sum,factor,t0,t1,t2,t3,t4,t5,t6;
	int ien,nterms,nt,nn[150],iemax,kv,ikv,genflag,iflag;
	int i,j,k,l,m,n,ians,ival,jj,j1,j2,npts;
/* ----------------------------------------------------
	If input range is outside limits, reject
---------------------------------------------------- */
	if( itube < 1 || itube > 3 ) goto ERROR_CONDITION;
	if( xkv < 18. || xkv > 42. ) goto ERROR_CONDITION;
/* ----------------------------------------------------
	Transfer Coefficients for Molybdenum Anode (itube=1)
	into working arrays
---------------------------------------------------- */
	if( itube==1 ) { /* MOLY */
		for( i=0; i<100; ++i ) {
			nn[i] = mnn[i];
			en[i] = menergy[i];
			for( j=0; j<5; ++j ) {
				aa[i][j] = maa[i][j];				
				}	
			}
		}
/* ----------------------------------------------------
	Transfer Coefficients for Rhodium Anode (itube=2)
	into working arrays
---------------------------------------------------- */
	if( itube==2 ) { /* RHODIUM */
		for( i=0; i<100; ++i ) {
			nn[i] = rnn[i];
			en[i] = renergy[i];
			for( j=0; j<5; ++j ) {
				aa[i][j] = raa[i][j];				
				}	
			}
		}
/* ----------------------------------------------------
	Transfer Coefficients for Tungsten Anode (itube=3)
	into working arrays
---------------------------------------------------- */
	if( itube==3 ) { /* TUNGSTEN */
		for( i=0; i<100; ++i ) {
			nn[i] = tnn[i];
			en[i] = tenergy[i];
			for( j=0; j<5; ++j ) {
				aa[i][j] = taa[i][j];				
				}	
			}
		}
/* ----------------------------------------------------
	If input range is outside limits, reject
---------------------------------------------------- */
	npts = 90;
	for( n=0; n<npts; ++n ) {
		energy = en[n];
		aspec[n] = 0.0;
		if( energy==0.0 || energy > xkv ) goto LOOP;
		arg = xkv;
		dsum = aa[n][0];
		for( j=1; j<4; ++j ) {
			dsum += aa[n][j] * arg;
			arg = arg * xkv;
			}
		if( dsum < 0.0 ) dsum = 0.0;
		aspec[n] = dsum;
LOOP:		;
		}	
	return( 0 );
ERROR_CONDITION:
	printf("Error in MAMSPEC\n");
	return( -1 );
}


