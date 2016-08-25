/* --------------------------------------------------------------------------
File Name:     GENSPEC1.C
Author:        John M. Boone, Ph.D.    (jmboone@ucdavis.edu)
Date:          September 3, 1997

This file contains two subroutines, genspec1.c calls dc_spectral_model.c.
In addition to the normal *.h files, genspec1.h and mual.h are also needed
to compile this code.  Float variables here assume 4 bytes per value (32 bit),
double variables assume 8 bytes (64 bit), and INT variables are also 4 byte
(32 bit) integer variables.
-----------------------------------------------------------------------------
Description: 
============
This routine generates spectra (photons/mm**2 per E) for one mAs and
the output is calibrated to the output of a constant potential generator in 
our laboratory (a Toshiba Model 2050 x-ray generator with a "Rotanode"
x-ray tube with housing model DRX-573HD-S and insert model DR-5735H).
These spectra were calibrated to output (mR/mAs) values and HVL values
that were measured down the central axis of the x-ray beam.  The output
data are given in the manuscript, "An accurate method for computer-
generating tunsten anode x-ray spectra from 30 to 140 kV", JM Boone and 
JA Seibert, Medical Physics 23 or 24 (don't know at this time).
--------------------------------------------------------------------------- */

#include "math.h"
#include "genspec1.h"
#include "mual.h"

/* -----------------------------------------------------------------------------
The following code generates an x-ray spectrum based in the input kV, added 
aluminum filtration, and voltage ripple (in percent). 
Inputs:
   kvolt:  The kilovoltage of the x-ray spectrum to be generated
   almm:   The added aluminum filtration thickness, in millimeters
   ripple: The kilovoltage waveform ripple factor, in percent (0.0 to 100.0)
Output:
    spec[]: An array containing the generated spectrum, with the keV
    corresponding to the index number.  For example, spec[50] contains the 
    the x-ray photon fluence (photons/mm**2) for that spectrum in the energy
    region from 49.5 to 50.5 keV.
----------------------------------------------------------------------------- */

int genspec1( kvolt,almm,ripple,spec )
float kvolt,almm,ripple,spec[];
{
	float xkv,xspec[200];
	float pi,kva,kvb,dkv,kv,x,factor,p;
	int i,j,k,iphase,np;
	if( ripple==0.0 ) {
		dc_spectral_model( kvolt,almm,spec );
		return( 0 );		
		}
	pi = 3.1415926535;
	kva = kvolt;
	dkv = kva * ripple * 0.01;
	kvb = kvolt - dkv;
	for( i=0; i<=150; ++i ) {
		spec[i] = 0.0;
		}
	np = 20;
	p = np;
	for( iphase=1; iphase<=np; ++iphase ) {
		x = (pi/p) * iphase;
		factor = fabs( sin( x ) );
		kv = kvb + dkv * factor;
		dc_spectral_model( kv,almm,xspec );
		for( i=0; i<=150; ++i ) {
			spec[i] += xspec[i];
			}
		}
	for( i=0; i<=150; ++i ) {
		spec[i] = spec[i] / p;
		}
	return( 0 );
}

/* -----------------------------------------------------------------------------
The following code, dc_spectral_model, generates raw spectra with DC waveform
characteristics.  
Inputs:
    kvolt:  The kilovoltage ranging from 30 to 140 kV
    almm:   The added aluminum filtration in mm
Output:
    spec[]: An array containing the generated spectrum, with the keV
    corresponding to the index number.  For example, spec[50] contains the 
    the x-ray photon fluence (photons/mm**2) for that spectrum in the energy
    region from 49.5 to 50.5 keV.
----------------------------------------------------------------------------- */

int dc_spectral_model( kvolt,almm,spec )
float kvolt,almm,spec[];
{
	double arg;
	float mr,xval,x,sum,factor;
	int ien,nterms,nt,iemax,kv,ikv;
	int i,j,k,l,m,n,ians,ival,npts,kpts,nx,ny,kans;
	x = kvolt;
	kv = kvolt;
	for( ien=0; ien<=150; ++ien ) spec[ien] = 0.0;
	for( ien=10; ien<=(kv+3); ++ien ) {
		sum = aa[ien][0];
		arg = x;
		for( j=1; j<nn[ien]; ++j ) {
			sum = sum + aa[ien][j] * arg;
			arg = arg * x;
			}
		if( sum < 0.0 ) sum = 0.0;
		spec[ien] = sum;
		}
	if( almm==0.0 ) return( 0 );
	for( ien=10; ien<150; ++ien ) {
		arg = almm * 0.1 * 2.7 * mual[ien];
		spec[ien] = spec[ien] * exp( -arg );
		}	
	return( 0 );
}

