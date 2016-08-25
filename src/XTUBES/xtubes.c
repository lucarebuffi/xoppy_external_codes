/* ------------------------------------------------------------
   File Name:   xtubes.c
   Author:      M. Sanchez del Rio based on similar file from
		John M. Boone, Ph.D.
   Date:        Dec-14-1998
   Description: Comand driven interface for MAMSPEC.C.
------------------------------------------------------------ */

#include "stdio.h"
#include "math.h"
#include "stdlib.h"
/* #include "dos.h" */

main()
{
	FILE *pt;
	float xkv,energy[100],spectrum[100];
	int i,ikv,itube;
/* ---------------------------------------------------
	Query for input values
--------------------------------------------------- */
	printf("Enter tube [1=Mo, 2=Rh, 3=W]     :   ");
	scanf("%i",&itube );
	printf("Enter voltage kV in the [18,42] interval :   ");
	scanf("%f",&xkv );

/* ---------------------------------------------------
	Check for errors
--------------------------------------------------- */
 	if (itube >= 4 | xkv <=0 ) {
          printf("Error: itube %i outside the valid interval [1,3]. \n",itube);
        exit(1);
        }
 	if (xkv >= 42 | xkv <=18 ) {
        printf("Error: Energy %2.1f outside the valid interval [18,42]. \n",xkv);
	pt = fopen( "xtubes_tmp.dat","w");
	fprintf(pt,"#F xtubes_tmp.dat \n");
        fprintf(pt,"#S 1 XTUBES result: Error: Voltage %2.1f kV outside the valid interval [18,42]. \n",xkv);
	fprintf(pt,"#C A program descibed in : \n");
        fprintf(pt,"#C J M Boone, R J Jennings, T R Fewell, Med Phys (Dec 1997), 24 (12) 1863-1997 \n");
        fprintf(pt,"#N 2 \n");
        fprintf(pt,"#L Energy [eV] Fluence [photons/sec/mm^2/0.5keV(bw)/mA] \n");
        fprintf(pt,"0 0");
	fclose(pt);
        exit(1);
        }


/* ---------------------------------------------------
	Generate spectrum
--------------------------------------------------- */
	mamspec( itube,xkv,energy,spectrum );
/* ---------------------------------------------------
	Output spectrum to file
--------------------------------------------------- */
	pt = fopen( "xtubes_tmp.dat","w");
	fprintf(pt,"#F xtubes_tmp.dat \n");
	if (itube == 1) 
        fprintf(pt,"#S 1 XTUBES result for Mo x-ray tube at E=%3.1f kV \n", xkv);
	if (itube == 2) 
        fprintf(pt,"#S 1 XTUBES result for Rh x-ray tube at E=%3.1f kV \n", xkv);
	if (itube == 3) 
        fprintf(pt,"#S 1 XTUBES result for W x-ray tube at E=%3.1f kV \n", xkv);
	fprintf(pt,"#C A program descibed in : \n");
        fprintf(pt,"#C J M Boone, R J Jennings, T R Fewell, Med Phys (Dec 1997), 24 (12) 1863-1997 \n");
        fprintf(pt,"#N 2 \n");
        fprintf(pt,"#L Energy [eV]  Fluence [photons/sec/mm^2/0.5keV(bw)/mA] \n");
	for( i=0; i<=89; ++i ) {
		fprintf(pt,"%e  %e\n",energy[i]*1e3,spectrum[i] );
		}
	fclose(pt);
/* ---------------------------------------------------
	Say Goodbye
--------------------------------------------------- */
	printf("\n\nProgram xtubes ended. File xtubes_tmp.dat written to disk.\n\n");
	exit(0);
}
