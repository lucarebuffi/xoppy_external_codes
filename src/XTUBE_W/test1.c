/* ------------------------------------------------------------
   File Name:   test1.c
   Author:      John M. Boone, Ph.D.
   Modified:    M. Sanchez del Rio 98-12-16 to be called by XOP
   Date:        Sep-3-1997
   Description: Test output for GENSPEC1.C, the TASMIP spectrum code
------------------------------------------------------------ */

#include "stdio.h"
#include "math.h"
#include "stdlib.h"
/* #include "dos.h" */

main()
{
	FILE *pt;
	float kv,almm,ripple,spec[200];
	int i,ikv;
/* ---------------------------------------------------
	Query for input values
--------------------------------------------------- */
	printf("Enter voltage kV in the interval [30,140] :   ");
	scanf("%f",&kv );
	printf("Enter added Al (mm): ");
	scanf("%f",&almm );
	printf("Enter ripple (%%)   : ");
	scanf("%f",&ripple );


/* ---------------------------------------------------
	Check for errors
--------------------------------------------------- */
 	if (kv > 140 | kv <30 ) {
        printf("Error: Energy %3.1f outside the valid interval [30,140]. \n",kv);
	pt = fopen( "tasmip_tmp.dat","w");
	fprintf(pt,"#F tasmip_tmp.dat \n");
        fprintf(pt,"#S 1 TASMIP result: Error: Voltage %3.1f kV outside the valid interval [30,140]. \n",kv);
	fprintf(pt,"#C A program descibed in : \n");
        fprintf(pt,"#C J M Boone & J A Seibert, Med Phys (Nov 1997), 24 (11) 1661-1670 \n");
        fprintf(pt,"#N 2 \n");
        fprintf(pt,"#L Energy [eV]  Flux \n");
        fprintf(pt,"0 0");
	fclose(pt);
        exit(1);
        }



/* ---------------------------------------------------
	Generate spectrum
--------------------------------------------------- */
	genspec1( kv,almm,ripple,spec );
/* ---------------------------------------------------
	Output spectrum to file
--------------------------------------------------- */
	pt = fopen( "tasmip_tmp.dat","w");
	fprintf(pt,"#F tasmip_tmp.dat \n");
        fprintf(pt,"#S 1 TASMIP result for W x-ray tube at E=%3.1f kV (%2.1f mm of Al filter and %2.0f%% ripple) \n", kv, almm, ripple);
	fprintf(pt,"#C A program descibed in : \n");
        fprintf(pt,"#C J M Boone & J A Seibert, Med Phys (Nov 1997), 24 (11) 1661-1670 \n");
        fprintf(pt,"#C Outputs: Energy [eV] and Flux [photons/1keV(bw)/mA/mm^2(@1m)/sec] \n");
        fprintf(pt,"#N 2 \n");
        fprintf(pt,"#L Energy [eV]  Flux [photons/1keV(bw)/mA/mm^2(@1m)/sec] \n");

	ikv = kv + 3.0;
	for( i=0; i<=(ikv); ++i ) {
		fprintf(pt,"%3.1f  %e\n",1e3*i,spec[i] );
		}
	fclose(pt);
/* ---------------------------------------------------
	Say Goodbye
--------------------------------------------------- */
	printf("\n\nProgram tasmip end. File tasmip_tmp.dat written to disk.\n\n");
	return 0;
}
