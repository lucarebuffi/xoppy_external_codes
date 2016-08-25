/* ------------------------------------------------------------
   File Name:   mainprog.c  (new file)
   Author:      John M. Boone, Ph.D.
   Date:        12-15-97  17:15
   Description: This is a main program which calls MAMSPEC
------------------------------------------------------------ */

#include "stdio.h"
#include "math.h"
#include "stdlib.h"
/* #include "dos.h" */

main()
{
	FILE *pt;
	float energy[100],spectrum[100],xkv;
	int ikv,itube,i;
	char fname[100];
	for( itube=1; itube<=3; ++itube ) {
		printf("===== tube: %d ======\n",itube );
		for( ikv=20; ikv<=40; ikv+=5 ) {
			xkv = ikv;
			printf("kV: %2d ....",ikv);
			if( itube==1 ) sprintf( fname,"moly_%d.aaa",ikv );
			if( itube==2 ) sprintf( fname,"rhod_%d.aaa",ikv );
			if( itube==3 ) sprintf( fname,"tung_%d.aaa",ikv );
			printf("Calc...");
			mamspec( itube,xkv,energy,spectrum );
			printf("Store...");
			pt = fopen( fname,"w" );
			for( i=0; i<=90; ++i ) {
				fprintf(pt,"%6.3f  %e\n",energy[i],spectrum[i] );
				}
			fclose(pt);
			printf("...done\n");
			}
		}
EXIT:	printf("\n\n[END]\n\n");
	exit(0);
}
