@echo off
rem ****
rem ****---- Compilation for Calculation_of_Powder_Patterns Program ----****
rem ****
rem **** Author: JRC
rem **** Revision: June-2009
rem ****

:G95
   g95 -c -O3  -std=f2003  -funroll-loops  -msse2   xpowder_fml.f90   -IC:\CrysFML\G95\LibC
   g95  *.o -o xpowder_fml -O3  -funroll-loops  -msse2  -LC:\CrysFML\G95\LibC -lcrysfml  -Wl,--heap=0x01000000
   goto END
rem

:END
   del *.obj *.mod *.o *.map *.bak > nul
