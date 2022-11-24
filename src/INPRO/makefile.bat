gfortran -c atomfacm.f
gfortran -c debyenewm.f
gfortran -c fresnelmk.f
gfortran -c iblank.f
gfortran -c iblank2.f
gfortran -c matparmo.f
gfortran -c parametermo.f
gfortran -c perfect_crystalmo.f
gfortran -c strufacnewm.f
gfortran -c inpro.f
gfortran -static-libgfortran -static-libgcc *.o -o inpro.exe
