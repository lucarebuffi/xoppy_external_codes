#
# compile in linux all applications used in xoppy EXCEPT xpowder_fml
#

rm ../../../bin/linux/*

cd DEJUS/TC
make -f makefile.unix clean
make -f makefile.unix 
sleep 0.5
cp linux/tc ../../../bin/linux/tc
sleep 0.5
make -f makefile.unix clean
cd ../../

cd DEJUS/US
make -f makefile.unix clean
make -f makefile.unix 
sleep 0.5
cp linux/us ../../../bin/linux/us
sleep 0.5
make -f makefile.unix clean
cd ../../

cd DEJUS/WS
make -f makefile.unix clean
make -f makefile.unix 
sleep 0.5
cp linux/ws ../../../bin/linux/ws
sleep 0.5
make -f makefile.unix clean
cd ../../

cd MLAYER
make -f makefile.unix clean
make -f makefile.unix 
sleep 0.5
cp mlayer ../../bin/linux/mlayer
sleep 0.5
make -f makefile.unix clean
cd ..

cd URGENT
make -f makefile.unix clean
make -f makefile.unix 
sleep 0.5
cp ./linux/urgent ../../bin/linux/urgent
sleep 0.5
make -f makefile.unix clean
cd ..

cd XTUBES
make -f Makefile clean
make -f Makefile
sleep 0.5
cp xtubes ../../bin/linux/xtubes
sleep 0.5
make -f Makefile clean
cd ..

cd XTUBE_W
make -f Makefile clean
make -f Makefile
sleep 0.5
cp tasmip ../../bin/linux/tasmip
sleep 0.5
make -f Makefile clean
cd ..

cd INPRO
make -f makefile.unix clean
make -f makefile.unix 
sleep 0.5
cp inpro ../../bin/linux/inpro
sleep 0.5
make -f makefile.unix clean
cd ..

cd XCOM31
make -f makefile.unix clean
make -f makefile.unix 
sleep 0.5
cp xcom ../../bin/linux/xcom
sleep 0.5
make -f makefile.unix clean
cd ..

rm -rf crystal
git clone https://github.com/srio/crystal
cd crystal
make
sleep 0.5
cp diff_pat ../../bin/linux/diff_pat
cd ..
rm -rf crystal




