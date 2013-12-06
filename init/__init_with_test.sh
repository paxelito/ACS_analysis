#!/bin/sh


rm -r exps.tar.gz test
cd ./exps/

nice python ../initializator.py -i 4 -C 2 -H 8 -I ./acsm2s.conf -c0.5 -R 0 -v 1 -a 2 -n2 -s5 -m6 -P2
cd ..
tar -cvzf exps.tar.gz exps
tar -cvzf carness_package.tar.gz exps.tar.gz luncher_background.sh luncher_screen.sh resetsims.sh

rm -r test
mkdir test
mv carness_package.tar.gz test
cd test
tar -xvzf carness_package.tar.gz 
