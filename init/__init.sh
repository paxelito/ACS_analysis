#!/bin/sh


rm -r exps.tar.gz
nice python ../initializator.py -o ./exps/ -i 4 -C 2 -H 8 -I ./exps/acsm2s.conf -c0.5 -R 0 -v 1 -a 2 -n2 -s5 -m6 -P0
tar -cvzf exps.tar.gz exps
tar -cvzf carness_package.tar.gz exps.tar.gz luncher_background.sh luncher_screen.sh resetsims.sh

