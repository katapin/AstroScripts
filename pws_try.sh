#!/bin/sh



reb="$6"
for k in `seq $3 $5 $4`
do

name=$$_${k}.ps
nbint=$k

powspec "$1" window="$2" dtnb="2.6" nbint=${nbint} nintfm=10 rebin=${reb} \
outfile=" " plot=yes plotdev=${name}/PS fast=no norm=-2 <<MARK
LAB 2 "Bin time 2.6 x $nbint x $reb"
LAB T $@
lab Y Power, (rms/mean)^2/Hz
R Y 0.01 500
p
date off
p
${name}/PS
p
exit
MARK

eps2png ${name} 

done;
