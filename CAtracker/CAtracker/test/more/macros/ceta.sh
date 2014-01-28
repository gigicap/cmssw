#!/bin/sh

Xvar=0.0001

for I in {0..13}
  do
	echo "eta cut value = $Xvar"
	#cmsRun CAMu10_arg_cfg.py $Xvar
	#mv histotrip.root ecutfiles/Mu10_$I.root
	cmsRun CATTbar_arg_cfg.py $Xvar
	mv histotrip.root ecutfiles/TTbar_$I.root
	#Xvar=$(($Xvar*2)) | bc
	Xvar=$(echo $Xvar*2 | bc)
   done
