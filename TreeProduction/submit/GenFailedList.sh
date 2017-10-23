#!/bin/bash
date

#. ./VecMesonTree.sh

if [ $# -eq 0 ]
then
  Energy=200GeV
  SM=SE

  OutPutList="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/failed_"$Energy".list"
  rm $OutPutList
  touch $OutPutList

  LogDirectory="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/Log"
  InPutList="./completed_"$Energy".log"
  grep -l "Work done" $LogDirectory/*.log > $InPutList

  sed -i 's/Log/Script/g' $InPutList
  sed -i 's/job/run/g' $InPutList
  sed -i 's/log/csh/g' $InPutList

  TempList="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/Temp_"$Energy".list"
  for item in `cat $InPutList`
  do
    cat $item | grep "list" >> $TempList
  done
  sed -i "s/root4star -b -q -x '//g" $TempList
  sed -i 's/VecMesonTree.C("//g' $TempList
  sed -i 's/\.list.*/.list/' $TempList

  CompletedList="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/completed_"$Energy".list"
  cat $TempList | sort > $CompletedList
  rm $TempList

  OriginList="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/origin_"$Energy".list"
  sort ./submit_${Energy}.list > $OriginList

  comm -13 $CompletedList $OriginList > $OutPutList
  rm $OriginList
  rm $CompletedList

  OutPutROOT="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/deleteROOT_"$Energy".list"
  OriginROOT="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/originROOT_"$Energy".list"
  OriginScript="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/Script"
  ls -d $OriginScript/*.csh | sort > $OriginROOT

  comm -13 $InPutList $OriginROOT > $OutPutROOT
  # sed -i 's/Script/SpinAlignment\/Resolution/g' $OutPutROOT # resolution mode
  # sed -i "s/runPhiSE${Energy}/file_${Energy}_Resolution_/g" $OutPutROOT
  sed -i 's/Script/SpinAlignment\/Phi\/Forest/g' $OutPutROOT # phi-meson production mode
  sed -i "s/runPhi${SM}${Energy}/file_${Energy}_Phi_${SM}_/g" $OutPutROOT
  sed -i "s/csh/root/g" $OutPutROOT

  rm $InPutList
  rm $OriginROOT
else
  echo "Wrong number of parameters"
fi
