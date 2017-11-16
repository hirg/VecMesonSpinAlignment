#!/bin/bash
date

#. ./VecMesonTree.sh

if [ $# -eq 0 ]
then
  Energy=200GeV
  SM=SE

  OutPutList="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/failed_"$Energy"_"$SM".list"
  rm $OutPutList
  touch $OutPutList

  LogDirectory="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/Log/SLURM"
  InPutList="./completed_"$Energy".log"
  grep -l "Work done" $LogDirectory/*${SM}*.log > $InPutList

  TempList="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/Temp_"$Energy".list"
  for item in `cat $InPutList`
  do
    # cat $item | grep "list" >> $TempList
    cat $item | grep "Processing VecMesonTree.C" >> $TempList
  done
  # sed -i "s/root4star -b -q -x '//g" $TempList
  sed -i 's/Processing VecMesonTree.C("//g' $TempList
  sed -i 's/\.list.*/.list/' $TempList

  CompletedList="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/completed_"$Energy".list"
  cat $TempList | sort > $CompletedList
  rm $TempList

  OriginList="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/origin_"$Energy".list"
  sort ./submit_${Energy}.list > $OriginList

  comm -13 $CompletedList $OriginList > $OutPutList
  rm $OriginList
  rm $CompletedList

  OutPutROOT="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/deleteROOT_"$Energy"_"$SM".list"

  CompletedROOT="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/completedROOT_"$Energy".list"
  TempROOT="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/TempCompletedROOT_"$Energy".list"
  cat $InPutList > $TempROOT
  sed -i 's/Log/SpinAlignment/g' $TempROOT
  sed -i 's/SLURM/Phi\/Forest/g' $TempROOT
  sed -i 's/log/root/g' $TempROOT
  sed -i 's/phi.E/file_/g' $TempROOT
  sed -i 's/GeV_\(.*\)_/GeV_Phi_SE_/g' $TempROOT
  cat $TempROOT | sort > $CompletedROOT
  rm $TempROOT

  OriginROOT="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/originROOT_"$Energy".list"
  OriginForest="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/SpinAlignment/Phi/Forest"
  ls -d $OriginForest/*${SM}*.root | sort > $OriginROOT

  comm -13 $CompletedROOT $OriginROOT > $OutPutROOT
  # # sed -i 's/Script/SpinAlignment\/Resolution/g' $OutPutROOT # resolution mode
  # # sed -i "s/runPhiSE${Energy}/file_${Energy}_Resolution_/g" $OutPutROOT
  # sed -i 's/Script/SpinAlignment\/Phi\/Forest/g' $OutPutROOT # phi-meson production mode
  # sed -i "s/runPhi${SM}${Energy}/file_${Energy}_Phi_${SM}_/g" $OutPutROOT
  # sed -i "s/csh/root/g" $OutPutROOT

  rm $InPutList
  rm $OriginROOT
  rm $CompletedROOT
else
  echo "Wrong number of parameters"
fi
