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
  InPutList="./failed_"$Energy".log"
  grep -l "Critical error" $LogDirectory/*.err > $InPutList

  sed -i 's/err/log/g' $InPutList

  for item in `cat $InPutList`
  do
    cat $item | grep "list" >> $OutPutList
  done
  sed -i 's/Processing VecMesonTree.C("//g' $OutPutList
  sed -i 's/\.list.*/.list/' $OutPutList

  rm $InPutList
else
  echo "Wrong number of parameters"
fi
