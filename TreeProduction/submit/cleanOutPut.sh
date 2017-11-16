#!/bin/bash
date

#. ./cleanOutPut.sh

if [ $# -eq 0 ]
  then
    Energy=200GeV
    SM=SE
    InPutList="./deleteROOT_${Energy}_${SM}.list"
    for item in `cat $InPutList`
    do
      rm $item
    done
    rm $InPutList

  else
    echo "Wrong number of parameters"
fi
