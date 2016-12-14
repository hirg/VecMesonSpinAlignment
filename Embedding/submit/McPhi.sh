#!/bin/bash
date

#. ./ResCorr_AMPT.sh

if [ $# -eq 0 ]
then
  PID=Phi
  Energy=200GeV
  # InPutList="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/List/${PID}_list/embedding_List/run_list/embedding.list"
  InPutList="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/List/${PID}_list/embedding_List/all.list"
  OutPutDir="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu200GeV/Embedding/Phi/StMcEvents/StMcEvents_"
  suffix=".root"
  counter=0
  for item in `cat $InPutList`
  do
    cp ./run.csh ./embedding_${PID}_$counter.csh
    OutPutFile=$OutPutDir$counter$suffix

    echo -n "root4star -b -q -x 'run_StMcPhiMaker.C(" >> ./embedding_${PID}_$counter.csh
    echo '"'$item'"'',''"'$OutPutFile'"'')'"'" >> ./embedding_${PID}_$counter.csh
    qsub -hard -l scratchfree=500,h_cpu=00:30:00,h_vmem=1.8G,projectio=1 -o /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/Log/Embedding/$PID/embedding_${PID}_$counter.log -e /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/Log/Embedding/$PID/embedding_${PID}_$counter.err ./embedding_${PID}_$counter.csh

    mv ./embedding_${PID}_$counter.csh /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/Script/Embedding/$PID/
    let counter=$counter+1
  done

else
  echo "Wrong number of parameters"
fi
