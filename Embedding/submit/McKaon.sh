#!/bin/bash
date

#. ./McKaon.sh

if [ $# -eq 0 ]
then
  PID=Kplus
  Energy=39GeV
  # InPutList="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/List/Kaon_list/run11/${PID}.list"
  InPutList="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/List/Kaon_list/${PID}.list"
  OutPutDir="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/Embedding/$PID/StMcEvents/StMcEvents_"
  suffix=".root"
  counter=0
  for item in `cat $InPutList`
  do
    cp ./run.csh ./embedding_${PID}_$counter.csh
    OutPutFile=$OutPutDir$counter$suffix

    echo -n "root4star -b -q -x 'run_StMcAnalysisMaker.C(" >> ./embedding_${PID}_$counter.csh
    echo '"'$item'"'',''"'$OutPutFile'"'')'"'" >> ./embedding_${PID}_$counter.csh
    qsub -hard -l scratchfree=500,h_cpu=00:30:00,h_vmem=1.8G,projectio=1 -o /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/Log/Embedding/Kaon/embedding_${PID}_$counter.log -e /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/Log/Embedding/Kaon/embedding_${PID}_$counter.err ./embedding_${PID}_$counter.csh

    mv ./embedding_${PID}_$counter.csh /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/Script/Embedding/Kaon/
    let counter=$counter+1
  done

else
  echo "Wrong number of parameters"
fi
