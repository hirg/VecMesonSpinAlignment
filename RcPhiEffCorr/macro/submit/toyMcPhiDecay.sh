#!/bin/bash
date

# . ./toyMcPhiDecay.sh
# run10 for 7.7, 11.5, 39 and 62.4 GeV | run11 for 19.6, 27 and 200 GeV | run14 for 14.5

if [ $# -eq 0 ]
then
  PID="Phi"
  Energy="62GeV"
  Name="toyMcPhiDecay"$Energy
  for((counter=11;counter<13;counter=counter+1))
  do
    cp ./run.csh ./$Name${counter}.csh

    echo "mkdir $Name$counter" >> ./$Name${counter}.csh
    echo "cp toyMcPhiDecay.C $Name$counter" >> ./$Name${counter}.csh
    echo "cd $Name$counter" >> ./$Name${counter}.csh

    echo "root4star -b -q -x toyMcPhiDecay.C++'('5','0','1','0','2000000','$counter')'" >> ./$Name${counter}.csh
    echo "cd .." >> ./$Name${counter}.csh
    echo "rm -rf $Name$counter" >> ./$Name${counter}.csh
    qsub -hard -l scratchfree=500,h_cpu=1:30:00,h_vmem=1.8G,projectio=1 -o /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/Log/job$Name$counter.log -e /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/Log/job$Name$counter.err ./$Name$counter.csh
    mv ./$Name${counter}.csh /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/Script/
  done
else
  echo "Wrong number of parameters!!"
fi

