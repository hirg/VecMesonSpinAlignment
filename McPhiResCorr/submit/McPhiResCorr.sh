#!/bin/bash
date

# . ./McPhiResCorr.sh

if [ $# -eq 0 ]
then
  Name="McPhi_"
  PID="Phi"
  Energy="200GeV"
  for((counter=0;counter<6;counter=counter+1))
  do
    cp ./run.csh ./run$Name${counter}.csh

    echo "mkdir $Name$counter" >> ./run$Name${counter}.csh
    echo "cp McPhiResCorr.C $Name$counter" >> ./run$Name${counter}.csh
    echo "cd $Name$counter" >> ./run$Name${counter}.csh

    for((rho=0;rho<10;rho=rho+1))
    do
      let Nrho=10*$counter+$rho
      echo "root4star -b -q -x McPhiResCorr.C++'('6','0','0','$Nrho','5000000')'" >> ./run$Name${counter}.csh
    done
    echo "cd .." >> ./run$Name${counter}.csh
    echo "rm -rf $Name$counter" >> ./run$Name${counter}.csh
    qsub -hard -l scratchfree=500,h_cpu=16:00:00,h_vmem=1.8G,projectio=1 -o /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/Log/SpinAlignment/$PID/ResCorr/PYTHIA/job$Name$counter.log -e /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/Log/SpinAlignment/$PID/ResCorr/PYTHIA/job$Name$counter.err ./run$Name$counter.csh
    mv ./run$Name${counter}.csh /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/Script/SpinAlignment/$PID/ResCorr/PYTHIA/
  done
else
  echo "Wrong number of parameters!!"
fi

