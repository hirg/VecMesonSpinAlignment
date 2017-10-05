#!/bin/bash
date

#. ./rVecMesonTree.sh

if [ $# -eq 0 ]
  then
    PID=Phi
    Energy=200GeV
    SM=SE
    InPutList="./failed_"$Energy".list"
    counter=10000;
    for item in `cat $InPutList`
    do
      Name=$PID$SM$Energy
      cp ./run.csh ./run$Name$counter.csh

      echo -n "root4star -b -q -x 'VecMesonTree.C(" >> run$Name$counter.csh
      echo -n '"'$item'",' >> run$Name$counter.csh
      echo -n $counter',' >> run$Name$counter.csh
###############################mode###################################
     # echo -n 0',' >> run$Name$counter.csh  # fill ReCenterPar mode
     # echo -n 1',' >> run$Name$counter.csh  # fill ShiftPar mode
     echo -n 2',' >> run$Name$counter.csh  # Resolution mode
     # echo -n 3',' >> run$Name$counter.csh  # Phi mode
###############################mode###################################

#############################energy###################################
      # echo -n 0',' >> run$Name$counter.csh  # 7.7GeV
     # echo -n 1',' >> run$Name$counter.csh  # 11.5GeV
     # echo -n 2',' >> run$Name$counter.csh  # 19.6GeV
     # echo -n 3',' >> run$Name$counter.csh  # 27GeV
     # echo -n 4',' >> run$Name$counter.csh  # 39GeV
     # echo -n 5',' >> run$Name$counter.csh  # 62.4GeV
     echo -n 6',' >> run$Name$counter.csh  # 200GeV
###############################energy###################################

###############################flag_ME###################################
      echo 0')'"'" >> run$Name$counter.csh  # Same Event 
     # echo 1')'"'" >> run$Name$counter.csh  # Mixed Event
###############################flag_ME###################################

      # qsub -hard -l projectio=1,scratchfree=500,h_cpu=8:00:00,h_vmem=1.8G -o /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Log/job$Name$counter.log -e /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Log/job$Name$counter.err ./run$Name$counter.csh

      mv run$Name$counter.csh /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Script/
      let counter=counter+1;
    done

  else
    echo "Wrong number of parameters"
fi
