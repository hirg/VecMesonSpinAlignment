#!/bin/bash
date

#. ./FillSpinAlignment.sh

if [ $# -eq 0 ]
  then
#    counter=0
    PID=Phi
    List_SM=SE
    SM=_${List_SM}_
    Energy=200GeV
    Name="_${Energy}_$PID$SM"
    suffix=".root"
    # for((counter=0;counter<3;counter=counter+1)) # 19GeV
    # for((counter=0;counter<5;counter=counter+1)) # 27GeV
    # for((counter=0;counter<7;counter=counter+1)) # 62GeV
    # for((counter=0;counter<162;counter=counter+1)) # 200GeV w/o ToF
    for((counter=0;counter<33;counter=counter+1)) # 200GeV with ToF
    do
      cp ./run.csh ./run$Name$counter.csh

      echo -n "root4star -b -q -x 'FillSpinAlignmentZdcSmd.C(" >> run$Name$counter.csh
###############################energy###################################
      # echo -n 0',' >> run$Name$counter.csh  # 7.7GeV
      # echo -n 1',' >> run$Name$counter.csh  # 11.5GeV
      # echo -n 2',' >> run$Name$counter.csh  # 19.6GeV
      # echo -n 3',' >> run$Name$counter.csh  # 27GeV
      # echo -n 4',' >> run$Name$counter.csh  # 39GeV
      # echo -n 5',' >> run$Name$counter.csh  # 62.4GeV
      echo -n 6',' >> run$Name$counter.csh  # 200GeV
###############################energy###################################

###############################X_flag###################################
      echo -n 0',' >> run$Name$counter.csh # Same Event 
      # echo -n 1',' >> run$Name$counter.csh  # Mixed Event
###############################X_flag###################################

      echo -n $counter',' >> run$Name$counter.csh # List

############################start_event#################################
      echo -n 0',' >> run$Name$counter.csh  # start_event
############################start_event#################################

#############################stop_event#################################
      echo -n 1000000000',' >> run$Name$counter.csh  # stop_event
      # echo -n 100024',' >> run$Name$counter.csh  # stop_event: test mode
#############################stop_event#################################

##############################Partilce##################################
      echo 0')'"'" >> run$Name$counter.csh  # phi-meson
      # echo 1')'"'" >> run$Name$counter.csh  # K*
      # echo 2')'"'" >> run$Name$counter.csh  # K0S
##############################Partilce##################################

      qsub -hard -l scratchfree=500,h_cpu=8:00:00,h_vmem=1.8G,projectio=1 -o /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/Log/job$Name$counter.log -e /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/Log/job$Name$counter.err ./run$Name$counter.csh

      mv run$Name$counter.csh /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/Script/
    done

  else
    echo "Wrong number of parameters"
fi
