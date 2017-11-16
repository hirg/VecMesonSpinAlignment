#!/bin/bash
date

# . ./submit.sh

if [ $# -eq 0 ]
then
  sbatch --account star --ntasks=1 --array=0-19 --time=8:00:00 VecMesonTree.slr # phi TTree mode
  # sbatch --account star --ntasks=1 --array=10000-10019 --time=8:00:00 VecMesonTree.slr # phi TTree resubmit mode
else
  echo "Wrong number of parameters!!"
fi

