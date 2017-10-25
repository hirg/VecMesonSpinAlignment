#!/bin/bash
date

# . ./submit.sh

if [ $# -eq 0 ]
then
  sbatch --account star --ntasks=1 --array=0-32 --time=08:00:00 FillSpinAlignment.slr
else
  echo "Wrong number of parameters!!"
fi

