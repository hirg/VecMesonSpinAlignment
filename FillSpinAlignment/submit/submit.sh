#!/bin/bash
date

# . ./submit.sh

if [ $# -eq 0 ]
then
  # sbatch --account rhstar --ntasks=1 --array=0-9 --time=01:00:00 FillSpinAlignment.slr # test mode
  sbatch --account rhstar --ntasks=1 --array=0-32 --time=08:00:00 FillSpinAlignment.slr # production mode
else
  echo "Wrong number of parameters!!"
fi

