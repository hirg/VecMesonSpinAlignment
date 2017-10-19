#!/bin/bash
date

# . ./submit.sh

if [ $# -eq 0 ]
then
  sbatch --account star --ntasks=1 --array=0-806 --time=03:00:00 VecMesonTree.slr
else
  echo "Wrong number of parameters!!"
fi

