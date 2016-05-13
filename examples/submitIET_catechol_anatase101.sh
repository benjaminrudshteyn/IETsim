#!/bin/bash
#BSUB -o .
#BSUB -J catechol_anatase101
#BSUB -q shared
#BSUB -W 06:00
#BSUB -n 1
#BSUB -M 5000
#BSUB -R "rusage[mem=5000]" 

#This script requests one processor on the shared queue for 6 hours using 5000 MB on each processor (the total as well)

echo "Slot count: $LSB_DJOB_NUMPROC"

~/dynamics-0.4.2/dynamics2/dynamics catechol_anatase101.bind

