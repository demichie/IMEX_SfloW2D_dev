#!/bin/bash
cd /home/user_sw/SW_RUNS

case $1 in
   plot_overlay)
      echo "plot_overlay"
      touch matplotlibrc
      echo "backend : agg" >> matplotlibrc      
      python3 /home/user_sw/IMEX_SfloW2D-master/UTILS/plot_overlay.py
      rm matplotlibrc;;
   run)   
      echo "run"
      /home/user_sw/IMEX_SfloW2D-master/bin/IMEX_SfloW2D;;
   *)
      echo "Please select: run, plot_overlay"
      echo"";;
esac

#
# pwd
