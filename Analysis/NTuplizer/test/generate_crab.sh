#!/bin/bash

config=JEC_MC_MINIAOD_cfg.py
publish=False
site=T2_IN_TIFR

fil_list=submit.sh
mon_list=monitor.sh


function loopOverArray(){

    jq -c '[].specialId' samples.json | while read i; do
        # Do stuff here
        echo "$i"
    done
}

loopOverArray
