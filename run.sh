#!/bin/bash
for i in {1..32};
do
    echo Sample$i
    python subtract_background.py -b range_maps/background/RIS40x40_OSF4x4_RDMest.mat -r range_maps/sample$i/RIS40x40_OSF4x4_RDMest.mat
done