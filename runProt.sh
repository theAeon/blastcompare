#!/bin/bash
echo BEGIN CHICKEN
python3 ../get_diff.py --compare-protein 1 -a -p chickenvfish
python3 ../get_diff.py --compare-protein 1 -a -p chickenvfrog
python3 ../get_diff.py --compare-protein 1 -a -p chickenvhuman
echo BEGIN FROG
python3 ../get_diff.py --compare-protein 1 -a -p frogvfish
python3 ../get_diff.py --compare-protein 1 -a -p frogvchicken
python3 ../get_diff.py --compare-protein 1 -a -p frogvhuman
echo BEGIN FISH
python3 ../get_diff.py --compare-protein 1 -a -p fishvfrog
python3 ../get_diff.py --compare-protein 1 -a -p fishvchicken
python3 ../get_diff.py --compare-protein 1 -a -p fishvhuman
echo BEGIN HUMAN
python3 ../get_diff.py --compare-protein 1 -a -p humanvfrog
python3 ../get_diff.py --compare-protein 1 -a -p humanvfish
python3 ../get_diff.py --compare-protein 1 -a -p humanvchicken
