#!/bin/bash


module load python/3.6.2

python -u process_est.py $1 --temprec $2

./slim -m -s 1234 instance_recipe.slim

