#!/bin/bash 

mkdir ./saved/$1_data

for file in *.dat; do
cp ./${file} ./saved/$1_data/${file}
done
