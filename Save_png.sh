#!/bin/bash 

for file in *.png; do
mv ./${file} ./animations/$1_${file}
done
