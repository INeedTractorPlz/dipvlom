#!/bin/bash 

for file in *.mp4; do
mv ./${file} ./animations/$1_${file}
done
