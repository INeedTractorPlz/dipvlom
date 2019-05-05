#!/bin/bash 

mkdir ./saved/$1_animations

for file in *.mp4; do
mv ./${file} ./saved/$1_animations/${file}
done

