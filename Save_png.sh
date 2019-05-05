#!/bin/bash 

mkdir ./saved/$1_pictures

for file in *.png; do
mv ./${file} ./saved/$1_pictures/${file}
done
