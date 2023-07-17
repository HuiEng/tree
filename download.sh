#!/bin/bash
mkdir data/$1
cd data/$1
folder=/home/lawhe/benchmarking-data/$1
cmd="ls ${folder}/parktree/*.txt ${folder}/sigClust/*.txt ${folder}/tree/*/prim*.txt ${folder}/tree/*/sec*.txt"
ssh lawhe@sef-bigdata-08 "${cmd}" >list
scp lawhe@sef-bigdata-08:{"$(cat list | tr '\n' ',' | sed 's/,$//')"} .
