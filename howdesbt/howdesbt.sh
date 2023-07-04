#!/bin/bash
path="../data/toy"
file="toy"
k=9
outfile="${file}-k${k}"
i=0
input=""

## Make BF
mkdir $outfile
mkdir $outfile/bf
while IFS= read -r line || [ -n "$line" ]; do
    input="${input}${line}"
    if (($((i%2))==1))
    then
        # echo $input
        howdesbt makebf /dev/fd/0 --out="./${outfile}/bf/$((i/2)).bf" --k=${k}
    else
        input=${line}
    fi
    ((i++))
done < "${path}/${file}.fasta"

cd ${outfile}

## Cluster
ls ./bf/*.bf > leafnames
howdesbt cluster --list=leafnames --bits=80K \
  --tree=union.sbt --nodename=node{number} --keepallnodes
rm leafnames

## Build Tree Topology
howdesbt build --HowDe --tree=union.sbt --outtree=howde.sbt
# mkdir rrr
# mv *.rrr.bf ./rrr/

# ## Query
# howdesbt query --tree=howde.sbt "../${path}/${file}.fasta" > ./${outfile}.dat


