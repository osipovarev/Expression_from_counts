#!/usr/bin/env/ bash

## dir with RNAseq results
RNADIR=$1

## file with list of genes
GENELIST=$2


while IFS= read -r gene;
do
    if [ ! -d "$gene" ];
    then
        mkdir -p $gene;
        for norm in $(ls $RNADIR/norm_counts.deseq2_res.*.tsv); 
        do 
            paste <(grep -w "^$gene" $norm | tr '\t' '\n' | grep -v $gene) <(for l in $(head -1 $norm); do grep -w $l $RNADIR/../ALL_info.tsv |cut -f3,4; done); 
        done >> $gene/norm_counts.tsv; 
        awk 'NF==3{print}' $gene/norm_counts.tsv > file; 
        cat <(echo -e "count\tbird\ttissue") file > $gene/norm_counts.tsv;
    fi 
done < $GENELIST
