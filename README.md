###################################################
#						  #
#  Describes how to prepare gene MSA for Figures  #
#                                                 #
###################################################

```
gene=MLXIPL
dic=mammals_birds.dict.csv
order=mammals_birds.phylo_ordered.lst


nectar=HLlepAsp1,HLcinPul1,HLlicMelCas1,HLlicPen1,HLphyNov1,HLfloFus1,HLphaSup1,HLcalAnn5,HLparPun1,HLdicExi1,HLamaAes1,HLaraSol1,HLtriMol2,HLlorGal1

gene=HK3
gene=SLC2A5
```

## Get prot alignments of full gene and BEB.nogaped
```
scp eosipova@10.10.50.202:/projects/project-eosipova/NectarivoryProject/$gene/prot* .
```

## Combine these files together and align
```
#cat prot.nectar.*nogaps.fa <(get_seq_by_name_fasta.py -ab -n $(echo $nectar) -f prot.*.$gene.fa) > combined.prot.$gene.fa
```
# better add only one nogaps nectar seq as reference and align 
```
mafft combined.prot.$gene.fa > mafft.combined.prot.$gene.fa
```

## Replace assembly codes wiht common names and sort MSA by phylogeny
```
dict=/Users/osipova/Documents/LabDocs/Birds_Phylogeny/45_birds.dict.csv
order=/Users/osipova/Documents/LabDocs/Birds_Phylogeny/45_birds.phylo_ordered.lst

rename_fasta_from_dict.py -f mafft.combined.prot.$gene.fa -d $dict > commonNames.combined.prot.$gene.fa

sort_fasta_by_list.py -f commonNames.combined.prot.$gene.fa -l $order -r chicken > ordered.commonNames.combined.prot.$gene.fa
```


###################################################
#                                                 #
#    Describes how to prepare expression plots	  #
#                                                 #
###################################################

```
RNADIR=/Users/osipova/Documents/LabDocs/Bird_transcriptomics/Transcriptome_for_NectarGenomics/DESeq2_results_Kallisto/

sp1=Annas_hummingbird; sp2=common_swift
sp1=New_Holland_honeyeater; sp2=zebra_finch
sp1=rainbow_lorikeet; sp2=cockatiel

# do for all 3 pairs of species
for gene in ANK2 ANO1 JADE1 FANCM FAT1 FAT4 GRIP1 GRIP2 SNX19 HK1 HK2 HK3 GCK HKDC1 GRIK1 GRIK3 GRIK4 HYAL1 HYAL2; \
do \
   if [ ! -d $gene ]; then
   echo $gene; \
   md $gene; \
   for norm in $(ls $RNADIR/norm_counts.deseq2_res.*.tsv); \
   do \
     paste <(grep -w "^$gene" $norm | tr '\t' '\n' | grep -v $gene) <(for l in $(head -1 $norm); do grep -w $l $RNADIR/../ALL_info.tsv |cut -f3,4; done); \
   done >> $gene/norm_counts.tsv; \
   awk 'NF==3{print}' $gene/norm_counts.tsv > file; \
   cat <(echo -e "count\tbird\ttissue") file > $gene/norm_counts.tsv; \
   fi;
done
```


### Check expression of target ChREBP genes!
gene lists are from: Iizuka_2017, https://doi.org/10.1016/j.bbadis.2016.11.029

ACACA has many2many orthology in all species
try to quantify one (highly-expressed) isoform of ACACA (or low-expressedd: rna-XM_030288106): 

```
for l in $(cat ALL_info.tsv | cut -f1); \
do \
  count=$(grep "ACACA_rna-XM_015295696\|rna-XM_030288108" Kallisto_quant_ALL/$l/abundance.tsv | cut -f4); \
  info=$(grep $l ALL_info.tsv | cut -f3,4); \
  echo -e "$count\t$info"; \
done | awk 'NF==3{print}' > ~/Documents/LabDocs/NectarivoryProject/absrel/PAML_BEB/ACACA_rna-XM_015295696/norm_counts.tsv
```

### Recover SI expression data for honeyeater and zebra finch
```
g=SI

db=phyNov
t=SI_rna-XM_015291762.2.410

db=taeGut
t=SI_rna-XM_002194436.3

# for both run:
for l in $(grep $db $INFO | cut -f1); do count=$(grep $t $RNADIR/Kallisto_quant_ALL/$l/validated.abundance.tsv | cut -f4 ); info=$(grep $l $INFO | cut -f3,4);  echo -e "$count\t$info"; done >> $g/norm_counts.tsv
```


### Recover AMY1A expression data for lorikeet and cockatiel
```
g=AMY1A

db=triMol
t=AMY1A_rna-XM_002186591.3.49

db=nymHol
t=AMY1A_mRNA16263

# for both run:
for l in $(grep $db $INFO | cut -f1); do count=$(grep $t $RNADIR/Kallisto_quant_ALL/$l/validated.abundance.tsv | cut -f4 ); info=$(grep $l $INFO | cut -f3,4);  echo -e "$count\t$info"; done >> $g/norm_counts.tsv
```


### Recover AMY2A expression data for lorikeet and cockatiel
```
g=AMY2A

for db in calAnn apuApu triMol nymHol phyNov taeGut; \
do \
  for l in $(grep $db $INFO | cut -f1); \
  do \
    count=$(grep $g $RNADIR/Kallisto_quant_ALL/$l/validated.abundance.tsv | cut -f4 | sum_stdin.py); \
    info=$(grep $l $INFO | cut -f3,4); \
    echo -e "$count\t$info"; \
  done >> $g/norm_counts.tsv;
done
```



