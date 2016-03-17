# Bacterial copy number pipeline

## Preprocessing

This steps need to be done once for each new genome

### 1: Create GC profile for genome

```{sh}
GENOMENAME=H37RvCO
mkNucFile.sh ${GENOMENAME}.genome ${GENOMENAME}.fasta
```

Where `${GENOMENAME}.genome` is the `BEDTOOLS` format genome file. The output of this command will be the base profile of the genome. Convert this to the format the R code wants with

```{sh}
Rscript --no-save mkRDA.R \
	${GENOMENAME}_1000by100_TileBin.nuc.gz
```

which will create the file `${GENOMENAME}gcpct.rda`. 

### 2: Make tiling windows bed file

Use bedtools

```{sh}
WINDOW=100
bedtools makewindows \
	-g ${GENOMENAME}.genome \
	-w ${WINDOW} \
	>${GENOMENAME}.genome.Windows_${WINDOW}.bed
```

## Count bins

Make a list of the BAMs to be processed in a file `bamList2` and then count the number of reads per bin per sample with:

```{sh}
bedtools multicov \
	-bed ${GENOMENAME}.genome.Windows_${WINDOW}.bed \
	-bams $(cat bamList2) \
	>counts.txt
```

## Normalize and process with DNAcopy

Use the R script `seqDNAcopyGCCorrection.R` to normalize and then segment the data with `DNAcopy` and plot the segmentation profile and write the segments file. 

Will need to edit the lines at the top of the file depending on your project files. 




