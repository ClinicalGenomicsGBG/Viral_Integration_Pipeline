# Virus integration Pipeline (SE)

## Background

This pipeline uses mapped reads to look for human (hg19) viral integrations directly from the alignment by extracting the SA tag. It also bin the breakpints by taking the highest amount of reads and select 10nt on either side. The breakpoints are annotated using annovar. 

Finally an additional PCR contamination step is added to the nanopore pipelines which filter out breakpoints found over multiple barcodes. The barcode with the highest amount of that fusion is kept while if detected in other barcodes removed (the bin is also 10 for this filtering step). 

## For NanoPore

The reads are aligned using minimap2. We extract the breakpoints using the screening for integration scripts that looks for the S flag. The reads are binned with 10 nt and annotated with annovar. The output from the annotation is filtered from PCR contamination.


### How to run 

The reads are mapped using the ```runMinimap2.sh``` script

```

qsub runMinimap2.sh sample1.fastq

```

The viral integrations are extracted using ```ScreeningFusionsForSE_Pipeline_vs4.py```

(the genome is the actual viral genome, if you mapped to another viral genome then HBV_D just make sure you include it in the reference hg19 and change the -g for the new viral genome)

```

python ScreeningFusionsForSE_Pipeline_vs4.py -b sample1.bam -g HBV_D

samtools view sample1.bam -h | grep -e $'\tSA:' -e ^@ | grep -e HBV_D -e ^@ | samtools view -Sb - > sample1_Chimeric.bam

```

Then we bin and annotate using annovar in script ```BinAndAnnotateVirusSEpipelien_10bin.py``` *Obs* I don't run this step in the nodes. It looks like it is overwriting the output in the annotation step so run it one by one in the login node! 


```

python BinAndAnnotateVirusSEpipelien_10bin.py sample1.OutSE_FusionPipe.txt HBV_D sample1.OutFusionPipe_Binned_10nt_Annotated.txt

```

Finally we remove the PCR contaminations by looking at fusions in the same bin shared across the barcodes. This is done with the script ```FilterViralOutForPCRSpilling_2.py``` If it is we save the one that is most common.

```

python FilterViralOutForPCRSpilling_2.py -i *.OutFusionPipe_Binned_10nt_Annotated.txt

```


## For short read sequencing (iontorrent SE)

The pipeline for the iontorrent is pretty much the same as the Nanopore. The only difference is that we are removing the common adapter GCCAGGTTCCAGTCAC and that we are aligning the reads with BWA. 

### How to run


The reads are filtered from the common tag GCCAGGTTCCAGTCAC


```

qsub Cutadapt.sh sample1.fastq.gz

```


The adapter filtered reads are mapped using the ```runBWA_HBV.sh``` script

```

qsub runBWA_HBV.sh sample1_Adapterfilt.fastq.gz

```

The viral integrations are extracted using ```ScreeningFusionsForSE_Pipeline_vs4.py```

(the genome is the actual viral genome, if you mapped to another viral genome then HBV_D just make sure you include it in the reference hg19 and change the -g for the new viral genome)

```

python ScreeningFusionsForSE_Pipeline_vs4.py -b sample1.bam -g HBV_D

samtools view sample1.bam -h | grep -e $'\tSA:' -e ^@ | grep -e HBV_D -e ^@ | samtools view -Sb - > sample1_Chimeric.bam

```

Then we bin and annotate using annovar in script ```BinAndAnnotateVirusSEpipelien_10bin.py``` *Obs* I don't run this step in the nodes. It looks like it is overwriting the output in the annotation step so run it one by one in the login node! 

```

python BinAndAnnotateVirusSEpipelien_10bin.py sample1.OutSE_FusionPipe.txt HBV_D sample1.OutFusionPipe_Binned_10nt_Annotated.txt

```

Finally we remove the PCR contaminations by looking at fusions in the same bins shared across the barcodes. This is done with the script ```FilterViralOutForPCRSpilling_2.py``` If it is we save the one that is most common.

```

python FilterViralOutForPCRSpilling_2.py -i *.OutFusionPipe_Binned_10nt_Annotated.txt

```

