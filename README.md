# Virus integration Pipeline (SE)

## Background

This pipeline uses mapped reads to look for human (hg19) viral integrations directly from the alignment by extracting the SA tag. It also bin the breakpints by taking the highest amount of reads and select 10nt on either side. The breakpoints are annotated using annovar. 

Finally an additional PCR contamination step is added to the nanopore pipelines which filter out breakpoints found over multiple barcodes. The barcode with the highest amount of that fusion is kept while if detected in other barcodes removed (the bin is also 10 for this filtering step). 

## For NanoPore

The reads are aligned using minimap2. We extract the breakpoints using the screening for integration scripts that looks for the S flag. The reads are binned with 10 nt and annotated with annovar. The output from the annotation is filtered from PCR contamination.


### How to run 

The reads are mapped using the ```runMinimap2.sh``` script

#### Alignment

```

qsub runMinimap2.sh sample1.fastq

```

#### Extracting Integrations

The viral integrations are extracted using ```ScreeningFusionsForSE_Pipeline_vs4.py```

(the genome is the actual viral genome, if you mapped to another viral genome then HBV_D just make sure you include it in the reference hg19 and change the -g for the new viral genome)

```

python ScreeningFusionsForSE_Pipeline_vs4.py -b sample1.bam -g HBV_D

samtools view sample1.bam -h | grep -e $'\tSA:' -e ^@ | grep -e HBV_D -e ^@ | samtools view -Sb - > sample1_Chimeric.bam

```

#### Within Sample binning

Then we bin and annotate using annovar in script ```BinAndAnnotateVirusSEpipelien_10bin.py``` *Obs* I don't run this step in the nodes. It looks like it is overwriting the output in the annotation step so run it one by one in the login node! 


```

python BinAndAnnotateVirusSEpipelien_10bin.py sample1.OutSE_FusionPipe.txt HBV_D sample1.OutFusionPipe_Binned_10nt_Annotated.txt

```

#### PCR contamination removal (Without or With SCOPE)

Finally we remove the PCR contaminations by looking at fusions in the same bin shared across the barcodes. This is done with the script ```FilterViralOutForPCRSpilling_2.py``` If it is we save the one that is most common.

```

python FilterViralOutForPCRSpilling_2.py -i *.OutFusionPipe_Binned_10nt_Annotated.txt

```

OBS, if the samples are SCOPE samples we cannot use the regular PCR filtering. This is due to multiple samples can be from different barcodes so we should take the barcode group into consiteration. This is done by script ```FilterViralOutForPCRSpilling_SCOPE.py```. Input to this is a metadata file, commaseperated containing the sample name and the group and a sample name for the output. Example of metadataformat below (run it in the same folder as where you have the output from BinAndAnnotateVirusSEpipelien_10bin.py script): 

```

barcode01_PorechopsOut_MergedFromPorechop_Binned_10nt_Annotated.txt,1
barcode02_PorechopsOut_MergedFromPorechop_Binned_10nt_Annotated.txt,1
barcode03_PorechopsOut_MergedFromPorechop_Binned_10nt_Annotated.txt,1
barcode04_PorechopsOut_MergedFromPorechop_Binned_10nt_Annotated.txt,1
barcode05_PorechopsOut_MergedFromPorechop_Binned_10nt_Annotated.txt,2
barcode06_PorechopsOut_MergedFromPorechop_Binned_10nt_Annotated.txt,2
barcode07_PorechopsOut_MergedFromPorechop_Binned_10nt_Annotated.txt,2
barcode08_PorechopsOut_MergedFromPorechop_Binned_10nt_Annotated.txt,2
barcode09_PorechopsOut_MergedFromPorechop_Binned_10nt_Annotated.txt,3
barcode10_PorechopsOut_MergedFromPorechop_Binned_10nt_Annotated.txt,3
barcode11_PorechopsOut_MergedFromPorechop_Binned_10nt_Annotated.txt,3
barcode12_PorechopsOut_MergedFromPorechop_Binned_10nt_Annotated.txt,3
barcode13_PorechopsOut_MergedFromPorechop_Binned_10nt_Annotated.txt,4
barcode14_PorechopsOut_MergedFromPorechop_Binned_10nt_Annotated.txt,4
barcode15_PorechopsOut_MergedFromPorechop_Binned_10nt_Annotated.txt,4
barcode16_PorechopsOut_MergedFromPorechop_Binned_10nt_Annotated.txt,4
barcode17_PorechopsOut_MergedFromPorechop_Binned_10nt_Annotated.txt,5
barcode18_PorechopsOut_MergedFromPorechop_Binned_10nt_Annotated.txt,5
barcode21_PorechopsOut_MergedFromPorechop_Binned_10nt_Annotated.txt,5
barcode22_PorechopsOut_MergedFromPorechop_Binned_10nt_Annotated.txt,5


```

Example of running it, obs it is written in python3! Tested in environment ```py3```

```

python /medstore/projects/P23-044/Code/Viral_Integration_Pipeline/NanoPore_pipeline/FilterViralOutForPCRSpilling_SCOPE.py --metadata metadata_scope.txt --sample SCOPE

```

Outputs from SCOPE PCR filtering are:

* 'samplename'_DiscardedHitsFromContamination_Bin10.txt
  - Containing the integrations removed from the PCR filtering step.  
* 'samplename'_MergedParsedFilteredFromContamination_Bin10.txt
  - Containing the integrations kept after the PCR filtering, with this file you can see which samples having this integration and which group they belong to.
* 'samplename'_PCRfilt.txt
  - Is just the samples after PCR filtering in the same format as the input files *Binned_10nt_Annotated.txt

Another important thing for the PCR filtering, if we are at the same amount of reads we are keeping the integration. Only remove if one is higher than the other, this is the case for the old PCR filtering as well! 


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

