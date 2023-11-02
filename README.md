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


## NanoPore integration Clustering (SCOPE)

Aim is to see if the integrations differs from non integrated virus genome. 
The idea is to create consensus for reads over the breakpoints and consensus for the non integrated HBV. These are then clustered together, distance matrix is calculated and distance based phylogenetic trees are generated. 

*What does it do in detail*

1. The integration coordinates for each barcode is extracted from the PCR filtering *MergedParsedFilteredFromContamination_Bin10.txt. Each barcode will be a key in coord dictionary. 
   a) Here we can introduce a depth treshold, good to start with perhaps 10. This means that depth is only applied on the integrations, not the unintegrated parts! 
2. Extract non integrated HBV and create a consensuns: 
   a) Extract all reads that does not contain a SA flag (unintegrated)
   b) Generates a consensus of all unintegrated HBV using samtools. Requires atleast 0.6 in evidence, no depth treshold. Print all positions
   c) We rerun *B* but printing only those with coverage. 
   d) Using *b, c* we can split the consensus if there are more than 5 Ns ina row. The sequences should be ateast 10 nt in length.
   e) The final resulting file is called **_Unintegrated_HBV_SplitOnBreaks.consensus.**
   f) Path to this file is saved to a dictionary, barcode is key 
3. Extract the integrated HBV: 
   a) We are looping through the sorted coord dictionary (by nreads) generated in 1.  
   b) The reads are extracted from *OutSE_FusionPipe.txt, only usese the viral path of the integration. 
   c) The extracted sequences are added to the same coord dictonary as 2. 
4. Create consenus for the breakpoints
   a) The seqeunce singletons are saved to **barcode*coord*fasta**
   b) Clustal omega is used to generate a MA for those with multiple reads across the same breakpoint, consenus with evidence of 0.6 (majority rule) is used to generate the consenus sequences that are saved to **_consensus.fasta**
   c) The consenus and the singletons for each barcode is merged to ***Merged.fasta**. This is the one to be used for Next MA. The path to the merge fasta is saved to the same dictoionary as the unintegrated consenus sequences. 
5. Sequence Clustering (Within each barcode)
   a) The Merged integrations and the non integrated sequences are saved into ***_Merged_unIntegrated_Integrated_consensus.fasta**
   b) Clustal omega is used to generate one MA for each barcode, saved to ***MA.fa**
   c) Distance matrix is calculated in R (DECIPHER), we are removing gaps in the beginning and ignore N. 
   d) A heatmap is generted using the distances 
   e) Two trees are generated using the distance matrix using APE (UPGMA and Neighour-joining), these are plotted using ggtree and saved in newick format.
   f) UPGMA uses sequenctial clustering starting with things that is most similair which results in a rooted tree. Neighbour joining uses average distance using the other leaves, it produces a unrooted tree. 
   
   

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

