#!/usr/bin/py

# Script that extracts the reads based on breakpoints from the binning, to know which reads should be part of which consensus,
# Extract those aligned that does not contain a breakpoint as well!
# Then for each breakpoint make a multiple alignment 

# Remember the break is from the end coord of Viral and start for Human, this is recorded in the binning as well!


# The script is using py3, most of the viral pipe is using py2

import sys
import argparse
from argparse import ArgumentParser, RawTextHelpFormatter
import os
import glob
import subprocess
import pysam 
from Bio.Align.Applications import ClustalOmegaCommandline 
from Bio import AlignIO
from Bio import Phylo
from Bio.Align import AlignInfo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
#from Bio import pairwise2
#from Bio.pairwise2 import format_alignment

import re

import pysam



# The problem with extracting the fusions, we are taking the ones with the highest evidence and goes from 10 on each side for those, Therefore one read can be calculated multiple times which is wrong! We need to order by amount of reads, starting from them to go through the line, and if a read was already matched ignore it.

# python /medstore/projects/P23-044/Code/Viral_Integration_Pipeline/NanoPore_pipeline/ExtractReadsOverBreakpoins.py --BinningFile barcode04_PorechopsOut_MergedFromPorechop_Binned_10nt_Annotated.txt --OutFromFusionPipe barcode04_PorechopsOut_MergedFromPorechop.OutSE_FusionPipe.txt  --Output barcode04

 #~/.conda/envs/ClustalOmega/bin/clustalo

# How should the consensus be generated? Now 0.7, ambigous sequence of N

 
def parseArgs():
    parser = argparse.ArgumentParser(description='Takes the binning output and the OUTSE to extract the reads across the fusions',formatter_class=RawTextHelpFormatter)
    parser.add_argument('--TargetFolder', dest='TargetFolder', help='path to the folder containing the Sample_MergedParsedFilteredFromContamination_Bin10.txt and the *.OutSE_FusionPipe.txt', required=True)
    parser.add_argument('--RawBams', dest='RawBamsFolder', help='path to the folder containing the alignmentFiles from minimap2', required=True)
    parser.add_argument('--ClustalOexec', dest='ClustalOexec', help='Path to ClustalOmega', required=True)
    parser.add_argument('--Output', dest='Output', help='Output filename, sets the basename', required=True)
    arguments=parser.parse_args(sys.argv[1:])
    return arguments



def extracCoords(TargetFolder):
    """
    
    """


    coords={}
    
    for f in glob.glob(TargetFolder+"/*_MergedParsedFilteredFromContamination_Bin10.txt"):
        print(f)
        with open(f, "r") as inf:
            next(inf)
            for l in inf:
                l=l.strip()
                groupnreads=l.split("\t")[0]
                groupcoords=l.split("\t")[1]
                groupinfo=l.split("\t")[2]
                SummaryInfo=l.split("\t")[3]

                for item in SummaryInfo.split("|"):
                    barcodereads=int(item.split(",")[0])
                    
                    barcodelchrom=item.split(",")[1]
                    barcodelstart=item.split(",")[2]
                    Viralcoords=str(barcodelchrom)+":"+str(barcodelstart)
                    
                    barcoderchrom=item.split(",")[3]
                    barcoderstart=item.split(",")[4]
                    Humancoords=str(barcoderchrom)+":"+str(barcoderstart)
                    
                    barcodeAnnno1=item.split(",")[5]
                    barcodeAnnno2=item.split(",")[6]

                    Anno=barcodeAnnno1+"|"+barcodeAnnno2
                    
                    barcode=item.split(",")[7].split(".OutSE_FusionPipe_Binned_10nt_Annotated.txt")[0]
                    if not barcode in coords:
                        coords[barcode]=[[barcodereads,Viralcoords,Humancoords,Anno]]
                    else:
                        coords[barcode].append([barcodereads,Viralcoords,Humancoords,Anno])
    return(coords)
    


def ExtractNonIntegratedHBV(coords,RawBamsFolder,Output):
    """
    What if we use the alignment, filter out non integrated reads across HBV and create a consensus from that? 
    """

    samtoolsexec= "/apps/bio/software/samtools/1.16.1/bin/samtools"


    unintegratedConsensusfiles=[]
    unintegratedConsensusfiles_split=[]
    
    for key, value in coords.items():
        rawbam=RawBamsFolder+"/"+key+".bam"
        unintegratedBam=Output+"/"+key+"_Unintegrated_HBV.bam"
        command="%s view %s HBV_D -h | grep -v $'\tSA:' | %s view -Sb - > %s" %(samtoolsexec, rawbam, samtoolsexec, unintegratedBam) # Extract the unintegrated reads by invert grepping the SA:Z
        subprocess.call(command, shell=True)        
        outConsensusUnintegrated_all = Output+"/"+key+"_Unintegrated_HBV_All.consensus"
        command = "%s consensus -l 5000 -c 0.6 -a %s > %s" %(samtoolsexec, unintegratedBam, outConsensusUnintegrated_all) # Generate the consensus with samtools consensus, print all so we get the range of unitegrated things,
        subprocess.call(command, shell=True)
        unintegratedConsensusfiles.append(outConsensusUnintegrated_all)

        outConsensusUnintegrated = Output+"/"+key+"_Unintegrated_HBV.consensus"
        command = "%s consensus -l 5000 -c 0.6 %s > %s" %(samtoolsexec, unintegratedBam, outConsensusUnintegrated) # Generate the consensus with samtools consensus, print all so we get the range of unitegrated things,
        subprocess.call(command, shell=True)

        seqall=""
        seq=""
        with open(outConsensusUnintegrated_all, "r") as inf:
            for l in inf:
                l=l.strip()
                if not l.startswith(">"):
                    seqall=l

        with open(outConsensusUnintegrated, "r") as inf:
            for l in inf:
                l=l.strip()
                if not l.startswith(">"):
                    seq=l
        # This part breaks the consensus fasta into smaller parts to improve the clustering, does not make sense to have one long consensus when performing the alignment
        # We are using the Ns reported from samtools consensus 
        #print(key)
        notNs=[]
        counter=-1
        for nucl in seqall:
            counter=counter+1
            if not nucl == "N":
                notNs.append(counter)
        counter=-1
        if seqall:
            outFastaSplitOnParts=Output+"/"+key+"_Unintegrated_HBV_SplitOnBreaks.consensus" 
            with open(outFastaSplitOnParts, "w") as o: 
                for i in notNs:
                    counter+=1
                    if counter == 0:
                        prev=i
                        first=i
                    else:
                        if not i-prev == 1: # We have an N within the nucleotice sequence
                            if i-prev > 5: # most 5 Ns on row, then we add the index of the Ns
                                start=first
                                end=prev+1
                                startanno=first+1 # I am extracting from 0 coord but for common usage one report as plus 1, that is why we will write it in the report
                                if end-start > 9:
                                    print(">HBV_unintegrated_"+str(startanno)+":"+str(end)+"\n"+seqall[start:end], file=o)
                                first=i
                        prev=i
                # we need to check if we have something in the last loop!, good example of this in barcode 12
                # check the first with our final prev, if the distance is more than 10 print it as well!
                if prev-first > 9:
                    start=first
                    startanno=first+1
                    end = prev
                    print(">HBV_unintegrated_"+str(startanno)+":"+str(end)+"\n"+seqall[start:end], file=o)
            unintegratedConsensusfiles_split.append(outFastaSplitOnParts)
                
    return(unintegratedConsensusfiles)
         
def ExtractFromOUTSE(TargetFolder,coords):
    """
    From the binnings extract all sequences, these will be used for the clustering
    """
    print("--- Extracting the sequences coupled to the bins ---")
    
    OUTSEfiles=glob.glob(TargetFolder+"/*.OutSE_FusionPipe.txt")

    coordswithseq={}
    sumofreads=0
    
    sumofdetec=0
    for key, values in coords.items():
        Detec=[] # We should only extract it once! 
        for c in values:
            counts=c[0]
            sumofreads+=counts

        ## Test this
        sortedvalues=sorted(values, key=lambda x: x[0], reverse=True)
        OutSE=TargetFolder+"/"+key+".OutSE_FusionPipe.txt"
    
        for c in sortedvalues:
            cwithseq=c
            VirusChromosome=c[1].split(":")[0]
            VirusCoord=int(c[1].split(":")[1] )
            HumanChromosome=c[2].split(":")[0]
            HumanCoord=int(c[2].split(":")[1])
            with open(OutSE, "r") as inf:
                next(inf)
                for l in inf:
                    identifier=l.split("\t")[0]
                    ViralFusChrom=str(l.split("\t")[1])
                    ViralBreak=int(l.split("\t")[3])
                    ViralSeq=str(l.split("\t")[5])
                    HumanFusChrom=str(l.split("\t")[8])
                    HumanBreak=int(l.split("\t")[9])
                    HumanSeq=str(l.split("\t")[12])
                    if ViralFusChrom == VirusChromosome and HumanChromosome == HumanFusChrom:
                        if (VirusCoord-10) <= ViralBreak <= (VirusCoord+10):
                            if (HumanCoord-10) <= HumanBreak <= (HumanCoord+10):
                                if not identifier in Detec:
                                    
                                    combinedseq=ViralSeq+HumanSeq
                                    cwithseq.append(ViralSeq)
                                    Detec.append(identifier)
                            
            if not key in coordswithseq:
                coordswithseq[key] = [cwithseq]
            else:
                coordswithseq[key].append(cwithseq)

        sumofdetec+=len(Detec)
        """
        sortedvalues=sorted(values, key=lambda x: x[0], reverse=True)        
        for OutSE in OUTSEfiles:
            if key in OutSE:
                with open(OutSE, "r") as inf:
                    next(inf)
                    for l in inf:
                        identifier=l.split("\t")[0]
                        ViralFusChrom=str(l.split("\t")[1])
                        ViralBreak=int(l.split("\t")[3])
                        ViralSeq=str(l.split("\t")[5])
                        HumanFusChrom=str(l.split("\t")[8])
                        HumanBreak=int(l.split("\t")[9])
                        HumanSeq=str(l.split("\t")[12])
                        for c in sortedvalues:
                            cwithseq=c                            
                            VirusChromosome=c[1].split(":")[0]
                            VirusCoord=int(c[1].split(":")[1] )
                            HumanChromosome=c[2].split(":")[0]
                            HumanCoord=int(c[2].split(":")[1])
                            if ViralFusChrom == VirusChromosome and HumanChromosome == HumanFusChrom:
                                if (VirusCoord-10) <= ViralBreak <= (VirusCoord+10):
                                    if (HumanCoord-10) <= HumanBreak <= (HumanCoord+10):
                                        if not identifier in Detec:
                                            combinedseq=ViralSeq+HumanSeq
                                            cwithseq.append(ViralSeq)
                                            Detec.append(identifier)
                                            
                                            if not key in coordswithseq:
                                                coordswithseq[key] = [cwithseq]
                                            else:
                                                coordswithseq[key].append(cwithseq)               
        """

        
    if sumofreads != sumofdetec:
        print("warning, not all reads were used in the clustering")
        print("There should be", sumofreads,"reads")
        print("We are getting",sumofdetec,"reads")

    return(coordswithseq)
                




def CreateConsensus(Output, ClustalOexec, coordswithseq):
    """
    This one performs the multiple sequence alignments and generates the consenus within each bin

    We are using clustal omega to generate the multiple alignment, 
    
    From this file we are taking dumb consensus 
    """
    try:
        os.makedirs(Output)
    except OSError: # folder exists, write to it.
        print("* Folder Exists, write to it")

    MergedFastaOutputs=[]
        
    for key, values in coordswithseq.items():
        Singletons=[]
        ForMA=[]
        uniqueBreaks=0
        #print(values)
        for c in values:
            uniqueBreaks+=1
            coord=c[1].replace(":","_")+"_"+c[2].replace(":","_")
            outfasta=Output+"/"+key+"_"+coord+".fasta"
            counter=0
            with open(outfasta, "w") as o:
                if len(c[4:]) == 1: # We only have one sequence, this will just be our Consens
                    print(">"+key+"_"+coord+"\n"+c[4], file=o)
                    Singletons.append(outfasta)
                else: # We have more than one, need to perform the MA for these
                    for seq in c[4:]: 
                        counter+=1
                        print(">"+str(counter)+"_"+coord+"\n"+seq, file=o)
                        if not outfasta in ForMA: # You should only append to it once! 
                            ForMA.append(outfasta)
                        
        
        print("--- Performing MA (For consenus generation)---")
        # Generates the multiple alignment for the consenus sequences
        
        Consenussequences=[]
        for f in ForMA:
            
            outf=f.split(".fasta")[0]+".phylip"
            command = "%s -i %s -o %s --outfmt phylip --force --threads 2" %(ClustalOexec, f, outf)
            #subprocess.call(command, shell=True)
            alignment = AlignIO.read(outf, "phylip")
            summary_align = AlignInfo.SummaryInfo(alignment)
            consensus = summary_align.dumb_consensus(threshold=0.6, ambiguous='X')  # ugly consensus just taking the majority, most common residue need to be atleast the treshold (60 % in this case)
            outf_consensus=f.split(".fasta")[0]+"_consensus.fasta"
            Consenussequences.append(outf_consensus)
            with open(outf_consensus, "w") as o:
                outheader=outf.split("/")[1].strip(".phylip")
                print(">"+outheader+"\n"+consensus, file=o)

        # Generates the merged fasta with singleTons and consensus
        OutMerged=Output+"/"+key+"_Merged.fasta"
        with open(OutMerged, "w") as o:
            # print the singletons to file
            for i in Singletons:
                with open(i, "r") as s:
                    for l in s:
                        l=l.strip()
                        print(l, file=o)
                        
            # print the consensus to file
            for i in Consenussequences:
                with open(i, "r") as cons:
                    for l in cons:
                        l=l.strip()
                        print(l, file=o)        
        print("For",key,"Reported", len(Singletons)+len(Consenussequences), "unique breakpoints of", uniqueBreaks)
        MergedFastaOutputs.append(OutMerged)

    return(MergedFastaOutputs)
        

def GenerateTree(Output, ClustalOexec, MergedFastaOutputs):
    """
    Here we generates the tree using the different consensus sequences 

    * The distance matrix is based on the multiple alignment, measures the identity distance between the sequences, for example first sequence towards first sequence is 0 as they are identitcal, if the alignment length is 13 and there are 3 differences between a pair of sequence the distance is 3/13=0.23

    * UPGMA is a unweigthed pair group with arithmetic mean. Sequential clustering that starts with things most similair. It results in a rooted tree, assumes that the rate of the evolution is the same among all organisms 

    * Neighbour-joining uses the average distances to other leaves as well, it produces an unrooted tree. Works fairly well in practice. As Neighbour joinning allows for unequal rates of evolution the branch lengths are proportional to the amount of change.  (I think Neighbour joining is the one people prefer to use!)
    
    
    """



    for f in MergedFastaOutputs:
        print(f)
    
    

        print("--- Performing MA (For tree generation)---")
        Outfa=f.split(".fasta")[0]+"_MA.fa"
        command = "%s -i %s -o %s --outfmt fasta --force --threads 4" %(ClustalOexec, f, Outfa)
        subprocess.call(command, shell=True)
        align = AlignIO.read(Outfa,'fasta')


        print("--- Calculates the distance matrix ---")
        # Calculate the distance matrix
        calculator = DistanceCalculator('identity')
        distMatrix = calculator.get_distance(align)

        print(distMatrix)

    """
    print("--- Constructs the phylogenetic trees---")
    # Create a DistanceTreeConstructor object
    constructor = DistanceTreeConstructor()
    # Construct the phlyogenetic tree using UPGMA algorithm

    print("--- UPGMA (rooted)---")
    UGMATree = constructor.upgma(distMatrix)
    #UPGMAout=OutMerged.strip(".fasta")+"_UPGMA.phyloxml"
    
    #Phylo.write(UGMATree, UPGMAout, "phyloxml")


    UPGMAout=OutMerged.strip(".fasta")+"_UPGMA.newick"
    Phylo.write(UGMATree, UPGMAout, "newick")

    
    # Construct the phlyogenetic tree using NJ algorithm

    print("--- Neighbor-Joining (unrooted) ---")
    NJTree = constructor.nj(distMatrix)

    for node in NJTree.get_nonterminals():
        node.name=None

    NJTreeout=OutMerged.strip(".fasta")+"_NJ.newick"
    Phylo.write(NJTree, NJTreeout, "newick")
    
    # Draw the phlyogenetic tree
    #Phylo.draw(UGMATree)
    # Draw the phlyogenetic tree using terminal
    #Phylo.draw_ascii(NJTree)
    
    """

def main(TargetFolder, RawBamsFolder, ClustalOexec, Output):
    coords=extracCoords(TargetFolder)
    ExtractNonIntegratedHBV(coords,RawBamsFolder,Output)
    #coordswithseq=ExtractFromOUTSE(TargetFolder,coords)
    #MergedFastaOutputs=CreateConsensus(Output, ClustalOexec, coordswithseq)
    #GenerateTree(Output, ClustalOexec, MergedFastaOutputs)

    
if __name__=='__main__':
    arguments=parseArgs()
    main(arguments.TargetFolder, arguments.RawBamsFolder, arguments.ClustalOexec, arguments.Output)
    


