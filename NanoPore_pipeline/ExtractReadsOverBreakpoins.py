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



# The problem with extracting the fusions, we are taking the ones with the highest evidence and goes from 10 on each side for those, Therefore one read can be calculated multiple times which is wrong! We need to order by amount of reads, starting from them to go through the line, and if a read was already matched ignore it.

# python /medstore/projects/P23-044/Code/Viral_Integration_Pipeline/NanoPore_pipeline/ExtractReadsOverBreakpoins.py --BinningFile barcode04_PorechopsOut_MergedFromPorechop_Binned_10nt_Annotated.txt --OutFromFusionPipe barcode04_PorechopsOut_MergedFromPorechop.OutSE_FusionPipe.txt  --Output barcode04

 #~/.conda/envs/ClustalOmega/bin/clustalo

# How should the consensus be generated? Now 0.7, ambigous sequence of N

 
def parseArgs():
    parser = argparse.ArgumentParser(description='Takes the binning output and the bamfile to extract the reads across the fusions',formatter_class=RawTextHelpFormatter)
    parser.add_argument('--BinningFile', dest='BinningFile', help='path to the binning file', required=True)
    parser.add_argument('--OutFromFusionPipe', dest='OFusionPipe', help='path to the output from fusion', required=True)
    parser.add_argument('--ClustalOexec', dest='ClustalOexec', help='Path to ClustalOmega', required=True)
    parser.add_argument('--Output', dest='Output', help='Output filename, sets the basename', required=True)
    arguments=parser.parse_args(sys.argv[1:])
    return arguments


def ExtractTheCoordinates(BinningFile):
    """
    The binning is for 10 at each side so +10 -10 of the actual coord to be included into the consensus
    This part extracts the binning coordinates already generated. 
    """

    print("--- Extracting the bin Coordinates ---")
    
    coords={}
    
    with open(BinningFile, "r") as B:
        next(B)
        key=0
        for l in B:
            key+=1
            l=l.strip()
            reads=int(l.split("\t")[0])
            #Viruscoord=l.split("\t")[1].split(":")[0] + ":" + str(int(l.split("\t")[1].split(":")[1])-10) + "-" + str(int(l.split("\t")[1].split(":")[1])+10)
            #Humancoord=l.split("\t")[2].split(":")[0] + ":" + str(int(l.split("\t")[2].split(":")[1])-10) + "-" + str(int(l.split("\t")[2].split(":")[1])+10)
            Viruscoord=l.split("\t")[1]
            Humancoord=l.split("\t")[2]
            Anno=l.split("\t")[3]+"-" + l.split("\t")[4]
            coords[key]=[reads,Viruscoord,Humancoord,Anno]
    return(coords)



def ExtractFromOUTSE(OFusionPipe,coords):
    """
    From the binnings extract all sequences, these will be used for the clustering
    """
    print("--- Extracting the sequences coupled to the bins ---")

    #print(coords)
    
    #print(sorted(coords.items(), key=lambda e: e[1][0], reverse=True))

    
    Detec=[]
    
    for key, values in sorted(coords.items(), key=lambda e: e[1][0], reverse=True):
        VirusChromosome=values[1].split(":")[0]
        VirusCoord=int(values[1].split(":")[1] )
        HumanChromosome=values[2].split(":")[0]
        HumanCoord=int(values[2].split(":")[1])
        with open(OFusionPipe, "r") as inf:
            next(inf)
            for l in inf:
                l=l.strip()
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
                                coords[key].append(combinedseq)
                                Detec.append(identifier)


    with open(OFusionPipe, "r") as inf:
        next(inf)
        nRows=0
        for l in inf:
            nRows+=1
                                
    print("Amount of fusions reads: ", nRows)
    print("Amount of reads used in clustering: ",len(Detec)) # These numbers should be the same
    return(coords)
                

def CreateConsensus(BinningFile, Output,ClustalOexec,coords):
    """
    This one performs the multiple sequence alignments and generates the consenus within each bin
    """
    try:
        os.makedirs(Output)
    except OSError: # folder exists, write to it.
        print("* Folder Exists, write to it")

    Singletons=[]
    ForMA=[]
    for key, values in coords.items():
        coord=values[1].replace(":","_")+"_"+values[2].replace(":","_")
        outfasta=Output+"/"+Output+"_"+coord+".fasta"
        
        with open(outfasta, "w") as o: 
            if len(values[4:]) == 1: # We only have one sequence, this will just be our Consens
                print(">"+Output+"_"+coord+"\n"+values[4], file=o)
                Singletons.append(outfasta)
            else: # We have more than one, need to perform the MA for these
                counter=0
                for seq in values[4:]: 
                    counter+=1
                    print(">"+str(counter)+"\n"+seq, file=o)
                    if not outfasta in ForMA: # You should only append to it once! 
                        ForMA.append(outfasta)

    print("Amount of singleton fusions",len(Singletons))

    print("--- Performing MA (For consenus generation)---")
    # Generates the multiple alignment for the consenus sequences

    Consenussequences=[]
    for f in ForMA:
        outf=f.strip(".fasta")+".phylip"

        command = "%s -i %s -o %s --outfmt phylip --force --threads 2" %(ClustalOexec, f, outf)
        #subprocess.call(command, shell=True)
        alignment = AlignIO.read(outf, "phylip")
        summary_align = AlignInfo.SummaryInfo(alignment)
        consensus = summary_align.dumb_consensus(threshold=0.7, ambiguous='X') 
        outf_consensus=f.strip(".fasta")+"_consensus.fasta"
        Consenussequences.append(outf_consensus)
        with open(outf_consensus, "w") as o:
            outheader=outf.split("/")[1].strip(".phylip")
            print(">"+outheader+"\n"+consensus, file=o)

    # Generates the merged fasta with singleTons and consensus

    OutMerged=Output+"/"+Output+"_Merged.fasta"


    
    with open(OutMerged, "w") as o:
        # print the singletons to file
        for i in Singletons:
            with open(i, "r") as s:
                for l in s:
                    l=l.strip()
                    print(l, file=o)
        # print the consensus to file
        for i in Consenussequences:
            with open(i, "r") as c:
                for l in c:
                    l=l.strip()
                    print(l, file=o)


    with open(BinningFile, "r") as inf:
        next(inf)
        uniqueBreaks=0
        for l in inf:
            uniqueBreaks+=1


    # These numbers should be the same
    print("Total unique breaks: " + str(uniqueBreaks))
    print("Total seq in fasta Merged: "+ str(len(Singletons)+len(Consenussequences)))
                    
    return(OutMerged)        


def GenerateTree(Output, ClustalOexec, OutMerged):
    """
    Here we generates the tree using the different consensus sequences 

    * The distance matrix is based on the multiple alignment, measures the identity distance between the sequences, for example first sequence towards first sequence is 0 as they are identitcal, if the alignment length is 13 and there are 3 differences between a pair of sequence the distance is 3/13=0.23

    * UPGMA is a unweigthed pair group with arithmetic mean. Sequential clustering that starts with things most similair. It results in a rooted tree, assumes that the rate of the evolution is the same among all organisms 

    * Neighbour-joining uses the average distances to other leaves as well, it produces an unrooted tree. Works fairly well in practice. As Neighbour joinning allows for unequal rates of evolution the branch lengths are proportional to the amount of change.  (I think Neighbour joining is the one people prefer to use!)
    
    
    """


    print("--- Performing MA (For tree generation)---")
    Outfa=OutMerged.strip(".fasta")+"_MA.fa"
    command = "%s -i %s -o %s --outfmt fasta --force --threads 4" %(ClustalOexec, OutMerged, Outfa)
    #subprocess.call(command, shell=True)
    align = AlignIO.read(Outfa,'fasta')


    print("--- Calculates the distance matrix ---")
    # Calculate the distance matrix
    calculator = DistanceCalculator('identity')
    distMatrix = calculator.get_distance(align)

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
    
def main(BinningFile,OFusionPipe, ClustalOexec, Output):
    coords=ExtractTheCoordinates(BinningFile)
    coords=ExtractFromOUTSE(OFusionPipe,coords)
    OutMerged=CreateConsensus(BinningFile, Output,ClustalOexec,coords)
    GenerateTree(Output, ClustalOexec, OutMerged)

    
if __name__=='__main__':
    arguments=parseArgs()
    main(arguments.BinningFile, arguments.OFusionPipe,arguments.ClustalOexec, arguments.Output)
    


