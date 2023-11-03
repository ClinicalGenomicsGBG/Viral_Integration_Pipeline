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
import pandas as pd
import re
import pysam

import matplotlib.pyplot as plt
import seaborn as sns

from pycirclize import Circos
from io import StringIO

from ete3 import Tree

import logging
logging.basicConfig(level=logging.INFO)

# The problem with extracting the fusions, we are taking the ones with the highest evidence and goes from 10 on each side for those, Therefore one read can be calculated multiple times which is wrong! We need to order by amount of reads, starting from them to go through the line, and if a read was already matched ignore it.

# python /medstore/projects/P23-044/Code/Viral_Integration_Pipeline/NanoPore_pipeline/ExtractReadsOverBreakpoins.py --BinningFile barcode04_PorechopsOut_MergedFromPorechop_Binned_10nt_Annotated.txt --OutFromFusionPipe barcode04_PorechopsOut_MergedFromPorechop.OutSE_FusionPipe.txt  --Output barcode04

 #~/.conda/envs/ClustalOmega/bin/clustalo

# How should the consensus be generated? Now 0.6, ambigous sequence of N

#For solving the qt.qpa.plugin error message run the following! 
#export QT_QPA_PLATFORM=offscreen

# /medstore/projects/P23-044/Code/Viral_Integration_Pipeline/NanoPore_pipeline/VisualizeDistanceMatrix.R

def parseArgs():
    parser = argparse.ArgumentParser(description='Takes the binning output and the OUTSE to extract the reads across the fusions',formatter_class=RawTextHelpFormatter)
    parser.add_argument('--TargetFolder', dest='TargetFolder', help='path to the folder containing the Sample_MergedParsedFilteredFromContamination_Bin10.txt and the *.OutSE_FusionPipe.txt', required=True)
    parser.add_argument('--RawBams', dest='RawBamsFolder', help='path to the folder containing the alignmentFiles from minimap2', required=True)
    parser.add_argument('--ClustalOexec', dest='ClustalOexec', help='Path to ClustalOmega', required=True)
    parser.add_argument('--Output', dest='Output', help='Output filename, sets the basename', required=True)
    parser.add_argument('--PathToRVizScript', dest='PhyloRScript', help='Path to Script that generates the distances and trees in R', required=True)
    parser.add_argument('--DepthTresh', dest='DepthTresh', help='Depth Treshold for the integrations, this is not applied for the non integrated! (default 1)', type=int, default=1)
    arguments=parser.parse_args(sys.argv[1:])
    return arguments



def extracCoords(TargetFolder, DepthTresh):
    """
    Get the breakouts for each barcode by looping the MergedParsedFilteredFromContamination_Bin10.txt file from the Breakpipeline
    """

    logging.info(" --- Extracting integrations based on depth filter " + str(DepthTresh) + " ---")
    
    coords={}
    
    for f in glob.glob(TargetFolder+"/*_MergedParsedFilteredFromContamination_Bin10.txt"):
        with open(f, "r") as inf:
            next(inf)
            for l in inf:
                l=l.strip()
                groupnreads=l.split("\t")[0]
                groupcoords=l.split("\t")[1]
                groupinfo=l.split("\t")[2]
                SummaryInfo=l.split("\t")[3]


                if int(groupnreads) >= DepthTresh:
                
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

    logging.info(" --- Extracting Non integrated reads ---")

    
    # Generate outputfolder!
    try:
        os.makedirs(Output)
        os.makedirs(Output+"/Intermediates")
        os.makedirs(Output+"/Phylo")
    except OSError: # folder exists, write to it.
        logging.info(" OutputFolder Exists, write to it")


    unintegratedConsensusfiles=[]
    unintegratedConsensusfiles_split={}
    
    for key, value in coords.items():
        rawbam=RawBamsFolder+"/"+key+".bam"
        unintegratedBam=Output+"/Intermediates/"+key+"_Unintegrated_HBV.bam"
        command="%s view %s HBV_D -h | grep -v $'\tSA:' | %s view -Sb - > %s" %(samtoolsexec, rawbam, samtoolsexec, unintegratedBam) # Extract the unintegrated reads by invert grepping the SA:Z
        subprocess.call(command, shell=True)        
        outConsensusUnintegrated_all = Output+"/Intermediates/"+key+"_Unintegrated_HBV_All.consensus"
        command = "%s consensus -l 5000 -c 0.6 -a %s > %s" %(samtoolsexec, unintegratedBam, outConsensusUnintegrated_all) # Generate the consensus with samtools consensus, print all so we get the range of unitegrated things,
        subprocess.call(command, shell=True)
        unintegratedConsensusfiles.append(outConsensusUnintegrated_all)

        outConsensusUnintegrated = Output+"/Intermediates/"+key+"_Unintegrated_HBV.consensus"
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
            outFastaSplitOnParts=Output+"/Intermediates/"+key+"_Unintegrated_HBV_SplitOnBreaks.consensus" 
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
                                    print(">"+key+"_HBV_unintegrated_"+str(startanno)+":"+str(end)+"\n"+seqall[start:end], file=o)
                                first=i
                        prev=i
                # we need to check if we have something in the last loop!, good example of this in barcode 12
                # check the first with our final prev, if the distance is more than 10 print it as well!
                if prev-first > 9:
                    start=first
                    startanno=first+1
                    end = prev
                    print(">"+key+"_HBV_unintegrated_"+str(startanno)+":"+str(end)+"\n"+seqall[start:end], file=o)
            unintegratedConsensusfiles_split[key]=[outFastaSplitOnParts]

    return(unintegratedConsensusfiles_split)
         
def ExtractFromOUTSE(TargetFolder,coords):
    """
    From the binnings extract all sequences, these will be used for the clustering
    """

    
    logging.info(" --- Extracting the sequences coupled to the bins ---")
    
    OUTSEfiles=glob.glob(TargetFolder+"/*.OutSE_FusionPipe.txt")

    coordswithseq={}
    sumofreads=0
    
    sumofdetec=0
    for key, values in coords.items():
        perbarcode=0
        Detec=[] # We should only extract it once! 
        for c in values:
            counts=c[0]
            sumofreads+=counts
            perbarcode+=counts
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
                                    


        
        #print(cwithseq)            
        """                        
            if not c[0] == len(c[4:]):
                print(key)
                print(c[0], c[4:])
        
        """


        sumofdetec+=len(Detec)
         

        if perbarcode != len(Detec):
            logging.warning(" not all reads were used in the clustering for: " + str(key))
            logging.warning(" There should be " + str(perbarcode) + " reads")
            logging.warning(" We are getting " + str(len(Detec)) + " reads")

        #print(key)
        #print(len(Detec))
        #print(perbarcode)
        #sumofdetec+=len(Detec)

        
        
    #if sumofreads != sumofdetec:
    #    logging.error("warning, not all reads were used in the clustering")
    #    logging.error("There should be " + str(sumofreads) + " reads")
    #    logging.error("We are getting " + str(sumofdetec) + " reads")

    return(coordswithseq)
                

def CreateConsensus(Output, ClustalOexec, coordswithseq, unintegratedConsensusfiles_split):
    """
    This one performs the multiple sequence alignments and generates the consenus within each bin

    We are using clustal omega to generate the multiple alignment, 
    
    From this file we are taking dumb consensus 
    """
    #try:
    #    os.makedirs(Output)
    #except OSError: # folder exists, write to it.
    #    print("* Folder Exists, write to it")

    MergedFastaOutputs=[]
        
    for key, values in coordswithseq.items():
        Singletons=[]
        ForMA=[]
        uniqueBreaks=0
        #print(values)
        for c in values:
            uniqueBreaks+=1
            coord=c[1].replace(":","_")+"_"+c[2].replace(":","_")
            outfasta=Output+"/Intermediates/"+key+"_"+coord+".fasta"
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
                            
        logging.info(" --- Performing MA (For consenus generation) for " + key + " --- ")
        # Generates the multiple alignment for the consenus sequences
        
        Consenussequences=[]
        for f in ForMA:
            
            outf=f.split(".fasta")[0]+".phylip"
            command = "%s -i %s -o %s --outfmt phylip --force --threads 4" %(ClustalOexec, f, outf)
            subprocess.call(command, shell=True)
            alignment = AlignIO.read(outf, "phylip")
            summary_align = AlignInfo.SummaryInfo(alignment)
            consensus = summary_align.dumb_consensus(threshold=0.6, ambiguous='N')  # ugly consensus just taking the majority, most common residue need to be atleast the treshold (60 % in this case)
            outf_consensus=f.split(".fasta")[0]+"_consensus.fasta"
            Consenussequences.append(outf_consensus)
            with open(outf_consensus, "w") as o:
                outheader=outf.split("/")[-1].strip(".phylip")
                print(">"+outheader+"\n"+consensus, file=o)

        # Generates the merged fasta with singleTons and consensus
        OutMerged=Output+"/Intermediates/"+key+"_Merged.fasta"
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
        logging.info(" For " + key +" Reported " + str(len(Singletons)+len(Consenussequences)) + " unique breakpoints of " + str(uniqueBreaks))
        #MergedFastaOutputs.append(OutMerged)
        
        if key in unintegratedConsensusfiles_split:
            unintegratedConsensusfiles_split[key].append(OutMerged)
        else:

            unintegratedConsensusfiles_split[key]=[OutMerged]
    return(unintegratedConsensusfiles_split)
        

def GenerateTree(Output, ClustalOexec, unintegratedConsensusfiles_split, PhyloRScript):
    """
    Here we generates the tree using the different consensus sequences 

    * The distance matrix is based on the multiple alignment, measures the identity distance between the sequences, for example first sequence towards first sequence is 0 as they are identitcal, if the alignment length is 13 and there are 3 differences between a pair of sequence the distance is 3/13=0.23

    * UPGMA is a unweigthed pair group with arithmetic mean. Sequential clustering that starts with things most similair. It results in a rooted tree, assumes that the rate of the evolution is the same among all organisms 

    * Neighbour-joining uses the average distances to other leaves as well, it produces an unrooted tree. Works fairly well in practice. As Neighbour joinning allows for unequal rates of evolution the branch lengths are proportional to the amount of change.  (I think Neighbour joining is the one people prefer to use!)
    
    
    """

    Rexec="/apps/bio/software/anaconda2/envs/R-4.3.1/bin/Rscript" # Path to R exec

    for key, values in unintegratedConsensusfiles_split.items():
        MergedIntegrationsAndunintegrated=Output+"/Intermediates/"+key+"_Merged_unIntegrated_Integrated_consensus.fasta"

        counter=0
        
        with open(MergedIntegrationsAndunintegrated, "w") as o: # Merge the consensuses of the integrations and the unintegrated ones
            for i in values:
                with open(i, "r") as inf:
                    for l in inf:
                        l=l.strip()
                        counter+=1
                        if ">" in l:
                            l=">"+l.split("MergedFromPorechop_")[-1]
                        print(l, file=o)

        
        if counter > 4: # We need to have more than 3 sequence to be able to perform the the tree
            logging.info("--- Performing MA (For tree generation) ---")
            Outfa=Output+"/Phylo/"+key+"_MA.fa"
            command = "%s -i %s -o %s --outfmt fasta --force --threads 4" %(ClustalOexec, MergedIntegrationsAndunintegrated, Outfa)
            subprocess.call(command, shell=True)

            # We keep this plotting things in the R scrip, specially important to allow for gaps in the distance matrix! 

            
            
            command = "%s %s %s " %(Rexec,PhyloRScript,Outfa)
            subprocess.call(command, shell=True)
            
            """
            align = AlignIO.read(Outfa,'fasta')
            print("--- Calculates the distance matrix ---")
            # Calculate the distance matrix
            # Additional Visuzalization from R script! 
            calculator = DistanceCalculator('identity')
            distMatrix = calculator.get_distance(align)
            newnames=[]
            for i in distMatrix.names: 
                newnames.append(i)
                
            dismatDF = pd.DataFrame(list(distMatrix), columns=newnames, index=newnames)
            outdistmat=Output+"/Distances/"+key+"_DistanceMatrix.csv"            
            dismatDF.to_csv(outdistmat)

            # Plotting using R script
            command = "%s %s %s"%(Rexec, DistancePlottingScript, outdistmat)
            subprocess.call(command, shell=True)
            
            print("--- Constructs the phylogenetic trees---")
            print("--- UPGMA (rooted)---")
            constructor = DistanceTreeConstructor()
            # Construct the phlyogenetic tree using UPGMA algorithm
            UGMATree = constructor.upgma(distMatrix)
            UPGMAout=Output+"/Phylo/"+key+"_UPGMA.newick"
            UPGMAout_treeplot=Output+"/Phylo/"+key+"_UPGMA.svg"
            Phylo.write(UGMATree, UPGMAout, "newick")
            tree = Phylo.read(UPGMAout, "newick")
            #Initialize circos sector with tree size
            circos = Circos(sectors={"Tree": tree.count_terminals()})
            sector = circos.sectors[0]
            # Plot tree
            track = sector.add_track((30, 100))
            track.tree(tree, leaf_label_size=6)
            fig = circos.plotfig()
            fig.savefig(UPGMAout_treeplot)

            # Construct the phlyogenetic tree using NJ algorithm
            print("--- Neighbor-Joining (unrooted) ---")
            NJTree = constructor.nj(distMatrix)
            
            for node in NJTree.get_nonterminals():
                node.name=None

            NJTreeout=Output+"/Phylo/"+key+"_NJ.newick"
            NJTreeout_treeplot=Output+"/Phylo/"+key+"_NJ.svg"
            Phylo.write(NJTree, NJTreeout, "newick")
            tree = Phylo.read(NJTreeout, "newick")
            #Initialize circos sector with tree size
            circos = Circos(sectors={"Tree": tree.count_terminals()})
            sector = circos.sectors[0]
            # Plot tree
            track = sector.add_track((30, 100))
            track.tree(tree, leaf_label_size=6)
            fig = circos.plotfig()
            fig.savefig(NJTreeout_treeplot)
            """
        else: 
            print("warning, less than 3 sequences in", key, "skipping!")
        

def main(TargetFolder, RawBamsFolder, ClustalOexec, PhyloRScript, Output, DepthTresh):
    coords=extracCoords(TargetFolder, DepthTresh)
    unintegratedConsensusfiles_split=ExtractNonIntegratedHBV(coords,RawBamsFolder,Output)
    coordswithseq=ExtractFromOUTSE(TargetFolder,coords)
    unintegratedConsensusfiles=CreateConsensus(Output, ClustalOexec, coordswithseq, unintegratedConsensusfiles_split)
    GenerateTree(Output, ClustalOexec, unintegratedConsensusfiles_split, PhyloRScript)

    
if __name__=='__main__':
    arguments=parseArgs()
    main(arguments.TargetFolder, arguments.RawBamsFolder, arguments.ClustalOexec, arguments.PhyloRScript, arguments.Output, arguments.DepthTresh)
    


