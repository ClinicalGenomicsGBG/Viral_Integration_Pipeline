#!/usr/bin/py

# Script that takes a nanopore alignment file, mapped using minimap2. Using the SA flag it finds chimeric reads and uses them to define a break point, look at the SA to check if iy is at the positive or negative strand

# If it is a mix take the negative strand read and substract the coords from the read length to make it comparable to the positive strand 

# Here we take a the HBV and set it to leftmost coord! Ignoring the identification in version3! We are expecting that HBV sets the start so the other fusion is calculated from HBV! 
# Be careful of the warnings, if we get warnings we fusion reads are discarded, need to change in the code to make sure it is captured


import pysam
import sys

import argparse
from Bio.Seq import Seq

import re

# Start by building it for HBV, We can change this later on!

def parseArgs(argv):
    '''
    Parsing the arguments
    '''
    parser = argparse.ArgumentParser(description='Takes an alignment file, nanopore SE, mapped with minimap2. Using the SA flag to detect fusions')
    parser.add_argument("-b", dest = 'bam',help ="bamfile for extraction")
    parser.add_argument("-g", dest = 'genome',help ="identifier that you use for genome")
    arguments = parser.parse_args(argv)
    return arguments


def extractChimerasINGenomeOfInterest(bam,genome):
    """
    This part takes the bam and extracts all chimeras 
    """
    samf=pysam.AlignmentFile(bam, "rb")
    chimeradict={}
    ReadsacrossGenomeOfInterest=samf.fetch(genome)
    outfile=bam.split(".bam")[0] + ".OutSE_FusionPipe.txt"
    for read in ReadsacrossGenomeOfInterest:
        if read.has_tag("SA"):
            if 5 in read.cigartuples[0]:
                startquery=read.query_alignment_start+read.cigartuples[0][-1]
                endquery=read.query_alignment_end+startquery
            else:
                startquery=read.query_alignment_start
                endquery=read.query_alignment_end
            if not read.is_reverse:
                chimeradict[read.query_name]=[startquery,endquery,"+", read.reference_name,read.reference_start,read.reference_end,read.query_alignment_sequence, read.infer_read_length()]
            else: 
                chimeradict[read.query_name]=[startquery,endquery,"-", read.reference_name,read.reference_start,read.reference_end,read.query_alignment_sequence,read.infer_read_length()]
    chimeradictnotingenomeofinterest={}
    allreads=samf.fetch()
    for read in allreads:
         if read.has_tag("SA"):
             if not read.reference_name == genome:
                 if 5 in read.cigartuples[0]: # The read is hardclipped and therefore we need to add it to the query alignment calculation otherwies the hardclipped is just ignored and start is at 0
                     startquery=read.query_alignment_start+read.cigartuples[0][-1]
                     endquery=read.query_alignment_end+startquery
                 else: # We are not touching the breakpoint or alignment at the reference position as we only using this position to see where it aligns
                     startquery=read.query_alignment_start
                     endquery=read.query_alignment_end
                 if read.query_name in chimeradict.keys():
                     if not read.is_reverse:
                         chimeradictnotingenomeofinterest[read.query_name]=[startquery,endquery, "+",read.reference_name,read.reference_start,read.reference_end,read.query_alignment_sequence, read.infer_read_length()]
                     else:
                        chimeradictnotingenomeofinterest[read.query_name]=[startquery,endquery, "-",read.reference_name,read.reference_start,read.reference_end,read.query_alignment_sequence, read.infer_read_length()]
    with open(outfile, "w") as o:
        
        print >> o, "Readname\tLeftmostChrom\tLeftmostStart\tLeftmostEnd\tLeftmoststrand\tLeftmostSequence\tLeftmostReadcoordstart\tLeftmostReadcoordend\tRightmostChrom\tRightmostStart\tRightmostEnd\tRightmoststrand\tRightmostSequence\tRightmostReadcoordstart\tRightmostReadcoordend\tComment"
        Fusions=[] # Save to list so you can track which ones that does not follow our set tresholds
        DetectedFusion=[]
        for i in chimeradictnotingenomeofinterest.keys():
            Fusions.append(i)
            start = chimeradict[i][-1]-chimeradict[i][1]
            end =  chimeradict[i][-1]-chimeradict[i][0]
            start_notchromofinterest = chimeradictnotingenomeofinterest[i][-1]-chimeradictnotingenomeofinterest[i][1]
            end_notchromofinterest =  chimeradictnotingenomeofinterest[i][-1]-chimeradictnotingenomeofinterest[i][0]
            if chimeradict[i][2] == "+" and chimeradictnotingenomeofinterest[i][2] == "+": # Both are + 
                if chimeradict[i][0] < chimeradictnotingenomeofinterest[i][0]: # Virus is the first integration, how it should be 
                    print >> o, str(i) + "\t" + chimeradict[i][3] + "\t" + str(chimeradict[i][4]) + "\t" + str(chimeradict[i][5]) + "\t" +  chimeradict[i][2] + "\t" + chimeradict[i][6] + "\t" + str(chimeradict[i][0]) + "\t" + str(chimeradict[i][1]) + "\t"+  chimeradictnotingenomeofinterest[i][3] + "\t" + str(chimeradictnotingenomeofinterest[i][4]) + "\t" + str(chimeradictnotingenomeofinterest[i][5]) + "\t" +  chimeradictnotingenomeofinterest[i][2] + "\t" + chimeradictnotingenomeofinterest[i][6] + "\t" + str(chimeradictnotingenomeofinterest[i][0]) + "\t" + str(chimeradictnotingenomeofinterest[i][1]) + "\tCorrect"
                    DetectedFusion.append(i)
                #else: # HBV is on plus but it is not smaller then the human fusion, human fusion is probably on minus calc to make it comparable to HBV 
                 #   print >> o, str(i) + "\t" + chimeradict[i][3] + "\t" + str(chimeradict[i][4]) + "\t" + str(chimeradict[i][5]) + "\t" +  chimeradict[i][2] + "\t" + chimeradict[i][6] + "\t" + str(chimeradict[i][0]) + "\t" + str(chimeradict[i][1]) + "\t"+  chimeradictnotingenomeofinterest[i][3] + "\t" + str(chimeradictnotingenomeofinterest[i][4]) + "\t" + str(chimeradictnotingenomeofinterest[i][5]) + "\t" +  chimeradictnotingenomeofinterest[i][2] + "\t" + chimeradictnotingenomeofinterest[i][6] + "\t" + str(chimeradictnotingenomeofinterest[i][0]) + "\t" + str(chimeradictnotingenomeofinterest[i][1])+ "\tCorrect"
                  #  DetectedFusion.append(i)
            if chimeradict[i][2] == "-" and chimeradictnotingenomeofinterest[i][2] == "-": # If both are on minus you need to calculate both from the length 
                if start < start_notchromofinterest: # Virus is the first coord when counting from the behind 
                    print >> o, str(i) + "\t" + chimeradict[i][3] + "\t" + str(chimeradict[i][4]) + "\t" + str(chimeradict[i][5]) + "\t" +  chimeradict[i][2] + "\t" + chimeradict[i][6] + "\t" + str(start) + "\t" + str(end) + "\t"+  chimeradictnotingenomeofinterest[i][3] + "\t" + str(chimeradictnotingenomeofinterest[i][4]) + "\t" + str(chimeradictnotingenomeofinterest[i][5]) + "\t" +  chimeradictnotingenomeofinterest[i][2] + "\t" + chimeradictnotingenomeofinterest[i][6] + "\t" + str(start_notchromofinterest) + "\t" + str(end_notchromofinterest)+ "\tCorrect"
                    DetectedFusion.append(i)
                
            # Then we go for one virus is positive and the other one is negative 
            if chimeradict[i][2] == "+" and chimeradictnotingenomeofinterest[i][2] == "-":
                if chimeradict[i][0] < start_notchromofinterest: 
                    print >> o, str(i) + "\t" + chimeradict[i][3] + "\t" + str(chimeradict[i][4]) + "\t" + str(chimeradict[i][5]) + "\t" +  chimeradict[i][2] + "\t" + chimeradict[i][6] + "\t" + str(chimeradict[i][0]) + "\t" + str(chimeradict[i][1]) + "\t"+  chimeradictnotingenomeofinterest[i][3] + "\t" + str(chimeradictnotingenomeofinterest[i][4]) + "\t" + str(chimeradictnotingenomeofinterest[i][5]) + "\t" +  chimeradictnotingenomeofinterest[i][2] + "\t" + chimeradictnotingenomeofinterest[i][6] + "\t" + str(start_notchromofinterest) + "\t" + str(end_notchromofinterest)+ "\tCorrect"
                    DetectedFusion.append(i)
            # Final option that virus is negative and the other one is positive
            if chimeradict[i][2] == "-" and chimeradictnotingenomeofinterest[i][2] == "+":
                if start < chimeradictnotingenomeofinterest[i][0]:
                    print >> o, str(i) + "\t" + chimeradict[i][3] + "\t" + str(chimeradict[i][4]) + "\t" + str(chimeradict[i][5]) + "\t" +  chimeradict[i][2] + "\t" + chimeradict[i][6] + "\t" + str(start) + "\t" + str(end) + "\t"+  chimeradictnotingenomeofinterest[i][3] + "\t" + str(chimeradictnotingenomeofinterest[i][4]) + "\t" + str(chimeradictnotingenomeofinterest[i][5]) + "\t" +  chimeradictnotingenomeofinterest[i][2] + "\t" + chimeradictnotingenomeofinterest[i][6] + "\t" + str(start_notchromofinterest) + "\t" + str(end_notchromofinterest)+ "\tCorrect"
                    DetectedFusion.append(i)
            # Final loop to correct the ones that are wrong... These does not behave as they should with the virus being to the leftmost coord...
        if not len(DetectedFusion) == len(Fusions):
            for i in Fusions:
                start = chimeradict[i][-1]-chimeradict[i][1]
                end =  chimeradict[i][-1]-chimeradict[i][0]
                start_notchromofinterest = chimeradictnotingenomeofinterest[i][-1]-chimeradictnotingenomeofinterest[i][1]
                end_notchromofinterest =  chimeradictnotingenomeofinterest[i][-1]-chimeradictnotingenomeofinterest[i][0]                
                if not i in DetectedFusion:
                    # start in human - end in virus, take the one with the largest distance! 
                    tester=[(end, start_notchromofinterest), (chimeradict[i][1], chimeradictnotingenomeofinterest[i][0]), (end, chimeradictnotingenomeofinterest[i][0])]# Distance from start and end 
                    tester_overlap =[(range(start,end), range(start_notchromofinterest,end_notchromofinterest)),(range(chimeradict[i][0],chimeradict[i][1]), range(chimeradictnotingenomeofinterest[i][0], chimeradictnotingenomeofinterest[i][1])), (range(start, end), range(chimeradictnotingenomeofinterest[i][0],chimeradictnotingenomeofinterest[i][1]))]
                    #print tester_overlap
                    res = []
                    res_overlap=[]
                    for t in tester: 
                        res.append(abs(t[-1]-t[0]))
                    for t in tester_overlap: # This one checks the amount of overlap, lest start extracting the coord were the distance is the smallest between the breaks! 
                        z=[x for x in t[0] if x in t[-1]]
                        res_overlap.append(len(z))
                    index_min = min(range(len(res)), key=res.__getitem__) # Get value closest to 0 
                    if index_min == 0:
                        print >> o, str(i) + "\t" + chimeradict[i][3] + "\t" + str(chimeradict[i][4]) + "\t" + str(chimeradict[i][5]) + "\t" +  chimeradict[i][2] + "\t" + chimeradict[i][6] + "\t" + str(start) + "\t" + str(end) + "\t"+  chimeradictnotingenomeofinterest[i][3] + "\t" + str(chimeradictnotingenomeofinterest[i][4]) + "\t" + str(chimeradictnotingenomeofinterest[i][5]) + "\t" +  chimeradictnotingenomeofinterest[i][2] + "\t" + chimeradictnotingenomeofinterest[i][6] + "\t" + str(start_notchromofinterest) + "\t" + str(end_notchromofinterest) + "\tExtraCheck" 
                    if index_min == 1: 
                        print >> o, str(i) + "\t" + chimeradict[i][3] + "\t" + str(chimeradict[i][4]) + "\t" + str(chimeradict[i][5]) + "\t" +  chimeradict[i][2] + "\t" + chimeradict[i][6] + "\t" + str(chimeradict[i][0]) + "\t" + str(chimeradict[i][1]) + "\t"+  chimeradictnotingenomeofinterest[i][3] + "\t" + str(chimeradictnotingenomeofinterest[i][4]) + "\t" + str(chimeradictnotingenomeofinterest[i][5]) + "\t" +  chimeradictnotingenomeofinterest[i][2] + "\t" + chimeradictnotingenomeofinterest[i][6] + "\t" + str(chimeradictnotingenomeofinterest[i][0]) + "\t" + str(chimeradictnotingenomeofinterest[i][1]) + "\tExtraCheck"
                    if index_min == 2:
                        print >> o, str(i) + "\t" + chimeradict[i][3] + "\t" + str(chimeradict[i][4]) + "\t" + str(chimeradict[i][5]) + "\t" +  chimeradict[i][2] + "\t" + chimeradict[i][6] + "\t" + str(start) + "\t" + str(end) + "\t"+  chimeradictnotingenomeofinterest[i][3] + "\t" + str(chimeradictnotingenomeofinterest[i][4]) + "\t" + str(chimeradictnotingenomeofinterest[i][5]) + "\t" +  chimeradictnotingenomeofinterest[i][2] + "\t" + chimeradictnotingenomeofinterest[i][6] + "\t" + str(chimeradictnotingenomeofinterest[i][0]) + "\t" + str(chimeradictnotingenomeofinterest[i][1]) + "\tExtraCheck"
                        
                    #else: 
                     #   print chimeradict[i], chimeradictnotingenomeofinterest[i]
                      #  print "Warning Missing Fusion!! Add extra step to capture it!" 
    samf.close()



"""
            # Start to see which read starts before the other
            if chimeradict[i][2]=="+" and chimeradictnotingenomeofinterest[i][2] =="+": # Both are on plus, very easy to treat! 
                if chimeradict[i][0] < chimeradictnotingenomeofinterest[i][0]: # Genome of interest is in the first coord
                    print >> o, str(i) + "\t" + chimeradict[i][3] + "\t" + str(chimeradict[i][4]) + "\t" + str(chimeradict[i][5]) + "\t" +  chimeradict[i][2] + "\t" + chimeradict[i][6] + "\t" + str(chimeradict[i][0]) + "\t" + str(chimeradict[i][1]) + "\t"+  chimeradictnotingenomeofinterest[i][3] + "\t" + str(chimeradictnotingenomeofinterest[i][4]) + "\t" + str(chimeradictnotingenomeofinterest[i][5]) + "\t" +  chimeradictnotingenomeofinterest[i][2] + "\t" + chimeradictnotingenomeofinterest[i][6] + "\t" + str(chimeradictnotingenomeofinterest[i][0]) + "\t" + str(chimeradictnotingenomeofinterest[i][1]) 
                else: 
                    print >> o, str(i) + "\t" +  chimeradictnotingenomeofinterest[i][3] + "\t" + str(chimeradictnotingenomeofinterest[i][4]) + "\t" + str(chimeradictnotingenomeofinterest[i][5]) + "\t" +  chimeradictnotingenomeofinterest[i][2] + "\t" + chimeradictnotingenomeofinterest[i][6] + "\t" + str(chimeradictnotingenomeofinterest[i][0]) + "\t" + str(chimeradictnotingenomeofinterest[i][1]) + "\t" + chimeradict[i][3] + "\t" + str(chimeradict[i][4]) + "\t" + str(chimeradict[i][5]) + "\t" +  chimeradict[i][2] + "\t" + chimeradict[i][6] + "\t" + str(chimeradict[i][0]) + "\t" + str(chimeradict[i][1])
            elif chimeradict[i][2]=="-" and chimeradictnotingenomeofinterest[i][2] =="-": # Both are on minus
                if chimeradict[i][0] < chimeradictnotingenomeofinterest[i][0]: # Genome of interest is in the first coord
                    print >> o, str(i) + "\t" + chimeradict[i][3] + "\t" + str(chimeradict[i][4]) + "\t" + str(chimeradict[i][5]) + "\t" +  chimeradict[i][2] + "\t" + chimeradict[i][6] + "\t" + str(chimeradict[i][0]) + "\t" + str(chimeradict[i][1]) + "\t"+  chimeradictnotingenomeofinterest[i][3] + "\t" + str(chimeradictnotingenomeofinterest[i][4]) + "\t" + str(chimeradictnotingenomeofinterest[i][5]) + "\t" +  chimeradictnotingenomeofinterest[i][2] + "\t" + chimeradictnotingenomeofinterest[i][6] + "\t" + str(chimeradictnotingenomeofinterest[i][0]) + "\t" + str(chimeradictnotingenomeofinterest[i][1]) 
                else: 
                    print >> o, str(i) + "\t" +  chimeradictnotingenomeofinterest[i][3] + "\t" + str(chimeradictnotingenomeofinterest[i][4]) + "\t" + str(chimeradictnotingenomeofinterest[i][5]) + "\t" +  chimeradictnotingenomeofinterest[i][2] + "\t" + chimeradictnotingenomeofinterest[i][6] + "\t" + str(chimeradictnotingenomeofinterest[i][0]) + "\t" + str(chimeradictnotingenomeofinterest[i][1])+"\t" + chimeradict[i][3] + "\t" + str(chimeradict[i][4]) + "\t" + str(chimeradict[i][5]) + "\t" +  chimeradict[i][2] + "\t" + chimeradict[i][6] + "\t" + str(chimeradict[i][0]) + "\t" + str(chimeradict[i][1])
            else: 
                # One is on plus and one is on the negative strand, you still use the coord, only first start. In the example where they overlapp se still 
                if chimeradict[i][2]=="-":
                    start = chimeradict[i][-1]-chimeradict[i][1]
                    end =  chimeradict[i][-1]-chimeradict[i][0]
                    if start <  chimeradictnotingenomeofinterest[i][0]:
                        print >> o, str(i) + "\t" + chimeradict[i][3] + "\t" + str(chimeradict[i][4]) + "\t" + str(chimeradict[i][5]) + "\t" +  chimeradict[i][2] + "\t" + chimeradict[i][6] + "\t"+ str(start) + "\t" + str(end) + "\t" +  chimeradictnotingenomeofinterest[i][3] + "\t" + str(chimeradictnotingenomeofinterest[i][4]) + "\t" + str(chimeradictnotingenomeofinterest[i][5]) + "\t" + chimeradictnotingenomeofinterest[i][2] + "\t" + chimeradictnotingenomeofinterest[i][6] + "\t" + str(chimeradictnotingenomeofinterest[i][0]) + "\t" + str(chimeradictnotingenomeofinterest[i][1])

                    else: 
                        print >>o, str(i) + "\t" +  chimeradictnotingenomeofinterest[i][3] + "\t" + str(chimeradictnotingenomeofinterest[i][4]) + "\t" + str(chimeradictnotingenomeofinterest[i][5]) + "\t" + chimeradictnotingenomeofinterest[i][2] + "\t" + chimeradictnotingenomeofinterest[i][6] + "\t" + str(chimeradictnotingenomeofinterest[i][0]) + "\t" + str(chimeradictnotingenomeofinterest[i][1]) + "\t" + chimeradict[i][3] + "\t" + str(chimeradict[i][4]) + "\t" + str(chimeradict[i][5]) + "\t" +  chimeradict[i][2] + "\t" + chimeradict[i][6] + "\t"+ str(start) + "\t" + str(end)
                else: # It is not the genome of interest that is on the negative strand
                    start=chimeradictnotingenomeofinterest[i][-1]-chimeradictnotingenomeofinterest[i][1]
                    end=chimeradictnotingenomeofinterest[i][-1]-chimeradictnotingenomeofinterest[i][0]
                    
                    if start < chimeradict[i][0]:
                        print >> o, str(i) + "\t" + chimeradictnotingenomeofinterest[i][3] + "\t" + str(chimeradictnotingenomeofinterest[i][4]) + "\t" + str(chimeradictnotingenomeofinterest[i][5]) + "\t" + chimeradictnotingenomeofinterest[i][2] + "\t" + chimeradictnotingenomeofinterest[i][6] + "\t" + str(start) + "\t" + str(end) + "\t" + chimeradict[i][3] + "\t" + str(chimeradict[i][4]) + "\t" + str(chimeradict[i][5]) + "\t" +  chimeradict[i][2] + "\t" + chimeradict[i][6] + "\t"+ str(chimeradict[i][0]) +  str(chimeradict[i][1])
                    else: 
                        print >> o, str(i) + "\t" + chimeradict[i][3] + "\t" + str(chimeradict[i][4]) + "\t" + str(chimeradict[i][5]) + "\t" +  chimeradict[i][2] + "\t" + chimeradict[i][6] + "\t" + str(chimeradict[i][0]) + "\t" + str(chimeradict[i][1]) + "\t" +  chimeradictnotingenomeofinterest[i][3] + "\t" + str(chimeradictnotingenomeofinterest[i][4]) + "\t" + str(chimeradictnotingenomeofinterest[i][5]) + "\t" +chimeradictnotingenomeofinterest[i][2] + "\t" + chimeradictnotingenomeofinterest[i][6] + "\t" + str(start) + "\t" + str(end)


"""

    

def main(bam,genome):
    extractChimerasINGenomeOfInterest(bam,genome)

if __name__ == '__main__':
    args=parseArgs(sys.argv[1:])
    main(args.bam, args.genome)



