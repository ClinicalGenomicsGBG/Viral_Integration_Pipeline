#!/usr/bin/py

# This script takes all final annotated outputs from binned and annotated files and filters breakpoints found in multiple samples, this is an contamination. The script saves the breakpoint with most reads, then outputs filtered files together with a table telling with breakpoints have been removed from the different samples. 


# Here we use a bin of + -10 from the samples, the easiest thing to do is to merge the files and then loop throught them sorted by amount of reads. Using the highest evidence as the reference level 


import sys
import argparse
import os



#source activate py3k

def parseArgs(argv):
    '''
    Parsing the arguments
    '''
    parser = argparse.ArgumentParser(description='Takes the input reports from the same run and filters insertions detected in more then one sample, these are most likely due to contamination. Saves the one with most reads')
    parser.add_argument("-i", dest = 'infiles', nargs="+",help ="Annotated outputfiles from the binn and annotated step")
    arguments = parser.parse_args(argv)
    return arguments
     

def mergeAndSortInputfiles(infiles):
    """
    Read the infiles, merge them to a mega file and sort by number of supporting reads 
    """
    with open("Merged.tmp","w") as o:
        for i in infiles: 
            with open(i, "r") as inf:
                next(inf)
                for l in inf:
                    l=l.strip()
                    print(i+"\t"+l, file=o)

    command = "sort -k 2,2nr Merged.tmp > Merged_Sorted.tmp"
    os.system(command)

def readMerged():
    hits={}
    with open("Merged_Sorted.tmp", "r") as inf: 
        for l in inf: 
            l=l.strip()
            lchrom=l.split("\t")[2].split(":")[0]
            rchrom=l.split("\t")[3].split(":")[0]
            lpos=int(l.split("\t")[2].split(":")[-1])
            rpos=int(l.split("\t")[3].split(":")[-1])
            sample=l.split("\t")[0]
            nreads=l.split("\t")[1]
            anno1=l.split("\t")[4]
            anno2=l.split("\t")[5]
            if sample in hits: 
                hits[sample].append([lchrom,lpos,rchrom,rpos,nreads,anno1,anno2])
            else: 
                hits[sample] = [[lchrom,lpos,rchrom,rpos,nreads,anno1,anno2]]
    discardlist=[]
    with open("FilteredHitsFromContamination_Bin10.txt", "w") as o:
        with open("Merged_Sorted_Saved.tmp", "w") as o2:
            with open("Merged_Sorted.tmp", "r") as inf:
                for l in inf:
                    l=l.strip()
                    lchrom=l.split("\t")[2].split(":")[0]
                    rchrom=l.split("\t")[3].split(":")[0]
                    lpos=int(l.split("\t")[2].split(":")[-1])
                    rpos=int(l.split("\t")[3].split(":")[-1])
                    sample=l.split("\t")[0]
                    nreads=int(l.split("\t")[1])
                    anno1=l.split("\t")[4]
                    anno2=l.split("\t")[5]
                    for key, values in hits.items(): 
                        if not key == sample: 
                            for i in values:
                                lchrom_s=i[0]
                                rchrom_s = i[2]
                                lpos_s=i[1]
                                rpos_s=i[3]
                                nreads_s=int(i[4])
                                anno1_s=i[5]
                                anno2_s=i[6]
                                if lchrom_s==lchrom and rchrom == rchrom_s:
                                    if abs(lpos - lpos_s) <= 10 and abs(rpos - rpos_s) <= 10:
                                        if nreads_s < nreads:
                                            discardstring=key+"\t" +str(nreads_s)+"\t"+lchrom_s+":"+str(lpos_s)+"\t"+rchrom_s+":"+str(rpos_s)+"\t"+anno1_s+"\t"+anno2_s
                                            discardlist.append(discardstring)
                                        else: 
                                            continue 
                    if l in discardlist:
                        print(l, file=o)
                    if not l in discardlist:
                        print(l, file=o2)
                    discardlist.append(l)

    command = "sed 's/.txt/_FiltForContaminationfromSameRun_Bin10.txt/g' Merged_Sorted_Saved.tmp -i"
    os.system(command)
    command="awk 'BEGIN{FS=\" \";OFS=\"\\t\"} {print $2,$3,$4,$5,$6>$1}' Merged_Sorted_Saved.tmp"
    os.system(command)
    command = "rm Merged_Sorted.tmp"
    os.system(command)
    command = "rm Merged.tmp"
    os.system(command)
    command = "rm Merged_Sorted_Saved.tmp"
    os.system(command)

                        
                
def main(infiles):
    mergeAndSortInputfiles(infiles)
    readMerged()
if __name__ == '__main__':
    args=parseArgs(sys.argv[1:])
    main(args.infiles)

