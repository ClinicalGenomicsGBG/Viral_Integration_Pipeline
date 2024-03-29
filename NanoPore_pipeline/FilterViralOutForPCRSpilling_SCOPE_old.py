#!/usr/bin/py

# This script takes all final annotated outputs from binned and annotated files and filters breakpoints found in multiple samples, this is an contamination. The script saves the breakpoint with most reads, then outputs filtered files together with a table telling with breakpoints have been removed from the different samples. 


# Here we use a bin of + -10 from the samples, the easiest thing to do is to merge the files and then loop throught them sorted by amount of reads. Using the highest evidence as the reference level 


import sys
import argparse
import os
import logging
import time
import psutil
import glob


#source activate py3k

def parseArgs(argv):
    '''
    Parsing the arguments
    '''
    parser = argparse.ArgumentParser(description='SCOPE PCR filtering, groups the bins 4 and 4. Takes the highest coverage in a sample within the bin, this is the breakpoint, saves the bins and its sum of depth. Compares it to the other bins')
    parser.add_argument("--metadata", dest = 'metadata', required=True,help ="Comma seperated file containing the samplename,group")
    parser.add_argument("--sample", dest = 'sample', required=True,help ="Sample name for the output")
    arguments = parser.parse_args(argv)
    return arguments

def ParseTheScopeFiles(metadata):
    """
    Read the metadata, for the same group merge the files to one merge,
    """
    logging.info('%s\tCreating one merged integration files per group, sorted by nReads',time.ctime().split(" ")[-2])

    ToRemove=[]
    
    outMerged=[]
    Grouped={}
    with open(metadata, "r") as inf:
        for l in inf:
            l=l.strip()
            filename=l.split(",")[0]
            group=l.split(",")[1]

            if not group in Grouped:
                Grouped[group]=[filename]
            else:
                Grouped[group].append(filename)
   
    Tot=0
    for key, values in Grouped.items():
        outmerged= key + "_Merged.tmp"
        ToRemove.append(outmerged)
        with open(outmerged, "w") as o:
            for i in values:
                with open(i, "r") as inf:
                    next(inf)
                    for l in inf:
                        l=l.strip()
                        t=int(l.split("\t")[0])
                        Tot+=t
                        print(i + "\t" + l, file=o)
        outmergedsorted= key + "_Merged_Sorted.tmp"
        ToRemove.append(outmergedsorted)
        command = "sort -k 2,2nr " + outmerged +  " > " + outmergedsorted
        os.system(command)
        outMerged.append(outmergedsorted)
    return(Tot,outMerged, ToRemove)



def readMergedFromScope(outMerged, sample, ToRemove):
    """
    Here we read the merged files to get a merged file sorted by the sum reads per group. It is a tab delimited file containing the sum of reads across the integration, integration coord. Group and all the integrations within this integration supporting it. 

    What happens if there are two ndepth 1s that support a break. We take the one sorted in the table first, Ok in the first round as we are not discarding anything, just parsin to get a coord for support! 

    """

    logging.info('%s\tParsing Files within the Sample groups',time.ctime().split(" ")[-2])    
    ParsedFiles=[]
    for i in outMerged:
        d={}
        coords=[]
        counter=0
        with open(i, "r") as firstround: # First open it to save the coords to have something to comapare it with
            detec=[]
            for l1 in firstround:
                counter+=1
                l1=l1.strip()
                l1coord=l1.split("\t")[2]+"-" + l1.split("\t")[3]
                l1lchrom=l1.split("\t")[2].split(":")[0]
                l1rchrom=l1.split("\t")[3].split(":")[0]
                l1lpos=int(l1.split("\t")[2].split(":")[-1])
                l1rpos=int(l1.split("\t")[3].split(":")[-1])
                l1sample=l1.split("\t")[0]
                l1nreads=int(l1.split("\t")[1])
                l1anno1=l1.split("\t")[4]
                l1anno2=l1.split("\t")[5]
                if not l1 in detec: 
                    detec.append(l1)
                    d[l1coord]=[[l1nreads,l1lchrom,l1lpos,l1rchrom,l1rpos,l1anno1,l1anno2, l1sample]]
                    with open(i, "r") as secondround: # Second loop through the same file to be able to compare the previous hit to all the others, if within rane 
                        for l2 in secondround:
                            l2=l2.strip()
                            l2coord=l2.split("\t")[2]+"-" + l2.split("\t")[3]
                            l2lchrom=l2.split("\t")[2].split(":")[0]
                            l2rchrom=l2.split("\t")[3].split(":")[0]
                            l2lpos=int(l2.split("\t")[2].split(":")[-1])
                            l2rpos=int(l2.split("\t")[3].split(":")[-1])
                            l2sample=l2.split("\t")[0]
                            l2nreads=int(l2.split("\t")[1])
                            l2anno1=l2.split("\t")[4]
                            l2anno2=l2.split("\t")[5]
                            if not l2 in detec: 
                                if l1lchrom==l2lchrom and l1rchrom == l2rchrom:
                                    if abs(l1lpos - l2lpos) <= 10 and abs(l1rpos - l2rpos) <= 10:
                                        detec.append(l2)
                                        d[l1coord].append([l2nreads,l2lchrom,l2lpos,l2rchrom,l2rpos,l2anno1,l2anno2, l2sample])
        
        tot=0
        group=i.split("_Merged_Sorted.tmp")[0]
        outfile=group+"_Parsed.tmp"
        ToRemove.append(outfile)
        with open(outfile, "w") as o:
            for key, values in d.items():
                tot+=len(values)
                suma=0
                for s in values:
                    suma+=s[0] # suma is the sum or nreads per integration bin
                #print(suma, key, group,";".join(values), file=o)
                IntegrationString=""
                for v in values:
                    if not IntegrationString:
                        IntegrationString = ",".join(str(x) for x in v)
                    else:
                        IntegrationString=IntegrationString +"|" +",".join(str(x) for x in v)
                print(str(suma) + "\t" + str(key) + "\t" + str(group) + "\t" + IntegrationString, file=o)
                
        logging.info('%s\tnIntegrations (unique) in %s: %s, nIntegrations recorded (unique):%s ', time.ctime().split(" ")[-2], i, counter, tot)
        ParsedFiles.append(outfile)


    out_parsedunfilt=sample+"_Mega_Merged_Parsed_unfilt.tmp"
        
    command = "cat " + " ".join(ParsedFiles) +" | sort -k 1,1nr >  %s" %out_parsedunfilt
    os.system(command)
    
    return(out_parsedunfilt, ToRemove)
        
        

def ParseAcrossTheSampleGroupsScope(Tot,out_parsedunfilt, sample, ToRemove):
    """
    
    We can loop through this merged sorted file again twice just as we did above, just keep track on the time.. it might be very slow!

    Loop trough the rows, if there is one within  +-10 coords range, if there is different groups within the same range but with the same nreads keep them both! Only discard when one is higher! 

    Outputs are the common FiltForContaminations_groups_Bin.txt nreads\tlcorrds\trcoord\tanno1\tanno2
    x The ParsedFiles output with removed contaminants! (New)
    x FilteredHitsFromContamination_Bin10.txt, which barcodes were filtered? sample\tnreads\tlcorrds\trcoord\tanno1\tanno2

    """

    logging.info('%s\tComparing the integrations between the groups',time.ctime().split(" ")[-2])
    OutParsed=sample+"_MergedParsedFilteredFromContamination_Bin10.txt"
    OutDiscarded=sample+"_DiscardedHitsFromContamination_Bin10.txt"
    
    with open(out_parsedunfilt, "r") as firstround, open(OutParsed, "w") as oMerged, open(OutDiscarded, "w") as oDiscard: # First open it to save the coords to have something to comapare it with
        detec=[]
        #discard=[]
        for l1 in firstround:
            l1=l1.strip()
            l1coord=l1.split("\t")[1]
            l1lchrom=l1.split("\t")[1].split(":")[0]
            l1rchrom=l1.split("\t")[1].split("-")[-1].split(":")[0]
            l1lpos=int(l1.split("\t")[1].split("-")[0].split(":")[-1])
            l1rpos=int(l1.split("\t")[1].split("-")[-1].split(":")[-1])           
            l1nreads=int(l1.split("\t")[0])
            if not l1 in detec: # As l1 comes first, if it is not in detect always print to the parsed file! 
                detec.append(l1)
                print(l1, file=oMerged)
                with open(out_parsedunfilt, "r") as secondround: # Second loop through the same file to be able to compare the previous hit to all the others, if within rane
                    for l2 in secondround:
                       l2=l2.strip()
                       l2coord=l2.split("\t")[1]
                       l2lchrom=l2.split("\t")[1].split(":")[0]
                       l2rchrom=l2.split("\t")[1].split("-")[-1].split(":")[0]
                       l2lpos=int(l2.split("\t")[1].split("-")[0].split(":")[-1])
                       l2rpos=int(l2.split("\t")[1].split("-")[-1].split(":")[-1])
                       l2nreads=int(l2.split("\t")[0])
                       if not l2 in detec:
                            if l1lchrom==l2lchrom and l1rchrom == l2rchrom:
                                if abs(l1lpos - l2lpos) <= 10 and abs(l1rpos - l2rpos) <= 10:
                                    detec.append(l2)
                                    if l2nreads < l1nreads: # If l2 is within the same range as l1 but smaller we discard it!
                                        filtstring=l2.split("\t")[3].replace("|","\n").replace(",","\t")
                                        print(filtstring, file=oDiscard)
                                    else:
                                        print(l2, file=oMerged) # if l2 is not smaller it can just mean it is the same size, then we still keep it!


    # Gathering statistics to make sure the input integrations is the same as the sum of reported (filtered and Discarded)
    TotFilt=0
    with open(OutParsed, "r") as inf:
        for l in inf:
            l=l.strip()
            c=int(l.split("\t")[0])
            TotFilt+=c

    TotDiscarded=0
    with open(OutDiscarded, "r") as inf:
        for l in inf:
            l=l.strip()
            c=int(l.split("\t")[0])
            TotDiscarded+=c
            
    logging.info('%s\tnFiltered integrations (Reads):%s \tnDiscarded Integrations (Reads): %s', time.ctime().split(" ")[-2],TotFilt, TotDiscarded)
    logging.info('%s\tnTotal output integrations (Reads): %s\t nTotal input integrations (Reads): %s', time.ctime().split(" ")[-2],TotFilt+TotDiscarded, Tot)

    # Parsing the report, one filtered per barcode, extracted from the OutParsed!
    # As we are parsing with echo (and append) we need to remove outputfiles if they are there
    """
    currentdir=os.getcwd()
    print(currentdir)
    res = glob.glob(currentdir+"/*_FiltForContaminationBetweenSameGroups_Bin10.txt")
    if res:
        logging.info('%s\tWarning, there is already outputs from SCOPE in Folder, removing these now!', time.ctime().split(" ")[-2])
        for i in res:
            os.remove(i)
    """

    Detected=[]
    with open(OutParsed, "r") as inf:
        for l in inf:
            l=l.strip()
            for r in l.split("\t")[3].split("|"):
                row=r.replace(",","\t")
                sample=row.split("\t")[-1].split("_Binned_10nt_Annotated.txt")[0]+"_FiltForContaminationBetweenSameGroups_Bin10.txt"
                outrow=row.split("\t")[0]+"\t"+ row.split("\t")[1] + ":" + row.split("\t")[2] + "\t" +  row.split("\t")[3] + ":" +  row.split("\t")[4] +"\t" +  row.split("\t")[5] + "\t" + row.split("\t")[6]
                if not sample in Detected: # This is the first time we reach it! Create the file and the header
                    command="echo \"nSupportedReads\tleftmostCoords\trightmostCoords\tAnno\tAnno2\" > " + sample
                    os.system(command)
                    command="echo \""+outrow+ "\" >> "+sample
                    os.system(command)
                    Detected.append(sample)
                else:
                    command="echo \""+outrow+ "\" >> "+sample
                    os.system(command)


    logging.info('%s\tCleaning', time.ctime())

    # Remove the tmps
    for r in ToRemove:
        os.remove(r)
                    

                    
def main(metadata, sample):
    # Create the logger
    process = psutil.Process(os.getpid())
    logging.basicConfig(level=logging.INFO)
    start = time.time()
    logging.info('%s\tStart PCR filtering in SCOPE', time.ctime())
    (Tot,outMerged,ToRemove)=ParseTheScopeFiles(metadata)
    (out_parsedunfilt,ToRemove)=readMergedFromScope(outMerged, sample,ToRemove)
    ParseAcrossTheSampleGroupsScope(Tot,out_parsedunfilt, sample, ToRemove)
    logging.info('%s\tEnd PCR filtering in SCOPE', time.ctime())

if __name__ == '__main__':
    args=parseArgs(sys.argv[1:])
    main(args.metadata, args.sample)

