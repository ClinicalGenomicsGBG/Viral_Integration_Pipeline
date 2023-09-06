#!/usr/bin/py

# This script takes all final annotated outputs from binned and annotated files and filters breakpoints found in multiple samples, this is an contamination. The script saves the breakpoint with most reads, then outputs filtered files together with a table telling with breakpoints have been removed from the different samples. 

# Here we use a bin of + -10 from the samples, the easiest thing to do is to merge the files and then loop throught them sorted by amount of reads. Using the highest evidence as the reference level 

# Difference to vs1
# Instead of looping through file twice loop through the dictinary keys we create. Will it be faster? 

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
    parser = argparse.ArgumentParser(description='SCOPE PCR filtering, groups the bins using metadatafile. Takes the highest coverage in a sample within the bin, this is the breakpoint, saves the bins and its sum of depth. Compares it to the other bins')
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
        with open(i, "r") as Groups: # First open it to save the coords to have something to comapare it with
            detec=[]
            for l in Groups:
                counter+=1
                l=l.strip()
                lcoord=l.split("\t")[2]+"-" + l.split("\t")[3]
                llchrom=l.split("\t")[2].split(":")[0]
                lrchrom=l.split("\t")[3].split(":")[0]
                llpos=int(l.split("\t")[2].split(":")[-1])
                lrpos=int(l.split("\t")[3].split(":")[-1])
                lsample=l.split("\t")[0]
                lnreads=int(l.split("\t")[1])
                lanno1=l.split("\t")[4]
                lanno2=l.split("\t")[5]
                lanno2=lanno2.replace(",","-") # So we stop having comma in anno2 for easier filtering on biology side

                if not d.keys(): # The dictionary is empty so we either way need to add to the dict
                    d[lcoord]=[[lnreads,llchrom,llpos,lrchrom,lrpos,lanno1,lanno2, lsample]]
                    detec.append(l)
                    
                else:
                    for dcoord in d.keys(): # For every line we check the directory if we are within the rage in the dictionary we append it, and we add it to the detec list
                        dlchrom=dcoord.split(":")[0]
                        drchrom=dcoord.split("-")[-1].split(":")[0]
                        dlpos=int(dcoord.split("-")[0].split(":")[-1])
                        drpos=int(dcoord.split(":")[-1])
                        if dlchrom == llchrom:
                            if drchrom == lrchrom:
                                if abs(llpos - dlpos) <= 10 and abs(lrpos - drpos) <= 10:
                                    if not l in detec: 
                                        detec.append(l)
                                        d[dcoord].append([lnreads,llchrom,llpos,lrchrom,lrpos,lanno1,lanno2, lsample])

                if not l in detec: # if the line is not in the detect list the break is new therefore we add it by itself in the dictory
                    d[lcoord]=[[lnreads,llchrom,llpos,lrchrom,lrpos,lanno1,lanno2, lsample]]
                    detec.append(l)
                    
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

    Important: if there is shared amount of reads in two groups with the same integration keeping them both! 
    
    Outputs are: 
    * DiscardedHitsFromContamination_Bin10.txt : Those were removed from the PCR filtering
    * MergedParsedFilteredFromContamination_Bin10.txt: Filtered with the PCR filtering in a merged file, In this file you can see the amount of reads in each sample as well and which samples that integration is detected in
    * SCOPE_PCRfilt.txt Contains the filtered files in the same format as the input

    Outputs are the common FiltForContaminations_groups_Bin.txt nreads\tlcorrds\trcoord\tanno1\tanno2
    x The ParsedFiles output with removed contaminants! (New)
    x FilteredHitsFromContamination_Bin10.txt, which barcodes were filtered? sample\tnreads\tlcorrds\trcoord\tanno1\tanno2

    """

    logging.info('%s\tComparing the integrations between the groups',time.ctime().split(" ")[-2])
    OutParsed=sample+"_MergedParsedFilteredFromContamination_Bin10.txt"
    OutDiscarded=sample+"_DiscardedHitsFromContamination_Bin10.txt"

    
    with open(out_parsedunfilt, "r") as firstround, open(OutParsed, "w") as oMerged, open(OutDiscarded, "w") as oDiscard: # First open it to save the coords to have something to comapare it with
        detec=[]
        d={}
        #discard=[]
        for l in firstround:
            l=l.strip()
            lcoord=l.split("\t")[1]
            llchrom=l.split("\t")[1].split(":")[0]
            lrchrom=l.split("\t")[1].split("-")[-1].split(":")[0]
            llpos=int(l.split("\t")[1].split("-")[0].split(":")[-1])
            lrpos=int(l.split("\t")[1].split("-")[-1].split(":")[-1])           
            lnreads=str(l.split("\t")[0])
            lgroup=str(l.split("\t")[2])
            lall=l.split("\t")[-1]
            dkey=lcoord+"_"+lgroup+"_"+lnreads
            if not d.keys(): # The dictionary is empty so we either way need to add to the dict
                d[dkey]=[[lgroup+","+lall]]
                detec.append(l)
                
            else:
                for key in d.keys(): # For every line we check the directory if we are within the rage in the dictionary we append it, and we add it to the detec list
                    
                    dlchrom=key.split(":")[0]
                    drchrom=key.split("-")[-1].split(":")[0]
                    dlpos=int(key.split("-")[0].split(":")[-1])
                    drpos=int(key.split(":")[-1].split("_")[0])
                    dnreads=int(key.split("_")[-1])
                    
                    if dlchrom == llchrom:
                        if drchrom == lrchrom:
                            if abs(llpos - dlpos) <= 10 and abs(lrpos - drpos) <= 10:
                                if not l in detec:
                                    if int(lnreads) < int(dnreads): # If the amount of reads is the same as previous we dont want to append it as we are not filtering them away. It will be a seperate hit in the dictionary instead
                                        d[key].append([lgroup+","+lall])                            
                                        detec.append(l)

                                        
            if not l in detec: # if the line is not in the detect list the break is new therefore we add it by itself in the dictory
                detec.append(l)
                d[dkey]=[[lgroup+","+lall]]

        print("nreads_tot\tcoords\tgroup\tsamples", file=oMerged)                
        for key, values in d.items(): # Looping through the key dict to save the outouts:
            ntotalreads=str(key.split("_")[-1])
            ngrouptosave=str(key.split("_")[2])
            lcoord=key.split("-")[0]
            rcoord=key.split("-")[-1].split("_")[0]
            coord=lcoord+"-"+rcoord
            if len(values) == 1: #  only one hit so no need to filter, just print it to the files
                for i in values:
                    details=",".join(i[0].split(",")[1:])
                    print(ntotalreads+"\t"+coord+"\t"+ngrouptosave+"\t"+details, file=oMerged)
                    
                    
            else: # We have multiple groups to the same barcode, here we need to perform filtering!
                # We know which group we are saving from the key, just extract it from the values
                for i in values:
                    group=str(i[0].split(",")[0])
                    if group == ngrouptosave:
                        details=",".join(i[0].split(",")[1:])
                        print(ntotalreads+"\t"+coord+"\t"+ngrouptosave+"\t"+details, file=oMerged)
                    else: # If there are multiple groups to the same bin but the hit group is not the same as the one saved in key we put it in discard!
                        # Obs, we should only discard when one is higher! Otherwise we keep them both!                        
                        for t in ",".join(i[0].split(",")[1:]).split("|"): # We want one sample on each row for quick filtering
                            reads=str(t.split(",")[0])
                            lchrom=str(t.split(",")[1])
                            lcoord=str(t.split(",")[2])
                            rchrom=str(t.split(",")[3])
                            rcoord=str(t.split(",")[4])
                            anno1=str(t.split(",")[5])
                            annot2=str(",".join(t.split(",")[6:-1]))
                            sample=str(t.split(",")[-1])

                            print(sample+"\t"+reads+"\t"+lchrom+":"+lcoord+"\t"+rchrom+":"+rcoord+"\t"+anno1+"\t"+annot2, file=oDiscard)
                            

    # Gathering statistics to make sure the input integrations is the same as the sum of reported (filtered and Discarded)
    TotFilt=0
    with open(OutParsed, "r") as inf:
        next(inf)
        for l in inf:
            l=l.strip()
            c=int(l.split("\t")[0])
            TotFilt+=c

    TotDiscarded=0
    with open(OutDiscarded, "r") as inf:
        for l in inf:
            l=l.strip()
            c=int(l.split("\t")[1])
            TotDiscarded+=c

    logging.info('%s\tnFiltered integrations (Reads):%s \tnDiscarded Integrations (Reads): %s', time.ctime().split(" ")[-2],TotFilt, TotDiscarded)
    logging.info('%s\tnTotal output integrations (Reads): %s\t nTotal input integrations (Reads): %s', time.ctime().split(" ")[-2],TotFilt+TotDiscarded, Tot)


    # Finally we loop through the Merged parsed file to split the integrations to their files
    Detected=[]
    with open(OutParsed, "r") as inf:
        next(inf)
        for l in inf:
            l=l.strip()
            for r in l.split("\t")[3].split("|"):
                row=r.replace(",","\t")
                sample=row.split("\t")[-1].split(".txt")[0]+"_SCOPE_PCRfilt.txt"
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

