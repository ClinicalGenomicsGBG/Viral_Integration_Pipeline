#!/usr/bin/py 

# This scripts takes the output from SingleEnd viral integration and bins the data and adding a human annotation (Read within gene, intron exonic intergenic)


# This is the script called binned for illumina in the beginning as we for nanoport only went for exact breakpoints. When the amount of reads in the nanopore increased we bin these as 10 as well so the annotation binning step is the same for short read and long reads

import sys 

import os 

import pandas

import itertools

# The annotaton 

def readdata():
    HumanGTF = "/jumbo/db/Homo_sapiens/Ensembl/GRCh37.87/Annotation/Homo_sapiens.GRCh37.87.CHR.gtf"
    AnnovarScript="/apps/bio/apps/annovar/20150322/annotate_variation.pl"
    Anndb="/apps/bio/apps/annovar/20150322/humandb"
    nanovirusfus=sys.argv[1]
    genometoskip=sys.argv[2] # We need this so we are not annotating something that is not in 
    outname=sys.argv[3]
    return(HumanGTF,AnnovarScript,Anndb,nanovirusfus,genometoskip,outname)



def binning(nanovirusfus,outname):
    """
    This part bins the fusions, this is a sperate output file as we only include the exact breaks here
    """

    tmp = outname + ".tmp"
    tmp2 = outname + "_S.tmp"
    with open(tmp, "w") as o:  
        with open(nanovirusfus, "r") as lines:
            next(lines)
            for l in lines:
                l=l.strip()
                leftBreak = l.split("\t")[1] +":"+l.split("\t")[3] 
                rigthBreak = l.split("\t")[8] +":"+l.split("\t")[9]
                print >> o, leftBreak + "\t" + rigthBreak
    command = "sort -n %s | uniq -c  |  sed 's/^ *//g' | sed 's/ /\t/g' | sort -k 1,1nr > %s" %(tmp,tmp2)
    os.system(command)
    command = "rm %s" %tmp
    os.system(command)
    return tmp2

def binningAdvance(tmp2, outname):
    """
    This part is unique to Illumina as it has more reads, it will take the highest amount of reads across and select 10 nt on each side of it. When it does this it also needs to keep track on the left side! 
    """
    tmp3 = outname + "_D.tmp"
    dataf=pandas.read_csv(tmp2,sep="\t", header=None)
    dataf[['chrom', 'pos']] = dataf[2].str.split(':', 1, expand=True)
    dataf[['Virus', 'Viruspos']] = dataf[1].str.split(':', 1, expand=True)
    #dataf.sort_values(['pos','Viruspos'], ascending=[True, True], inplace=True)
    dataf.sort_values([0],inplace=True, ascending=False)
    with open(tmp3, "w") as o: 
        grouped=dataf.groupby('chrom')
        for name, group in grouped:
            viruslist=[]
            humanlist=[]
            countlist=[]
            detected=[]
            #print group
            for row_index, row in group.iterrows():
                vp=int(row["Viruspos"])
                hp=int(row["pos"])
                c=int(row[0])
                viruslist.append(vp)
                humanlist.append(hp)
                countlist.append(c)
            counter=-1
            for vp,hp,c in itertools.izip(viruslist,humanlist, countlist):
                #print "HBV" + ":" +str(vp) + "\t" +  str(row['chrom'])+":"+ str(hp) + "\t" + str(c)
                counter+=1
                indexsecondloop=0 
                indexsecondlooplist=[]
                #if counter == len(group) and searchstring in detected:
                #   if len(indexsecondlooplist)==
                #  str(countlist[counter-1]) + "\t"+ row['Virus'] + ":" + str(vp) + "\t" + str(row['chrom']) + ":" + str(hp)
                for vp2, hp2,c2 in itertools.izip(viruslist[counter:],humanlist[counter:], countlist[counter:]):
                    indexsecondloop+=1
                    if vp2-10 <= vp <= vp2+10 and hp2-10 <= hp <= hp2+10:
                        indexsecondlooplist.append(indexsecondloop+counter)
                del indexsecondlooplist[0]
                searchstring=row['Virus'] + ":" + str(vp) + "-" + str(row['chrom']) + ":" + str(hp)
                if searchstring in detected:
                    #print "!"
                    continue
                else:
                    detected.append(searchstring)
                    if not indexsecondlooplist:
                        print >>o, str(c) + "\t"+ row['Virus'] + ":" + str(vp) + "\t" + str(row['chrom']) + ":" + str(hp)
                        #detected.append(searchstring)
                    else: 
                        if len(indexsecondlooplist) == 1:

                            searchstring = row['Virus'] +":"+ str(viruslist[indexsecondlooplist[0]-1]) + "-" + str(row['chrom']) + ":" + str(humanlist[indexsecondlooplist[0]-1])
                            if not searchstring in detected:
                                print >>o, str(c+countlist[indexsecondlooplist[0]-1]) + "\t"+ row['Virus'] + ":" + str(vp) + "\t" + str(row['chrom']) + ":" + str(hp)

                            else:
                               print >>o, str(c) + "\t"+ row['Virus'] + ":" + str(vp) + "\t" + str(row['chrom']) + ":" + str(hp)
                               
                            detected.append(searchstring)
                        if len(indexsecondlooplist) > 1: 
                            sums=0
                            for i in indexsecondlooplist:
                                i=i-1
                                searchstring = row['Virus'] +":"+ str(viruslist[i]) + "-" + str(row['chrom']) + ":" + str(humanlist[i])
                                if not searchstring in detected: 
                                    sums+=int(countlist[i])
                                detected.append(searchstring)
                            print >>o, str(c+sums)  + "\t"+ row['Virus'] + ":" +  str(vp) + "\t" + str(row['chrom']) + ":" + str(hp)            
    return tmp3
            

def Annotating(tmp3,outname,genometoskip,AnnovarScript,Anndb):
    """
    The part annotates the bins 
    """
    AnnotatedFile=outname + "_Binned_10nt_Annotated.txt"
    tmp = outname + "Annovar.tmp"
    tpmanno=outname + "Annovar2.tmp"
    with open(AnnotatedFile, "w") as o: 
        print >>o,"nSupportedReads\tleftmostCoords\trightmostCoords\tAnno\tAnno2"
        with open(tmp3, "r") as inf: 
            for l in inf: 
                with open(tmp, "w") as otmp: 
                    l=l.strip()
                    leftmost=l.split("\t")[1]
                    rightmost=l.split("\t")[2]
                    if genometoskip in leftmost:
                        coord=rightmost
                        print >> otmp, coord.split(":")[0] + "\t" + coord.split(":")[1] + "\t" + coord.split(":")[1] + "\t0\t0"
                        # The virus is not in this coord
                    elif genometoskip in rightmost: # If the virus is in leftmost we look at the rightmost   
                        coord=leftmost
                        print >> otmp, coord.split(":")[0] + "\t" + coord.split(":")[1] + "\t" + coord.split(":")[1] + "\t0\t0"
                
                CommandRunAnnovar = "perl %s --geneanno --buildver hg19 %s %s -out %s" %(AnnovarScript, tmp, Anndb, tpmanno)
                os.system(CommandRunAnnovar)
                annovarout=tpmanno+".variant_function"
                with open(annovarout,"r") as variantf:
                    for lines in variantf:
                        strippedl=lines.strip()
                        generalinsert = strippedl.split("\t")[0]
                        geneinsertanno= strippedl.split("\t")[1]
                        print >>o, l + "\t" + generalinsert + "\t" + geneinsertanno
    command="rm *.tmp*"
    os.system(command)

def main():
    (HumanGTF,AnnovarScript,Anndb,nanovirusfus,genometoskip,outname)=readdata()
    tmp2=binning(nanovirusfus,outname)
    tmp3=binningAdvance(tmp2, outname)
    Annotating(tmp3,outname,genometoskip,AnnovarScript,Anndb)

main()
