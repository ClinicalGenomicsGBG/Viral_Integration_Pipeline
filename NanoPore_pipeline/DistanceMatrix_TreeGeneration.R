#!/usr/bin/R

suppressMessages(library(DECIPHER))
suppressMessages(library(pheatmap))
suppressMessages(library(ape))
suppressMessages(library(phangorn))
suppressMessages(library(ggtree))


# Reads the MA 
# Generates the distance matrix 
# Creates trees in UPGMA and NJ 
# We need to create the distance matrix using DECIPHER as in python i could not find a way of handle gaps!

args = commandArgs(trailingOnly=TRUE)
MAFile=args[1]

# Important, we are not including 
MA<-readDNAStringSet(MAFile, format="fasta", use.names=TRUE)
d<-DistanceMatrix(MA,type='matrix', includeTerminalGaps=FALSE, penalizeGapLetterMatches = TRUE)
d[is.na(d)]<-1 # Set NA which means no shared identity to 1


# The warning is from the X from the AlignmentFile, set them as N in the consensus generation!


if (nrow(d) < 100){
sizes=15
}else{
sizes=nrow(d)/5
}

outfile=gsub(".fa","_Distance.txt",MAFile)
write.table(d, file=outfile,sep="\t", quote=F)


outpdf=gsub(".fa","_Distance.pdf",MAFile)
pdf(outpdf,height=sizes, width = sizes)
pheatmap(d)
dev.off()



if (nrow(d) < 1000){
sizes=50
}else{
sizes=nrow(d)/5
}


# NJ tree
outfile=gsub(".fa","_NJ.newick", MAFile)
tr<-nj(d)
write.tree(tr, file=outfile)
tree<-read.tree(outfile)

# Fix Meta
#Meta<-replace(tree$tip.label, grep("unintegrated",tree$tip.label, invert=TRUE), "Integrated")
#Meta<-replace(Meta, grep("unintegrated",tree$tip.label, invert=FALSE), "Unintegrated")
#Meta<-as.data.frame(Meta)
#Meta$Tip<-tree$tip.label



outpdf=gsub(".fa","_NJ_tree.pdf", MAFile)
p<-ggtree(tree, branch.length='none', layout='circular') + geom_tiplab(size=3, color="black")
pdf(outpdf,width=sizes, height=sizes)
p
dev.off()


# UPGMA Will be rooted
my_upgma <- phangorn::upgma(d)
outfile=gsub(".fa","_UPGMA.newick", MAFile)
write.tree(my_upgma, file=outfile)
tree<-read.tree(outfile)

outpdf=gsub(".fa","_UGPMA_tree.pdf", MAFile)
p<-ggtree(tree, branch.length='none', layout='circular') + geom_tiplab(size=3, color="black")
pdf(outpdf,width=sizes, height=sizes)
p
dev.off()


#ggsave("testTree.svg", plot=p, width = 60, height = 20, limitsize=FALSE)
