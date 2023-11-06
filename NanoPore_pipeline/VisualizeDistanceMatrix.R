#!/usr/binx


# Visualizing the distance matrixes! 
# Clustering is observed amount of mutations, 
# Phylogeny estimates the expected mutations, Specially an issue with reversion mutations where one position mutatated twice and went back to the original state, Uncorrected distance will measure this as 0 (observed), while the expected should be 2! 


library(pheatmap)
args = commandArgs(trailingOnly=TRUE)
DistanceMatrix=args[1]
outpdf=gsub(".csv",".pdf",DistanceMatrix)
dist<-read.csv(DistanceMatrix, header=T, row.names=1)

if (nrow(dist) < 100){
sizes=15
}else{
sizes=nrow(dist)/5
}


pdf(outpdf,height=sizes, width = sizes)
pheatmap(dist)
dev.off()


#library(qgraph)
#d<-as.dist(dist)
#mds.coor <- cmdscale(d)
#head(mds.coor)
# Multidimensional scaling!
#plot(mds.coor[,1], mds.coor[,2], type="n", xlab="", ylab="")
#text(jitter(mds.coor[,1]), jitter(mds.coor[,2]),
#     rownames(mds.coor), cex=0.8)
#abline(h=0,v=0,col="gray75")

