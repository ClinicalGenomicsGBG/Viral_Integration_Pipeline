#!/usr/bin/py

library("treeio")
library("ggtree")

#nwk<-system.file("extdata","/medstore/projects/P23-044/Intermediate/SCOPE-Batch1_HBVintegrationPipeline/VirusPipelineOut/barcode04/barcode04_Merged_NJ.newick", package="treeio")


tree<-read.tree("/medstore/projects/P23-044/Intermediate/SCOPE-Batch1_HBVintegrationPipeline/VirusPipelineOut/barcode04/barcode04_Merged_NJ.newick")

ggtree(tree, layout='circular') + xlim(-10, NA)

D
