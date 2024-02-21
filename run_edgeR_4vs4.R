#!/usr/bin/env Rscript


library("edgeR")

setwd('/mnt/projects/users/aalayeva/smallRNA/smRNA_STAR_RESULT/aligned_deduplicated/featureCounts_results/optional')

fc <- read.delim("count_edited_matrix_toGenome.txt", stringsAsFactors = FALSE)
ncol(fc)
genes <- fc[,1:6]
counts <- data.matrix(fc[,7:30])
row.names(counts) <- paste(genes$Geneid, genes$Start, genes$End, sep=".")
group <- factor( c("FL","FL","HL","HL","FL","FL","HL","HL","FL","FL","HL","HL","FL","FL","HL","HL") )
ncol(group)
ncol(counts)


y <- DGEList(counts, genes=genes, group=group)
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)

#write.table(topTags(et, , file='edjer_output.txt', sep='\t', quote=FALSE, col.names=FALSE)
