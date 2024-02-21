#!/usr/bin/env Rscript

library(edgeR)

count_file <- "/mnt/projects/usersfeatureCounts_results/count_matrix_toGenome.txt"
count_data <- read.table(count_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

print(ncol(count_data))

#here we will take columns 7-18 as experimental data and 19-30 as norm data.
gene_id <- count_data$Geneid
samples_normal <- count_data[, c(7:18)]
samples_experiment <- count_data[, c(19:30)]

dge_list <- DGEList(counts = cbind(samples_normal, samples_experiment), genes = gene_id)

group <- factor(rep(c("Normal", "Experiment"), each = ncol(samples_normal)))

dge_list <- calcNormFactors(dge_list)
design_matrix <- model.matrix(~0 + group)
dge_list <- estimateDisp(dge_list, design_matrix)
fit <- glmQLFit(dge_list, design_matrix)
contrast <- makeContrasts(groupExperiment - groupNormal, levels = colnames(design_matrix))
qlf <- glmQLFTest(fit, contrast = contrast)

results <- topTags(qlf, n = 25)

# Save all results to a file
write.table(results$table, file = "EDGER_full_results.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Save summary information about top 25 differentially expressed genes
top_genes_summary <- results$table[1:25, c("genes", "logFC", "logCPM", "PValue", "FDR")]
write.table(top_genes_summary, file = "EDGER_topgenes_output.txt", sep = "\t", quote = FALSE, row.names = FALSE)
