library(ggplot2)
library(umap)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("Missing parameters. Usage - 'Rscript plot_umap_kmer_freq.R [kmer_count_file.tsv] [output_dir] [output_prefix]'")
}

kmer_file = args[1]
output_directory = file.path(args[2])
output_prefix = args[3]

print(paste("Reading input file:", kmer_file), quote = FALSE)
kmer_data <- read.delim(kmer_file, row.names = 1)

print("Running UMAP", quote = FALSE)
umap_res <- umap(kmer_data)
umap_layout <- as.data.frame(umap_res$layout)

colnames(umap_layout)[colnames(umap_layout) == "V1"] = "UMAP1"
colnames(umap_layout)[colnames(umap_layout) == "V2"] = "UMAP2"

print(paste("Writing output file: ", output_directory, .Platform$file.sep, output_prefix, "_UMAP.tsv", sep = ""), quote = FALSE)
write.table(umap_layout, sep = "\t", quote = FALSE, col.names = NA, row.names = TRUE, file = file.path(output_directory, paste(output_prefix, "_UMAP.tsv", sep = "")))

ggplot(data = umap_layout, aes(x = UMAP1, y = UMAP2)) + geom_point() + xlab("UMAP1") + ylab("UMAP2") + ggtitle(paste("UMAP Kmer Clustering -", output_prefix)) + theme(plot.title = element_text(face = "bold", hjust = 0.5))
print(paste("Writing output figure: ", output_directory, .Platform$file.sep, output_prefix, "_UMAP.pdf", sep = ""), quote = FALSE)
ggsave(paste(output_prefix, "_UMAP.pdf", sep = ""), path = output_directory, width = 12, height = 12)



message("script finished")
