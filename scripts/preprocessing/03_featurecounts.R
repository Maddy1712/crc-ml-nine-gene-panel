# Gene-level counting with Rsubread::featureCounts (Methods, Section 2.3).
# Usage: Rscript scripts/preprocessing/03_featurecounts.R <bam_dir> <Homo_sapiens.GRCh38.106.gtf> <out.tsv>
suppressPackageStartupMessages(library(Rsubread))
args <- commandArgs(trailingOnly = TRUE)
bams <- list.files(args[1], pattern = "\\.bam$", full.names = TRUE)
fc <- featureCounts(files = bams, annot.ext = args[2],
                    isGTFAnnotationFile = TRUE, GTF.featureType = "gene", GTF.attrType = "gene_id",
                    useMetaFeatures = TRUE, isPairedEnd = TRUE,
                    strandSpecific = 2,          # reverse-stranded, as stated in the Methods
                    nthreads = 16)
cts <- fc$counts
colnames(cts) <- sub("\\.bam$", "", basename(colnames(cts)))   # make IDs match the metadata sample IDs
out <- data.frame(Geneid = rownames(cts), cts, check.names = FALSE)
write.table(out, args[3], sep = "\t", quote = FALSE, row.names = FALSE)
