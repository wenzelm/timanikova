library(dada2)

args = commandArgs(trailingOnly=TRUE)
indir=args[1]
outdir=args[2]
taxdb=args[3]
maxEE=args[4]
dir.create(outdir)

samp <- list.files(indir, "*.fastq.gz")
raw <- file.path(indir, samp)

# quality plots
pdf(width=45, height=15, file=file.path(outdir, paste0("quality_plots.EE", maxEE, ".pdf")))
	plotQualityProfile(raw[sapply(raw, file.size)>100])
dev.off()

# filter and trim
filtered <- file.path(paste0(outdir, "/trimmed.EE", maxEE), samp)
out <- filterAndTrim(raw, filtered, maxN=0, maxEE=maxEE, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
out <- cbind(out, Percent=100*out[,2]/out[,1])
write.table(out, file=file.path(outdir, paste0("filtered.EE", maxEE, ".txt")), sep="\t", row.names=T, col.names=NA, quote=F)

# error rates
filtered <- filtered[!is.na(sapply(filtered, file.size))]
err <- learnErrors(filtered, multithread = TRUE)
pdf(file=file.path(outdir, paste0("error_rates.EE", maxEE, ".pdf")))
	plotErrors(err, nominalQ=TRUE)
dev.off()

# inference
asv <- dada(filtered, err=err, multithread=TRUE, pool=TRUE, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)

# sequence table
seqtab <- makeSequenceTable(asv)
sink(file.path(outdir, paste0("seqtab.EE", maxEE, ".log")))
dim(seqtab)
table(nchar(getSequences(seqtab)))

# chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
saveRDS(seqtab.nochim, file=file.path(outdir, paste0("seqtab.EE", maxEE, ".RDS")))

# taxonomy
taxa <- assignTaxonomy(seqtab.nochim, taxdb, multithread=T)
saveRDS(taxa , file=file.path(outdir, paste0("taxa.EE", maxEE, ".RDS")))

apply(taxa, 2, function(x) 100*sum(!is.na(x))/length(x))
sink()

# make ASV FASTA file
fasta <- c(paste0(">ASV", 1:ncol(seqtab.nochim)), colnames(seqtab.nochim))[order(rep(1:ncol(seqtab.nochim), times=2))]
writeLines(fasta, con=file.path(outdir, paste0("seqtab.EE", maxEE, ".fasta")))

# tidy up
unlink(dirname(filtered), recursive = TRUE)