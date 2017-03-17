library("edgeR")

files <- c(

"ERR458507.fastq.gz.quant.counts",
"ERR458508.fastq.gz.quant.counts",
"ERR458509.fastq.gz.quant.counts",    
"ERR458500.fastq.gz.quant.counts",
"ERR458501.fastq.gz.quant.counts",
"ERR458502.fastq.gz.quant.counts",
"ERR458493.fastq.gz.quant.counts",
"ERR458494.fastq.gz.quant.counts",
"ERR458495.fastq.gz.quant.counts",
"ERR458878.fastq.gz.quant.counts",
"ERR458879.fastq.gz.quant.counts",
"ERR458880.fastq.gz.quant.counts"   
)

labels=c("2.MT.507","2.MT.508","2.MT.509","1.MT.500","1.MT.501","1.MT.502","1.WT.493","1.WT.494","1.WT.495","2.WT.878","2.WT.879","2.WT.880")

data <- readDGE(files)

print(data)

###
#6 mutants go first and then the 6 WT

group <- c(rep("mut", 6), rep("wt", 6))

dge = DGEList(counts=data, group=group)
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)

# make an MA-plot

et <- exactTest(dge, pair=c("mut", "wt"))
etp <- topTags(et, n=100000)
etp$table$logFC = -etp$table$logFC
pdf("yeast-edgeR-MA-plot.pdf")
plot(
  etp$table$logCPM,
  etp$table$logFC,
  xlim=c(-3, 20), ylim=c(-12, 12), pch=20, cex=.3,
  col = ifelse( etp$table$FDR < .2, "red", "black" ) )
dev.off()

# plot MDS
pdf("yeast-edgeR-MDS.pdf")
plotMDS(dge, labels=labels)
dev.off()

# output CSV for 0-6 hr
write.csv(etp$table, "yeast-edgeR.csv")
