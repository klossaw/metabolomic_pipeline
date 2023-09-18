args <- commandArgs(T)
peakFile <- args[1]   # "raw/CD_pos.xlsx"
POS_NEG <- args[2] # "pos" / "neg"
OUTDIR <- args[3]
INFO <- args[4] # "sample_info.xls"
VSFILE <- args[5]  # vs_group.xls
THREAD <- args[6] # 2

library(metaX)
library(openxlsx)
#library(readxl)
library(stringr)

dir.create(OUTDIR, showWarnings = FALSE)


# sampleList file
sampleList <- read.table(INFO, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
sfile <- sampleList[, c("SampleName", "batch", "class", "order")]
colnames(sfile) <- c("sample", "batch", "class", "order")
write.table(sfile, file = paste0(OUTDIR, "/", "sampleList.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


# peak file
peak <- read.xlsx(peakFile, check.names=FALSE)
name_col1 <- peak$`RT.[min]`
name_col2 <- peak$`Molecular.Weight`
name_col <- paste(name_col1, name_col2, sep = "_")
samples_col_names <- grep("^Area:", colnames(peak), value = TRUE)
samples_col <- peak[, samples_col_names]
pfile <- cbind(name_col, samples_col, stringsAsFactors = FALSE)
samples_names <- str_match(samples_col_names, "Area:.(.*?).raw")[, 2]
colnames(pfile) <- c("name", sampleList[match(samples_names, sampleList$SampleID), "SampleName"])
write.table(pfile, file = paste0(OUTDIR, "/", POS_NEG, "_peak.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


# comp group
comps <- read.table(VSFILE, header = TRUE, sep = "\t", comment.char = "")
comps_groups <- paste(comps$X.treat, comps$ctrl, sep = ":")
compGroups <- paste(comps_groups, collapse = ";")

para <- new("metaXpara")
pFile <- paste0(OUTDIR, "/", POS_NEG, "_peak.txt")
sFile <- paste0(OUTDIR, "/", "sampleList.txt")
rawPeaks(para) <- read.delim(pFile, check.names = FALSE)
sampleListFile(para) <- sFile
ratioPairs(para) <- compGroups
plsdaPara <- new("plsDAPara")
plsdaPara@kfold = 5
para@outdir <- OUTDIR
prefix(para) <- POS_NEG
p <- metaXpipe(para, plsdaPara = plsdaPara, cvFilter = 0.3, 
	remveOutlier = FALSE, outTol = 1.2, doQA = FALSE, 
	doROC = FALSE, qcsc = FALSE, nor.method = "pqn", 
	pclean = FALSE, t = 1, scale = "uv", cpu = THREAD)
