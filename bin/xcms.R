args <- commandArgs(T)
SAMPLES <- args[1] # string
POS_NEG <- args[2] # "pos" / "neg"
OUTDIR <- args[3]
INFO <- args[4] # "sample_info.txt"

library(xcms)
library(scales)
library(ggplot2)
library(reshape2)

# dda_files <- c("Input/pos/LPS_1.mzXML", "Input/pos/LPS_2.mzXML", "Input/pos/LPS_3.mzXML", "Input/pos/LPS_4.mzXML", "Input/pos/LPS_5.mzXML", "Input/pos/LPS_6.mzXML", "Input/pos/LPS_Nilo_1.mzXML", "Input/pos/LPS_Nilo_2.mzXML", "Input/pos/LPS_Nilo_3.mzXML", "Input/pos/LPS_Nilo_4.mzXML", "Input/pos/LPS_Nilo_5.mzXML", "Input/pos/LPS_Nilo_6.mzXML", "Input/pos/QC1.mzXML", "Input/pos/QC2.mzXML", "Input/pos/QC3.mzXML", "Input/pos/QC4.mzXML", "Input/pos/QC5.mzXML")
#dir.create(OUTDIR,  showWarnings=FALSE, recursive=TRUE)
samples_list <- unlist(strsplit(SAMPLES, split = ','))
dda_files <- paste("Input/", POS_NEG, "/", samples_list, ".mzXML", sep = "")

info <- read.table(INFO, sep="\t", quote="", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
pd <- info[, c("SampleName", "Group")]
raw_data <- readMSData(files=dda_files, pdata=new("NAnnotatedDataFrame", pd), mode="onDisk")


# ##保留时间查看
# rtime(raw_data)

# ##查看质荷比
# mz(raw_data)

# ##查看强度
# intensity(raw_data)

##按文件分割数据
mzs <- mz(raw_data)
mzs_by_file <- split(mzs, f = fromFile(raw_data))
length(mzs_by_file)

## 离子流图信息统计总表
rt_list <- rtime(raw_data)
indensity_list <- tic(raw_data)
smp_index_list <- raw_data@featureData@data$fileIdx
smp_ids <- unique(smp_index_list)
smp_names <- raw_data$SampleName
smp_names_list_tmp <- match(smp_index_list, smp_ids)
smp_names_list <- smp_names[smp_names_list_tmp]
smp_group <- raw_data$Group
smp_group_list_tmp <- match(smp_index_list, smp_ids)
smp_group_list <- smp_group[smp_group_list_tmp]
df <- data.frame(rt_list, indensity_list, 
	             smp_names_list, smp_group_list, stringsAsFactors = FALSE)
names(df) <- c("time", "intensity", "sample", "group")


# ## 总离子流图
# bpis <- chromatogram(raw_data, aggregationFun = "sum")
# num = dim(bpis)[2]
# group_colors <- hue_pal()(num)
# names(group_colors) <- bpis@phenoData@data$SampleName
# ## Plot all chromatograms.
# pdf("Output/pos.tic.all.pdf")
# plot(bpis, col = group_colors[raw_data$SampleName])
# dev.off()

bpis <- chromatogram(raw_data, aggregationFun = "sum")
num = dim(bpis)[2]
comp <- c()
for (i in 1:num){
	comp_tmp <- names(bpis[1,i]@rtime)
	comp <- c(comp, comp_tmp)
}
df4Plot <- df[comp, ]
write.table(df4Plot, paste0(OUTDIR, "/", POS_NEG, ".tic.xls"), quote = FALSE, sep = "\t", 
			col.names = TRUE, row.names = FALSE)

P <- ggplot(df4Plot, aes(x = time, y = intensity, color = sample, group = sample)) +
  geom_line() + 
  theme_classic() +
  ylab("Intensity") + xlab("Retention time") +
  scale_x_continuous(breaks = seq(0, 600, 150)) +
  scale_y_continuous(breaks = seq(0, 3 * (10**9), 10**9))
ggsave(filename = paste0(OUTDIR, "/", POS_NEG, ".tic.all.png"), P, width = 30, height = 21, units = "cm", dpi = 600)
ggsave(filename = paste0(OUTDIR, "/", POS_NEG, ".tic.all.pdf"), P, width = 30, height = 21, units = "cm")


P <- ggplot(df4Plot, aes(x = time, y = intensity, color = group, group = group)) +
  geom_line() + 
  facet_grid(group ~ .) + 
  theme_bw() +
  ylab("Intensity") + xlab("Retention time") +
  scale_x_continuous(breaks = seq(0, 600, 150)) +
  scale_y_continuous(breaks = seq(0, 3 * (10**9), 10**9))
ggsave(filename = paste0(OUTDIR, "/", POS_NEG, ".tic.all.facet.png"), P, width = 30, height = 21, units = "cm", dpi = 600)
ggsave(filename = paste0(OUTDIR, "/", POS_NEG, ".tic.all.facet.pdf"), P, width = 30, height = 21, units = "cm")


for (i in unique(df4Plot$group)){
  P <- ggplot(df4Plot[which(df4Plot$group == i), ], aes(x = time, y = intensity, color = sample, group = sample)) +
	geom_line() + 
	theme_classic() +
    ylab("Intensity") + xlab("Retention time") +
    scale_x_continuous(breaks = seq(0, 600, 150)) +
    scale_y_continuous(breaks = seq(0, 3 * (10**9), 10**9))
  # ggsave(filename = paste0("Output/pos.tic.class", i, ".png"), P, width = 30, height = 21, units = "cm", dpi = 600)
  ggsave(filename = paste0(OUTDIR, "/", POS_NEG, ".tic.class", i, ".pdf"), P, width = 30, height = 21, units = "cm")
}

## 单个样本的离子流图
for (i in 1:num){
	smp <- bpis@phenoData@data$SampleName[i]
	group_colors <- hue_pal()(1)
	names(group_colors) <- smp
	pdf(paste0(OUTDIR, "/", POS_NEG, ".tic.eachsample.", smp, ".pdf"))
	plot(bpis, col = group_colors[raw_data$SampleName])
	dev.off()
}


## centWave进行峰检测
cwp <- CentWaveParam(ppm=30, peakwidth=c(5, 25), snthresh=6, mzdiff=0.01)
xdata <- findChromPeaks(raw_data, param = cwp)
write.table(chromPeaks(xdata), file = paste0(OUTDIR, "/", POS_NEG, ".peak_raw.info.xls"), row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

# 统计检测到的峰
summary_fun <- function(z) {
	c(peak_count = nrow(z), rt = quantile(z[, "rtmax"] - z[, "rtmin"]))
}
T <- lapply(split.data.frame(chromPeaks(xdata), 
	        f = chromPeaks(xdata)[, "sample"]), 
	        FUN = summary_fun)
T <- do.call(rbind, T)
rownames(T) <- basename(fileNames(xdata))
# pandoc.table(T, 
# 	         caption = paste0("Summary statistics on identified chromatographic",
#              " peaks. Shown are number of identified peaks per",
#              " sample and widths/duration of chromatographic ",
#              "peaks."))
# write.table(T, file = "Output/pos.peak_raw.stat.xls", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

# # 保留时间校正
# xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = 0.1))
# num <- length(unique(xdata$Group))
# group_colors <- hue_pal()(num)
# names(group_colors) <- unique(xdata$Group)
# png("Output/pos-aligned_rt.png")
# plotAdjustedRtime(xdata, col = group_colors[xdata$Group])
# dev.off()
# #xdata <- dropAdjustedRtime(xdata)

## 按照组进行保留时间校正
pdp <- PeakDensityParam(sampleGroups = xdata$Group, minFraction = 0.8)
#根据组来提取数据
xdata <- groupChromPeaks(xdata, param = pdp)
#根据组来校正时间
pgp <- PeakGroupsParam(minFraction = 0.85)
xdata <- adjustRtime(xdata, param = pgp)
#对校正前和校正后的结果作图
num <- length(unique(xdata$Group))
group_colors <- hue_pal()(num)
names(group_colors) <- unique(xdata$Group)
png(paste0(OUTDIR, "/", POS_NEG, "-aligned_rt.png"))
plotAdjustedRtime(xdata, col = group_colors[xdata$Group],peakGroupsCol = "grey", peakGroupsPch = 1)
dev.off()

## 按照组进行峰提取
pdp <- PeakDensityParam(sampleGroups = xdata$Group,
                        minFraction = 0.5, bw = 5)
xdata <- groupChromPeaks(xdata, param = pdp)


## 峰补齐
xdata <- fillChromPeaks(xdata)


## 保存结果
# 获取峰检测得到的feature信息
feature.info <- featureDefinitions(xdata)
# 获取feature定量信息（峰面积）
feature.quant <- featureValues(xdata, value = "into")
# 数据导出
result <- cbind(as.data.frame(feature.info), feature.quant)
feature <- rownames(result)
out <- cbind(feature, result)
colnames(out) <- c("feature", colnames(result))
write.table(out, file = paste0(OUTDIR, "/", POS_NEG, "_result.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


# png("Output/pos-mz-rt.png")
# plot(xdata, type = "XIC")
# dev.off()


library(LSD)
png(paste0(OUTDIR, "/", POS_NEG, "-mz-rt.png"))
heatscatter(result$rtmed,  result$mzmed, 
	colpal = "bl2gr2rd", cor = FALSE,
	main = "",
	xlab = "Average retention time", ylab = "Average m/z")
dev.off()

