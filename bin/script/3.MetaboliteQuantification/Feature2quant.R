suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(argparser))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(htmltools))


draw_boxplot = function(intensity_df,outdir){
  intensity_box_plot_data = intensity_df[c(2:length(row.names(intensity_df))),c(2:length(colnames(intensity_df)))]
  Log10_intensity = log10(as.double(as.vector(t(as.matrix(intensity_box_plot_data)))))
  samples = rep(colnames(intensity_box_plot_data),each = length(row.names(intensity_box_plot_data)))
  class = rep(as.matrix(intensity_df[1,c(2:length(colnames(intensity_df)))]),each=length(row.names(intensity_box_plot_data)))
  intensity_boxplot_df = data.frame(Log10_intensity = Log10_intensity,samples=samples,class=class)
  
  
  p=ggplot(data=intensity_boxplot_df,aes(samples,Log10_intensity,fill=class))+
    geom_boxplot()+
    theme_bw()+
    ylab("Log10(intensity)")+
    xlab("Sample")+
    theme(
      axis.text.x = element_blank(), # Remove x axis tick labels
      axis.ticks = element_blank(),
      axis.title=element_text(size=14,face="bold"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold"))
  ggsave(filename = paste0(outdir,"/intensity.boxplot.png"),p)
  ggsave(filename = paste0(outdir,"/intensity.boxplot.pdf"),p)
}

draw_heatmap_plot=function(intensity_df,outdir){
  intensity_hm_df_original = intensity_df[c(2:length(row.names(intensity_df))),c(2:length(colnames(intensity_df)))]
  
  
  intensity_hm_matrix = as.matrix(intensity_hm_df_original)
  intensity_hm_matrix_num = apply(intensity_hm_matrix, 2, as.numeric)
  intensity_hm_matrix_num = apply(intensity_hm_matrix_num, 2, log10)
  intensity_hm_df = as.data.frame(intensity_hm_matrix_num)
  row.names(intensity_hm_df) = intensity_df[c(2:length(row.names(intensity_df))),1]
  colnames(intensity_hm_df) = colnames(intensity_hm_df_original)
  Sampleclass= c(as.matrix(intensity_df[1,c(2:length(colnames(intensity_df)))]))
  
  samples = colnames(intensity_hm_df)

  deta_samples = length(samples)*0.3
  if (deta_samples<12) {
    deta_samples=12
  }
  annotation_df = data.frame(Sampleclass)
  rownames(annotation_df) <- colnames(intensity_hm_df)
  
  pheatmap(intensity_hm_df,border_color=NA,show_rownames = FALSE,
           annotation_col =annotation_df,colorRampPalette(c("blue", "white", "red"))(50),
           width = deta_samples,
           height = 15,
           scale = "row",
           filename = paste0(outdir,"/intensity.heatmap.cluster.pdf"))
  pheatmap(intensity_hm_df,border_color=NA,show_rownames = FALSE,
           annotation_col =annotation_df,colorRampPalette(c("blue", "white", "red"))(50),
           width = deta_samples,
           height = 15,
           scale = "row",
           filename = paste0(outdir,"/intensity.heatmap.cluster.png"))
  
  pheatmap(intensity_hm_df,border_color=NA,show_rownames = FALSE,
           annotation_col =annotation_df,cluster_cols = FALSE,colorRampPalette(c("blue", "white", "red"))(50),
           width = deta_samples,
           height = 15,
           scale = "row",
           filename = paste0(outdir,"/intensity.heatmap.pdf"))
  pheatmap(intensity_hm_df,border_color=NA,show_rownames = FALSE,
           annotation_col =annotation_df,cluster_cols = FALSE,colorRampPalette(c("blue", "white", "red"))(50),
           width = deta_samples,
           height = 15,
           scale = "row",
           filename = paste0(outdir,"/intensity.heatmap.png"))
}

draw_corrplot=function(intensity_df,outdir){
  intensity_df_original = intensity_df[c(2:length(row.names(intensity_df))),c(2:length(colnames(intensity_df)))]
  
  
  intensity_matrix = as.matrix(intensity_df_original)
  intensity_matrix_num = apply(intensity_matrix, 2, as.numeric)
  res = rcorr(intensity_matrix_num)
  png(height=1980, width=1980, file=paste0(outdir,"/sample.correlation.png"))
  corrplot::corrplot.mixed(res$r,lower = "ellipse",upper = "number",
                           order = 'alphabet',tl.pos = "lt",tl.cex =2,
                           tl.col="black",col.lim=c(0,1),number.cex=2,)
  dev.off()
  pdf(file=paste0(outdir,"/sample.correlation.pdf"))
  corrplot::corrplot.mixed(res$r,lower = "ellipse",upper = "number",
                           order = 'alphabet',tl.pos = "lt",tl.cex =0.6,
                           tl.col="black",col.lim=c(0,1),number.cex=0.6,)
  dev.off()
}

draw_pca_loading_plot = function(pca_loading_file,outdir){
  
  colourCount = length(unique(pca_loading_file$Superclass))
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  
  p = ggplot(data = pca_loading_file,aes(x=PC1,y=PC2,colour=Superclass))+
    geom_point()+
    theme_bw()+
    scale_color_manual(values=getPalette(colourCount))+
    theme(legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold",size=14),
          axis.title = element_text(face = "bold"),
          panel.grid=element_blank(),
    )
  ggsave(filename = paste0(outdir,"/pca.loading.png"),p,width = 2100,height = 1200,units = "px")
  ggsave(filename = paste0(outdir,"/pca.loading.pdf"),p,width = 2100,height = 1200,units = "px")
  p_html = ggplotly(p,width = 2100,height = 1200)
  save_html(p_html,paste0(outdir,"/pca.loading.html"))
}

draw_pca_score_plot = function(pca_score_file,outdir){
  pca_score_noqc = pca_score_file[pca_score_file["class"] != "QC",]
  
  colourCount = length(unique(pca_score_file$class))
  colourCount_noqc = length(unique(pca_score_file$class))
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  
  p_2d = ggplot(data = pca_score_file,aes(x=pca_score_file[,2],y=pca_score_file[,3],colour=class,label=sample))+
    geom_point()+
    theme_bw()+
    scale_color_manual(values=getPalette(colourCount))+
    stat_ellipse(level = 0.95)+
    xlab("PC1")+
    ylab("PC2")+
    #geom_text(aes(x=pca_score_file[,2],y=pca_score_file[,3]))+ #是否显示样本名
    theme(legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold",size=14),
          axis.title = element_text(face = "bold"),
          panel.grid=element_blank(),
    )
  
  p_2d_no_qc = ggplot(data = pca_score_noqc,aes(x=pca_score_noqc[,2],y=pca_score_noqc[,3],colour=class,label=sample))+
    geom_point()+
    theme_bw()+
    scale_color_manual(values=getPalette(colourCount_noqc))+
    stat_ellipse(level = 0.95)+
    xlab("PC1")+
    ylab("PC2")+
    #geom_text(aes(x=pca_score_file[,2],y=pca_score_file[,3]))+ #是否显示样本名
    theme(legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold",size=14),
          axis.title = element_text(face = "bold"),
          panel.grid=element_blank(),
    )
  
  ggsave(filename = paste0(outdir,"/pca.score.png"),p_2d,width = 1500,height = 1200,units = "px")
  ggsave(filename = paste0(outdir,"/pca.score.pdf"),p_2d,width = 1500,height = 1200,units = "px")
  ggsave(filename = paste0(outdir,"/pca.score_no_qc.png"),p_2d_no_qc,width = 1500,height = 1200,units = "px")
  ggsave(filename = paste0(outdir,"/pca.score_no_qc.pdf"),p_2d_no_qc,width = 1500,height = 1200,units = "px")
  p_html = ggplotly(p_2d,width = 1500,height = 1200)
  save_html(p_html,paste0(outdir,"/pca.score.2d.html"))
  
  
  # 3D plot
  p_html_3d <- plot_ly(data = pca_score, x = pca_score[,2], y =  pca_score[,3], z = pca_score[,4], 
                       color = ~class, colors = getPalette(colourCount))
  p_html_3d <- p_html_3d %>% add_markers()
  p_html_3d <- p_html_3d %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                                 yaxis = list(title = 'PC2'),
                                                 zaxis = list(title = 'PC3')))
  save_html(p_html_3d,paste0(outdir,"/pca.score.3d.html"))
}

############# main ###############
# data prepare
# argv <- list(
#   intensity="Output/metaX/neg/data/neg-norm-metaboAnalystInput.csv",
#   pca_rds="Output/metaX/neg/data/neg-pca.rds",
#   sample_file="Output/metaX/neg/sampleList.txt",
#   output_dir="Output/3.MetaboliteQuantification/Feautre2quant/neg",
#   annotation="Output/2.MetaboliteIdentification/identification.raw.intensity.csv",
#   type="neg"
# );
p <- arg_parser("section3 Feature2quant plot")
p <- add_argument(p, "--intensity", help="metaX intensity norm data")
p <- add_argument(p, "--pca_rds", help="pca rds file")
p <- add_argument(p, "--sample_file", help="sample_file")
p <- add_argument(p, "--annotation", help="section2 annotation file")
p <- add_argument(p, "--type", help="pos or neg")
p <- add_argument(p, "--output_dir", help="output dir")

argv <- parse_args(p)

outdir <- argv$output_dir
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
intensity_df = read.csv(file = argv$intensity,sep = ",")
pca_object <- readRDS(argv$pca_rds)

sample_list <- read.csv(file=argv$sample_file, sep="\t", stringsAsFactors=FALSE)
sample_list[is.na(sample_list$class), 'class'] <- "QC"
rownames(sample_list) <- sample_list$sample
# pca score
pca_score = as.data.frame(pca_object@scores)
pca_score <- merge(pca_score, sample_list, by=0)
# loading
type <- argv$type
annotation_data <- read.csv(argv$annotation, stringsAsFactors=FALSE)
annotation_data <- annotation_data[grep(type, annotation_data$ID), ]
write.csv(annotation_data, file=paste0(outdir, "/intensity.csv") ,row.names=F, sep=",")

annotation_data <- annotation_data[, c("ID", "MZ", "RT", "Metabolite", "Superclass")]
annotation_data$Sample <- lapply(annotation_data$ID, function(X){strsplit(X, "-", fixed=T)[[1]][2]})
rownames(annotation_data) <- annotation_data$Sample
annotation_data <- annotation_data[order(rownames(annotation_data)), ]
pca_loading = as.data.frame(pca_object@loadings)
pca_loading = merge(pca_loading, annotation_data, by=0)
pca_loading <- pca_loading[pca_loading$Superclass !="", ]
# draw plot
#draw_boxplot(intensity_df = intensity_df,outdir = outdir)
draw_heatmap_plot(intensity_df = intensity_df,outdir=outdir)
draw_corrplot(intensity_df,outdir)
draw_pca_score_plot(pca_score,outdir)
draw_pca_loading_plot(pca_loading,outdir)
