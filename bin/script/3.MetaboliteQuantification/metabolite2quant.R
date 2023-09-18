suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(htmltools))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(argparser))

draw_corr_plot = function(meta_intensity, sample_class, outdir, type="neg") {

  meta_intensity <- meta_intensity[grep(type, meta_intensity$ID), ]
  
  meta_intensity$Metabolite <- lapply(meta_intensity$Metabolite, function(X){
    tryCatch({
      if(nchar(X) > 70){
        substring(X, 1, 70)
      }else{
        X
      }
    },error=function(e){
      ""
    }
    )
  })
  meta_intensity <- meta_intensity[meta_intensity$Metabolite != "", ]
  meta_intensity <- meta_intensity[!duplicated(meta_intensity$Metabolite), ]
  
  meta_intensity.metabolite <- t(meta_intensity[, (ncol(meta_intensity) - length(sample_class) + 1): ncol(meta_intensity)])



  metacorr_df = rcorr(as.matrix(meta_intensity.metabolite))$r
  rownames(metacorr_df) <- meta_intensity$Metabolite
  colnames(metacorr_df) <- meta_intensity$Metabolite

  png(
    height = 1200,
    width = 1200,
    file = paste0(outdir, "/Metabolite.correlation.png")
  )
  corrplot(
    metacorr_df,
    method = "color",
    type = "upper",
    tl.col = "black",
    tl.cex = 0.5,
    tl.srt= 45,
    order = "hclust",
    col.lim = c(-1, 1),
    col = colorRampPalette(c("red", "white", "blue"))(50)
  )
  dev.off()
  pdf(file = paste0(outdir, "/Metabolite.correlation.pdf"))
  corrplot(
    metacorr_df,
    method = "color",
    type = "upper",
    tl.col = "black",
    tl.cex = 0.2,
    tl.srt=45,
    order = "hclust",
    col.lim = c(-1, 1),
    col = colorRampPalette(c("red", "white", "blue"))(50)
  )
  dev.off()
}

draw_metabolite_intensity_plot = function(meta_intensity, sample_class, outdir, type="neg") {
  meta_intensity <- meta_intensity[grep(type, meta_intensity$ID), ]
  # make annotation col
  annotation_df_col = as.data.frame(sample_class)
  row.names(annotation_df_col) = colnames(meta_intensity)[(ncol(meta_intensity) - length(sample_class) + 1): ncol(meta_intensity)]

  # make annotation row
  data <- meta_intensity[, (ncol(meta_intensity) - length(sample_class) + 1):ncol(meta_intensity)]
  data[data==0] <- NA
  data[is.na(data)] <- min(data,na.rm = T)*0.01
  intensity_df = log10(data)
  row.names(intensity_df) = as.character(row.names(intensity_df))
  annotation_df_row = as.data.frame(meta_intensity$Superclass)
  colnames(annotation_df_row) = "MetaboliteClass"
  row.names(annotation_df_row) = row.names(intensity_df)

  pheatmap(
    intensity_df,
    border_color = NA,
    show_rownames = FALSE,
    annotation_col = annotation_df_col,
    annotation_row = annotation_df_row,
    annotation_names_col = F,
    annotation_names_row = F,
    scale="row",
    fontsize_col=10-0.5*sqrt(length(sample_class)),
    colorRampPalette(c("blue", "white", "red"))(50),
    filename = paste0(outdir, "/Metabolite.intensity.heatmap.cluster.pdf")
  )

  pheatmap(
    intensity_df,
    border_color = NA,
    show_rownames = FALSE,
    annotation_col = annotation_df_col,
    annotation_row = annotation_df_row,
    annotation_names_col = F,
    annotation_names_row = F,
    scale="row",
    fontsize_col=10-0.5*sqrt(length(sample_class)),
    colorRampPalette(c("blue", "white", "red"))(50),
    filename = paste0(outdir, "/Metabolite.intensity.heatmap.cluster.png")
  )

  pheatmap(
    intensity_df,
    border_color = NA,
    show_rownames = FALSE,
    annotation_col = annotation_df_col,
    annotation_row = annotation_df_row,
    annotation_names_col = F,
    annotation_names_row = F,
    cluster_cols = FALSE,
    scale="row",
    colorRampPalette(c("blue", "white", "red"))(50),
    filename = paste0(outdir, "/Metabolite.intensity.heatmap.pdf")
  )
  pheatmap(
    intensity_df,
    border_color = NA,
    show_rownames = FALSE,
    annotation_col = annotation_df_col,
    annotation_row = annotation_df_row,
    annotation_names_col = F,
    annotation_names_row = F,
    cluster_cols = FALSE,
    scale="row",
    colorRampPalette(c("blue", "white", "red"))(50),
    filename = paste0(outdir, "/Metabolite.intensity.heatmap.png")
  )
}

draw_single_metabolite_intensity_boxplot = function(meta_intensity,
                                                    meta_intensity_data,
                                                    id,
                                                    outdir) {
  single_data = meta_intensity[which(meta_intensity_data$ID == id), ]
  Sampleclass = c()
  for (i in colnames(meta_intensity)) {
    names = strsplit(i, ".", fixed = T)[[1]][1]
    if (length(grep("QC", names)) > 0) {
      Sampleclass = c(Sampleclass, "QC")
    } else{
      Sampleclass = c(Sampleclass, names)
    }
  }
  intensity = t(as.matrix(single_data))
  df = data.frame(Sampleclass, intensity)
  colnames(df) = c("sampleclass", "intensity")
  df = df[df["sampleclass"] != "QC", ]
  p = ggplot(data = df, aes(x = sampleclass, y = intensity)) + geom_boxplot(aes(fill =
                                                                                  sampleclass),
                                                                            show.legend = FALSE,
                                                                            size = 1) +
    geom_dotplot(
      binaxis = 'y',
      stackdir = 'center',
      dotsize = 0.6,
      fill = "pink",
      stackratio = 1
    ) +
    theme_bw() +
    xlab("Group") +
    ylab("Intensity") +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        size = 13,
        hjust = 1
      ),
      axis.text.y = element_text(size = 13),
      axis.title = element_text(size = 14, face = "bold"),
      legend.position = 'none',
    )
  if (!("boxplot" %in% dir(outdir))) {
    dir.create(paste0(outdir, "/boxplot"))
  }
  ggsave(filename = paste0(outdir, "/boxplot/", id, ".pdf"), p)
}


draw_circle_plot=function(meta_intensity, sample_class, outdir, type){
  meta_intensity <- meta_intensity[grep(type, meta_intensity$ID), ]
  meta_intensity <- meta_intensity[meta_intensity$Superclass != "", ]
  meta_intensity_data <- meta_intensity[,(ncol(meta_intensity) - length(sample_class) + 1): ncol(meta_intensity)]
  meta_intensity_data = apply(meta_intensity_data, 1, mean)
  meta_intensity_data <- log10(meta_intensity_data)

  superclass =as.factor(meta_intensity$Superclass)
  metabolite <- meta_intensity$Metabolite
  id = meta_intensity$ID
  id_meta = paste0(metabolite,"@",id)

  circle_data = data.frame(superclass,id_meta,meta_intensity_data)
  empty_bar = 10
  to_add = data.frame(matrix(NA,empty_bar*nlevels(circle_data$superclass),ncol(circle_data)))
  colnames(to_add) = colnames(circle_data)
  to_add$superclass = rep(levels(circle_data$superclass),each=empty_bar)
  circle_data=rbind(circle_data,to_add)

  circle_data$superclass = as.factor(circle_data$superclas)
  circle_data = circle_data[order(circle_data$superclass),]
  circle_data$id_meta = factor(circle_data$id_meta,levels = circle_data$id_meta)

  idx = c(1:nrow(circle_data))
  circle_data$idx = idx
  angle = -90 -360 * idx/nrow(circle_data)
  hjust = rep(1,length(idx))
  hjust[which(angle< (-270))] =  0
  angle[which(angle< (-270))] = angle[which(angle< -270)] + 180

  based_data = circle_data %>%
    group_by(superclass) %>%
    dplyr::summarize(start=min(idx),end=max(idx)-empty_bar) %>%
    rowwise() %>%
    mutate(title=mean(c(start,end)))

  p=ggplot(data=circle_data,aes(x=idx,y=meta_intensity_data,fill=superclass))+geom_bar(stat = "identity")+
    coord_polar()+
    theme_minimal()+
    geom_text(aes(x=idx,y=meta_intensity_data+2,label=id_meta),angle=angle,size=1.5,hjust=hjust,vjust=1)+
#     geom_segment(data = based_data,aes(x=start,y=-1,xend=end,yend=-1),color="#000000",fontface ="bold",size =1.5)+
    ylim(-10,20)+
    scale_fill_manual(values=colorRampPalette(brewer.pal(11, "Spectral"))(16))+
    theme(axis.title = element_blank(),
          panel.grid = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank()
    )
  ggsave(paste0(outdir,"/Metabolite.intensity.mean.circle.png"),p,width = 4800,height = 4000,units = "px")
  ggsave(paste0(outdir,"/Metabolite.intensity.mean.circle.pdf"),p,width = 4800,height = 4000,units = "px")
}

draw_mean_intensity_rank_plot = function(meta_intensity_data_all, sample_class, outdir, type) {
  meta_intensity <- meta_intensity[grep(type, meta_intensity$ID), ]
  meta_intensity <- meta_intensity[meta_intensity$Superclass != "", ]
  meta_intensity_data <- meta_intensity[,(ncol(meta_intensity) - length(sample_class) + 1): ncol(meta_intensity)]
  meta_intensity_data = apply(meta_intensity_data, 1, mean)
  meta_intensity_data <- log10(meta_intensity_data)
  Rank = c(
    rank(-meta_intensity_data, ties.method = "first")
  )
  superClass = meta_intensity$Superclass

  mean_rank_df = data.frame(meta_intensity_data, Rank, superClass)
  colnames(mean_rank_df) = c("mean_intensity_all", "Rank", "superClass")

  p = ggplot(data = mean_rank_df) + geom_point(aes(x = Rank, y = meta_intensity_data, color =
                                                     superClass)) +
    theme_bw() +
    ylab("Log10(Mean Intensity") +
    xlab("Rank") +
    theme(panel.grid = element_blank(), )
  ggsave(
    filename = paste0(outdir, "/Metabolite.intensity.mean.rank.png"),
    p,
    width = 3000,
    height = 1800,
    units = "px"
  )
  ggsave(filename = paste0(outdir, "/Metabolite.intensity.mean.rank.pdf"), p,
         width = 3000,
         height = 1800,
         units = "px")
  p_html = ggplotly(p, width = 3000, height = 1800)
  save_html(p_html,
            paste0(outdir, "/Metabolite.intensity.mean.rank.html"))
}

# main --------------------------------------------------------------------

# argv <- list(
#   intensity="../已完成数据/metaX/neg/data/neg-norm-metaboAnalystInput.csv",
#   output_dir="../已完成数据/3.MetaboliteQuantification/Metabolite2quant/neg",
#   annotation="../已完成数据/2.MetaboliteIdentification/identification.raw.intensity.csv",
#   type="neg"
# );
p <- arg_parser("section3 Feature2quant plot")
p <- add_argument(p, "--intensity", help="metaX intensity norm data")
p <- add_argument(p, "--annotation", help="section2 annotation file")
p <- add_argument(p, "--type", help="pos or neg")
p <- add_argument(p, "--output_dir", help="output dir")

argv <- parse_args(p)
# data prepare
outdir = argv$output_dir
#sample_class <- as.character(read.csv(file=argv$intensity, sep=",", skip=1, nrows=1, header=F)[1, ])
sample_class <- unlist(read.csv(file=argv$intensity, sep=",", skip=1, nrows=1, header=F), use.names=F)
sample_class <- sample_class[2:length(sample_class)]
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

# sample_list <- read.csv(file=argv$sample_file, sep="\t")
#
# annotation_data <- read.csv(argv$annotation, encoding="utf8")
# annotation_data <- annotation_data[grep(type, annotation_data$ID), ]
# # sample_cols = colnames(annotation_data)[(ncol(annotation_data) - nrow(sample_list) + 1):ncol(annotation_data)]
# intensity_data <- annotation_data[annotation_data$Metabolite != "", sample_cols]
#
# rownames(intensity_data) <- annotation_data[annotation_data$Metabolite != "", "Metabolite"]
# metacorr_df = rcorr(intensity_data)

meta_intensity = read.csv(argv$annotation, header=TRUE, stringsAsFactors=F, encoding="utf-8")
# drawplot
meta_intensity[is.na(meta_intensity)] = 0
draw_metabolite_intensity_plot(meta_intensity, sample_class, outdir, argv$type)
draw_corr_plot(meta_intensity, sample_class, outdir, argv$type)
draw_circle_plot(meta_intensity, sample_class, outdir, argv$type)
draw_mean_intensity_rank_plot(meta_intensity, sample_class, outdir, argv$type)
# draw_single_metabolite_intensity_boxplot(
#   meta_intensity = meta_intensity_neg,
#   meta_intensity_data = meta_intensity_data_neg,
#   id = "neg-M88T61",
#   outdir = outdir
# )





