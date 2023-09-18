suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(htmltools))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(argparser))
# options(warn = -1)
# stat plot
rscript_bin = " /opt/anaconda3/envs/proteomics/bin/Rscript"
draw_stat_data = function(data, outdir) {
  # main plot
  for(col in colnames(data)){
    if(col == "comparison"){

    }else{
      data[, col] = as.numeric(data[, col])
    }
  }
  Comparison = rep(data$comparison, 2)
  num_feature = c(data$pos_up + data$neg_up, data$pos_down + data$neg_down)
  group = c(rep("up", length(data$comparison)), rep("down", length(data$comparison)))
  p = ggplot() + geom_bar(
    aes(x = Comparison, y = num_feature, fill = group,),
    stat = "identity",
    position = 'dodge'
  ) +
    scale_fill_manual(values = colorRampPalette(c("#2c7fb8", "#f1a340"))(2)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(
            angle = 45,
            size = 10,
            vjust = 0.5
          ))
  ggsave(paste0(outdir, "/comparison_stat.png"), p)
  
  
  # facet plot
  Comparison_face = rep(data$comparison, 4)
  group_face =  c(rep("up", length(Comparison)), rep("down", length(Comparison)))
  num_feature_face = c(data$pos_up, data$neg_up, data$pos_down, data$neg_down)
  face_pos_neg = c(
    rep("positive", length(data$comparison)),
    rep("negative", length(data$comparison)),
    rep("positive", length(data$comparison)),
    rep("negative", length(data$comparison))
  )
  
  data = data.frame(Comparison_face, num_feature_face, group_face, face_pos_neg)
  p1 = ggplot(data = data) + geom_bar(
    aes(x = Comparison_face, y = num_feature_face, fill = group_face),
    stat = "identity",
    position = 'dodge'
  ) +
    theme_bw() +
    scale_fill_manual(values = colorRampPalette(c("#2c7fb8", "#f1a340"))(2)) +
    facet_wrap( ~ face_pos_neg, nrow = 2, strip.position = "right") +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(
            angle = 45,
            size = 10,
            vjust = 0.5
          ))
  ggsave(paste0(outdir, "/comparison_stat_facet.png"), p1)
}

draw_sig_metabolite_heatmap = function(data,outdir,name, sample_class, target_sample_cols, pvalue_cutoff=0.05, log2ratio_cutoff=1, pvalue_type='t.test_p.value', vip=1){
  data.metabolite = data[which(data$Metabolite!=""),]
  data.metabolite = data.metabolite[which(data.metabolite[pvalue_type]<pvalue_cutoff & abs(log2(data.metabolite$ratio))>log2ratio_cutoff & data.metabolite$VIP > vip),]
  data.metabolite <- data.metabolite[!duplicated(data.metabolite$Metabolite), ]
  meta_intensity = data.metabolite[, target_sample_cols]
  # make annotation col
  annotation_df_col = as.data.frame(sample_class)

  row.names(annotation_df_col) = colnames(meta_intensity)
  
  # make annotation row
  meta_intensity_numeric = apply(meta_intensity, 2, as.numeric)
  if(NROW(meta_intensity) == 1){
    intensity_df = data.frame(t(meta_intensity_numeric))
  }else{
      intensity_df = as.data.frame(meta_intensity_numeric)
  }
  row.names(intensity_df) = data.metabolite$Metabolite
  # row.names(intensity_df) = lapply(data.metabolite$Metabolite, function(X){
  #   if(nchar(X) > 50){
  #     substring(X, 1, 50)
  #   }else{
  #     X
  #   }
  # })
  if (dim(intensity_df)[1] > 1){
      deta_row = nrow(intensity_df) *0.1
      if (deta_row<8) {
        deta_row=8
      }
      if(deta_row > 100){
        deta_row=100
      }
      deta_col = ncol(intensity_df) * 0.4
      if (deta_col<12) {
        deta_col=12
      }


      pheatmap(
        intensity_df,
        border_color = "#666666",
        show_rownames = T,
        annotation_col = annotation_df_col,
        annotation_names_col = F,
        annotation_names_row = F,
        fontsize_row=5,
        cluster_cols = T,
        scale = "row",
        colorRampPalette(c("blue", "white", "red"))(50),
        filename = paste0(outdir, "/", name, ".heatmap_Metabolite.png"),
        width = deta_col,
        height = deta_row
      )

      pheatmap(
        intensity_df,
        border_color = "#666666",
        show_rownames = T,
        annotation_col = annotation_df_col,
        annotation_names_col = F,
        annotation_names_row = F,
        cluster_cols = T,
        fontsize_row=5,
        scale = "row",
        colorRampPalette(c("blue", "white", "red"))(50),
        filename = paste0(outdir, "/", name, ".heatmap_Metabolite.pdf"),
        width = deta_col,
        height = deta_row
      )

      pheatmap(
        intensity_df,
        border_color = "#666666",
        show_rownames = T,
        annotation_col = annotation_df_col,
        annotation_names_col = F,
        annotation_names_row = F,
        cluster_cols = F,
        fontsize_row=5,
        scale = "row",
        colorRampPalette(c("blue", "white", "red"))(50),
        filename = paste0(outdir,"/", name,".heatmap_Metabolite_nocluster.pdf"),
        width = deta_col,
        height = deta_row
      )
      pheatmap(
        intensity_df,
        border_color = "#666666",
        show_rownames = T,
        annotation_col = annotation_df_col,
        annotation_names_col = F,
        annotation_names_row = F,
        cluster_cols = F,
        fontsize_row=5,
        scale = "row",
        colorRampPalette(c("blue", "white", "red"))(50),
        filename = paste0(outdir, "/", name,".heatmap_Metabolite_nocluster.png"),
        width = deta_col,
        height = deta_row
      )
  }else{
    print("no enough significant data found")
    cat(NULL, file=paste0(outdir, "/", name, ".heatmap_Metabolite.png"))
    cat(NULL, file=paste0(outdir, "/", name, ".heatmap_Metabolite.pdf"))
    cat(NULL, file=paste0(outdir,"/", name,".heatmap_Metabolite_nocluster.pdf"))
    cat(NULL, file=paste0(outdir, "/", name,".heatmap_Metabolite_nocluster.png"))
  }

}

draw_sig_heatmap = function(data,outdir,name,sample_class,target_sample_cols, pvalue_cutoff=0.05, log2ratio_cutoff=1, pvalue_type='t.test_p.value', vip=1){

  data = data[which(data[pvalue_type]<pvalue_cutoff & abs(log2(data$ratio))>log2ratio_cutoff & data$VIP > vip),]
  
  meta_intensity = data[, target_sample_cols]
  # prepare sample class
  
  # make annotation col
  annotation_df_col = as.data.frame(sample_class)
  row.names(annotation_df_col) = colnames(meta_intensity)
  
  # make annotation row
  meta_intensity_numeric = apply(meta_intensity, 2, as.numeric)
  if(NROW(meta_intensity) == 1){
      intensity_df = data.frame(t(meta_intensity_numeric))
  }else{
      intensity_df = as.data.frame(meta_intensity_numeric)
  }

  if(dim(intensity_df)[1] > 1){
      deta_col = ncol(intensity_df) * 0.4
      if (deta_col<12) {
        deta_col=12
      }

      pheatmap(
        intensity_df,
        border_color = "#666666",
        show_rownames = F,
        annotation_col = annotation_df_col,
        annotation_names_col = F,
        annotation_names_row = F,
        cluster_cols = T,
        scale = "row",
        colorRampPalette(c("blue", "white", "red"))(50),
        filename = paste0(outdir, "/", name,".heatmap.pdf"),
        width = deta_col,
        height = 12
      )

      pheatmap(
        intensity_df,
        border_color = "#666666",
        show_rownames = F,
        annotation_col = annotation_df_col,
        annotation_names_col = F,
        annotation_names_row = F,
        cluster_cols = T,
        scale = "row",
        colorRampPalette(c("blue", "white", "red"))(50),
        filename = paste0(outdir, "/", name, ".heatmap.png"),
        width = deta_col,
        height = 12
      )

      pheatmap(
        intensity_df,
        border_color = "#666666",
        show_rownames = F,
        annotation_col = annotation_df_col,
        annotation_names_col = F,
        annotation_names_row = F,
        cluster_cols = F,
        scale = "row",
        colorRampPalette(c("blue", "white", "red"))(50),
        filename = paste0(outdir, "/", name,"_intensity.heatmap.pdf"),
        width = deta_col,
        height = 12
      )
      pheatmap(
        intensity_df,
        border_color = "#666666",
        show_rownames = F,
        annotation_col = annotation_df_col,
        annotation_names_col = F,
        annotation_names_row = F,
        cluster_cols = F,
        scale = "row",
        colorRampPalette(c("blue", "white", "red"))(50),
        filename = paste0(outdir, "/", name,"_intensity.heatmap.png"),
        width = deta_col,
        height = 12
      )
  }else{
    print("no enough significant data found")
    cat(NULL, file=paste0(outdir, "/", name,".heatmap.pdf"))
    cat(NULL, file=paste0(outdir, "/", name, ".heatmap.png"))
    cat(NULL, file=paste0(outdir, "/", name,"_intensity.heatmap.pdf"))
    cat(NULL, file=paste0(outdir, "/", name,"_intensity.heatmap.png"))
  }
}

draw_volcano = function(data,outdir,name, pvalue_cutoff=0.05, log2ratio_cutoff=1, pvalue_type='t.test_p.value', vip=1){
  cut_y=-log10(pvalue_cutoff)
  cut_x = log2ratio_cutoff
  FC = log2(data$ratio)
  log_Q_value = c(as.matrix(-log10(data[pvalue_type])))
  
  
  DEM = data.frame(FC,log_Q_value)
  regulated = rep("none",length(FC))
  regulated[which(FC>cut_x & log_Q_value>cut_y & data$VIP > vip)]="up"
  regulated[which(FC< -cut_x & log_Q_value>cut_y & data$VIP > vip)]="down"
  
  p = ggplot(data=DEM, aes(x=FC, y=log_Q_value,color=regulated)) + geom_point(alpha=1, size=1) +
    geom_hline(yintercept = cut_y, linetype = 4) +
    geom_hline(yintercept = cut_y, linetype = 4) +
    geom_vline(xintercept = cut_x, linetype = 4) +
    geom_vline(xintercept = -cut_x, linetype = 4) +
    theme_bw()+
    theme(
          panel.grid = element_blank(),
          legend.title = element_text(size =5),
          legend.text = element_text(size =5)
          ) + 
    xlab(paste0("Log2",name))+
    ylab("-log10 P-value") + 
    xlim(-10,10)+
    scale_colour_manual(values = c('#008B00','#999999','red')) 	## corresponding to the levels(res$change)
  
  ggsave(paste0(outdir, "/", name, ".volcanoplot.png"),p,width = 1800,height = 1200,units = "px")
  ggsave(paste0(outdir, "/", name, ".volcanoplot.pdf"),p,width = 1800,height = 1200,units = "px")
  
}

draw_volcano_sig_meta = function(data,outdir,name, pvalue_cutoff=0.05, log2ratio_cutoff=1, pvalue_type='t.test_p.value', vip=1){
  data=data[which(data$Superclass != ""),]
  data = data[which(data[pvalue_type]<pvalue_cutoff & abs(log2(data$ratio))>log2ratio_cutoff & data$VIP > vip),]
  
  cut_y=-log10(pvalue_cutoff)
  cut_x = log2ratio_cutoff
  FC = log2(data$ratio)
  log_Q_value = c(as.matrix(-log10(data[pvalue_type])))
  superClass = data$Superclass
  
  DEM = data.frame(FC,log_Q_value,superClass)
  # regulated = rep("none",length(FC))
  # regulated[which(FC>cut_x & log_Q_value>cut_y)]="up"
  # regulated[which(FC< -cut_x & log_Q_value>cut_y)]="down"
  
  p = ggplot(data=DEM, aes(x=FC, y=log_Q_value,color=superClass,shape=superClass)) + geom_point(alpha=1, size=2) +
    geom_hline(yintercept = cut_y, linetype = 4) +
    geom_hline(yintercept = cut_y, linetype = 4) +
    geom_vline(xintercept = cut_x, linetype = 4) +
    geom_vline(xintercept = -cut_x, linetype = 4) +
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      legend.title = element_text(size =5),
      legend.text = element_text(size =5)
    ) + 
    xlab(paste0("Log2",name))+
    ylab("-log10 P-value") + 
    xlim(-10,10)+
    scale_colour_manual(values = colorRampPalette(brewer.pal(8,"Spectral"))(50))	## corresponding to the levels(res$change)
  ggsave(paste0(outdir,"/", name, ".volcanoplot.Metabolite.png"),p,width = 1800,height = 1200,units = "px")
  ggsave(paste0(outdir,"/", name, ".volcanoplot.Metabolite.pdf"),p,width = 1800,height = 1200,units = "px")
}

draw_plsad_loading = function(data_merge,outdir,name, ion_type){
  VIP_type = rep(0.5,length(data_merge$ID))
  VIP_type[which(data_merge$VIP>1)] = 1.5
  data_merge$VIP_type = as.factor(VIP_type)
  
  
  p = ggplot(data = data_merge,aes(x=PC1,y=PC2,colour=VIP,size = VIP_type))+
    geom_point()+
    theme_bw()+
    scale_color_gradient(low = "blue",high = "red")+
    scale_size_manual(values = c(0.5,1.5),labels=c("<1",">1"),name="VIP type")+
    theme(
          legend.title = element_text(face = "bold",size=10),
          axis.title = element_text(face = "bold"),
          panel.grid=element_blank(),
    )
  ggsave(filename = paste0(outdir,"/", name,".", ion_type,".plsda.loading.png"),p,width = 1500,height = 1200,units = "px")
  ggsave(filename = paste0(outdir,"/", name,".", ion_type,".plsda.loading.pdf"),p,width = 1500,height = 1200,units = "px")
}

draw_pca_loading_plot = function(pca_loading_file,outdir, name, ion_type){
  pca_loading_file = pca_loading_file[pca_loading_file$Superclass != "", ]
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
  ggsave(filename = paste0(outdir, "/", name ,".", ion_type,".pca.loading.png"),p,width = 2100,height = 1200,units = "px")
  ggsave(filename = paste0(outdir,"/", name ,".", ion_type,".pca.loading.pdf"),p,width = 2100,height = 1200,units = "px")
  p_html = ggplotly(p,width = 2100,height = 1200)
  save_html(p_html,paste0(outdir,"/", name ,".", ion_type,".pca.loading.html"))
}

draw_pca_score_plot = function(pca_score_file,outdir, name, ion_type){
  pca_score_noqc = pca_score_file[which(pca_score_file["class"]!="QC"),]
  
  colourCount = length(unique(pca_score_file$class))
  colourCount_noqc = length(unique(pca_score_file$class))
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  
  p_2d = ggplot(data = pca_score_file,aes(x=pca_score_file[,2],y=pca_score_file[,3],colour=class,label=sample))+
    geom_point()+
    theme_bw()+
    scale_color_manual(values=getPalette(colourCount))+
    stat_ellipse(level = 0.9)+
    xlab("PC1")+
    ylab("PC2")+
    theme(legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold",size=14),
          axis.title = element_text(face = "bold"),
          panel.grid=element_blank(),
    )
  
  p_2d_no_qc = ggplot(data = pca_score_noqc,aes(x=pca_score_noqc[,2],y=pca_score_noqc[,3],colour=class,label=sample))+
    geom_point()+
    theme_bw()+
    scale_color_manual(values=getPalette(colourCount_noqc))+
    stat_ellipse(level = 0.9)+
    xlab("PC1")+
    ylab("PC2")+
    theme(legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold",size=14),
          axis.title = element_text(face = "bold"),
          panel.grid=element_blank(),
    )
  
  ggsave(filename = paste0(outdir,"/", name ,".", ion_type,".pca.score.png"),p_2d,width = 1500,height = 1200,units = "px")
  ggsave(filename = paste0(outdir,"/", name ,".", ion_type,".pca.score.pdf"),p_2d,width = 1500,height = 1200,units = "px")
  ggsave(filename = paste0(outdir,"/", name ,".", ion_type,".pca.score_no_qc.png"),p_2d_no_qc,width = 1500,height = 1200,units = "px")
  ggsave(filename = paste0(outdir,"/", name ,".", ion_type,".pca.score_no_qc.pdf"),p_2d_no_qc,width = 1500,height = 1200,units = "px")
  p_html = ggplotly(p_2d,width = 1500,height = 1200)
  save_html(p_html,paste0(outdir,"/", name ,".", ion_type,".pca.score.2d.html"))
  
  
  # 3D plot
  p_html_3d <- plot_ly(data = pca_score_file, x = pca_score_file[,2], y =  pca_score_file[,3], z = pca_score_file[,4], 
                       color = ~class, colors = getPalette(colourCount))
  p_html_3d <- p_html_3d %>% add_markers()
  p_html_3d <- p_html_3d %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                                 yaxis = list(title = 'PC2'),
                                                 zaxis = list(title = 'PC3')))
  save_html(p_html_3d,paste0(outdir,"/", name ,".", ion_type,".pca.score.3d.html"))
}

draw_plsda_validation_plot = function(result, out_file){
  plotdat <- as.matrix(cbind(unlist(result$plsda$res),result$plsda$perm))
    plotdat <- as.data.frame(t(plotdat))
    x1 <- plotdat$cor[1]
    y1 <- plotdat$R2[1]
    y2 <- plotdat$Q2[1]
    plotdat$cor <- abs(plotdat$cor)
    plotdat <- plotdat[order(plotdat$cor),]
    png(file = paste0(out_file, ".png"))
    par(mar=c(3,3,2,1),mgp=c(1.6,0.6,0),cex.lab=1.2,cex.main=0.9)
    plot(plotdat$cor,plotdat$R2,ylim=c(min(plotdat$R2,plotdat$Q2),1),pch=16,
         xlab="Cor",ylab="Value",col="blue")
    points(plotdat$cor,plotdat$Q2,pch=15,col="red")

    lm.r <- lm(I(R2-y1)~I(cor-x1)+0,data=plotdat)
    lm.q <- lm(I(Q2-y2)~I(cor-x1)+0,data=plotdat)
    #lines(plotdat$cor,predict(lm.r,data=plotdat$cor),col="blue",lty=2)
    #lines(plotdat$cor,predict(lm.q,data=plotdat$cor),col="red",lty=2)
    int.R <- predict(lm.r,newdata=list(cor=0))+y1
    int.Q <- predict(lm.q,newdata=list(cor=0))+y2

    abline(int.R,coef(lm.r),lty=2,col="blue")
    abline(int.Q,coef(lm.q),lty=2,col="red")
    #abline(lm.q,lty=2,col="red")
    legend("bottomright",pch=c(16,15),legend = c("R2","Q2"),col=c("blue","red"))

    title(main = paste("Intercepts:","R2=(0.0,",sprintf("%.4f",int.R),
                       "), Q2=(0.0,",sprintf("%.4f",int.Q),")"))
    dev.off()

    pdf(file = paste0(out_file, ".pdf"))
    par(mar=c(3,3,2,1),mgp=c(1.6,0.6,0),cex.lab=1.2,cex.main=0.9)
    plot(plotdat$cor,plotdat$R2,ylim=c(min(plotdat$R2,plotdat$Q2),1),pch=16,
         xlab="Cor",ylab="Value",col="blue")
    points(plotdat$cor,plotdat$Q2,pch=15,col="red")

    lm.r <- lm(I(R2-y1)~I(cor-x1)+0,data=plotdat)
    lm.q <- lm(I(Q2-y2)~I(cor-x1)+0,data=plotdat)
    #lines(plotdat$cor,predict(lm.r,data=plotdat$cor),col="blue",lty=2)
    #lines(plotdat$cor,predict(lm.q,data=plotdat$cor),col="red",lty=2)
    int.R <- predict(lm.r,newdata=list(cor=0))+y1
    int.Q <- predict(lm.q,newdata=list(cor=0))+y2

    abline(int.R,coef(lm.r),lty=2,col="blue")
    abline(int.Q,coef(lm.q),lty=2,col="red")
    #abline(lm.q,lty=2,col="red")
    legend("bottomright",pch=c(16,15),legend = c("R2","Q2"),col=c("blue","red"))

    title(main = paste("Intercepts:","R2=(0.0,",sprintf("%.4f",int.R),
                       "), Q2=(0.0,",sprintf("%.4f",int.Q),")"))
    dev.off()
}

##############----main-----##################
argv <- list(
  quant_neg="Output/metaX/neg/data/neg-quant.txt",
  quant_pos="Output/metaX/pos/data/pos-quant.txt",
  annotation="Output/2.MetaboliteIdentification/identification.raw.intensity.csv",
  intensity="Output/metaX/neg/data/neg-norm-metaboAnalystInput.csv",
  sample_file="Output/metaX/neg/sampleList.txt",
  output_dir="Output/4.MetaboliteComparison",
  pvalue=0.05,
  log2ratio=1,
  vip=1,
  pvalue_type="t.test_p.value",
  mode="all"
);

p <- arg_parser("section4 metabolite comparison")
p <- add_argument(p, "--quant_neg", help="metaX neg quant data")
p <- add_argument(p, "--quant_pos", help="metaX pos quant data")
p <- add_argument(p, "--intensity", help="intensity file")
p <- add_argument(p, "--sample_file", help="sample file")
p <- add_argument(p, "--annotation", help="annotation file")
p <- add_argument(p, "--output_dir", help="output dir")
p <- add_argument(p, "--pvalue", help="pvalue cutoff", default=0.05)
p <- add_argument(p, "--log2ratio", help="log2ratio cutoff", default=1.5)
p <- add_argument(p, "--vip", help="vip cutoff", default=1)
p <- add_argument(p, "--compound", help="compound file")
p <- add_argument(p, "--pvalue_type", help="pvalue_type, default=t.test_p.value", default='t.test_p.value')
p <- add_argument(p, "--mode", help="mode, default=all", default='all')
p <- add_argument(p, "--kgml", help="kgml path")

argv <- parse_args(p)

# data prepare
outdir = argv$output_dir
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
# metaX quant data
quant_data_neg = read.csv(argv$quant_neg, sep="\t")
quant_data_neg$ID = paste0("neg-", quant_data_neg$ID)
quant_data_pos = read.csv(argv$quant_pos, sep="\t")
quant_data_pos$ID = paste0("pos-", quant_data_pos$ID)

quant_data <- rbind(quant_data_neg, quant_data_pos)
groups <- unique(quant_data$sample)
# sample list data
sample_list <- read.csv(file=argv$sample_file, sep="\t", stringsAsFactors=FALSE)
sample_list[is.na(sample_list$class), 'class'] <- "QC"
rownames(sample_list) <- sample_list$sample
# annotation data
annotation <- read.csv(file=argv$annotation)
rownames(annotation) <- annotation$ID
# intensity data
intensity <- as.list(read.csv(argv$intensity, nrows=1, stringsAsFactors=FALSE, colClasses = c("character")))

# kegg analysis script

get_current_script_path <- function() {
        cmdArgs <- commandArgs(trailingOnly = FALSE)
        needle <- "--file="
        match <- grep(needle, cmdArgs)
        if (length(match) > 0) {
                # Rscript
                return(normalizePath(sub(needle, "", cmdArgs[match])))
        } else {
                # 'source'd via R console
                return(normalizePath(sys.frames()[[1]]$ofile))
        }
}

kegg_analysis_script <- paste(dirname(dirname(get_current_script_path())), "kegg_analysis.py", sep="/")
kegg_plot_script <- paste(dirname(dirname(get_current_script_path())), "2.MetaboliteIdentification/KEGG.R", sep="/")
kegg_color_script <- paste(dirname(dirname(get_current_script_path())), "color_pathway.R", sep="/")

# init stat
stats <- data.frame(comparison=character(0), pos_all=numeric(0), pos_up=numeric(0), pos_down=numeric(0), neg_all=numeric(0), neg_up=numeric(0), neg_down=numeric(0), all_regulate=numeric(0), stringsAsFactors = F)
for(group in groups){
  print(group)
  group_data = quant_data[quant_data$sample==group, ]
  group_new_name <- str_replace(group, ":", "_")
  group1 = unlist(str_split(group, ":"))[1]
  group2 = unlist(str_split(group, ":"))[2]
  group_out_dir = file.path(outdir, group_new_name)
  dir.create(group_out_dir, recursive =TRUE, showWarnings=FALSE)
  for(ion_type in c("neg", "pos")){
      if(ion_type == "neg"){
        pca_object <- readRDS(paste0(dirname(argv$quant_neg), "/", ion_type, "-", group_new_name, "-pca.rds"))
        plsda_model <- readRDS(paste0(dirname(argv$quant_neg), "/", ion_type, "-", group_new_name, "-plsDAmodel.rds"))
        file.copy(paste0(dirname(argv$quant_neg), "/", ion_type, "-", group_new_name, "-PLSDA-score.png"), paste0(group_out_dir, "/", group_new_name, ".", ion_type, ".plsda.score.png"))
      }else{
        pca_object <- readRDS(paste0(dirname(argv$quant_pos), "/", ion_type, "-", group_new_name, "-pca.rds"))
        plsda_model <- readRDS(paste0(dirname(argv$quant_pos), "/", ion_type, "-", group_new_name, "-plsDAmodel.rds"))
        file.copy(paste0(dirname(argv$quant_pos), "/", ion_type, "-", group_new_name, "-PLSDA-score.png"), paste0(group_out_dir, "/", group_new_name, ".", ion_type, ".plsda.score.png"))
      }

      pca_score = as.data.frame(pca_object@scores)
      pca_score <- merge(pca_score, sample_list, by=0)

      pca_loading = as.data.frame(pca_object@loadings)

      pca_data1 = data.frame(rownames(pca_loading),pca_loading$PC1,pca_loading$PC2)
      colnames(pca_data1) = c("ID","PC1","PC2")
      pca_data1$ID = paste0(ion_type, "-", pca_data1$ID)
      pca_data2 = data.frame(group_data$ID,group_data$VIP)
      colnames(pca_data2) = c("ID","VIP")

      pca_data_merge = merge(pca_data1,pca_data2,by=c("ID"))
      # merge pca loading with annotation
      rownames(pca_loading) = paste0(ion_type, "-", rownames(pca_loading))
      pca_loading_total = merge(pca_loading, annotation, by=0)

      ## 区分 neg / pos
      draw_pca_score_plot(pca_score, group_out_dir, group_new_name, ion_type)
      draw_pca_loading_plot(pca_loading_total, group_out_dir, group_new_name, ion_type)
      draw_plsad_loading(pca_data_merge, group_out_dir, group_new_name, ion_type)
      draw_plsda_validation_plot(plsda_model, paste0(group_out_dir, "/", group_new_name, ".", ion_type, ".plsda.validation"))
  }

  #merge annotation with quant
  target_sample_cols = names(which(intensity == group1 | intensity == group2))
  group_annotation = annotation[, c("ID", target_sample_cols, "MZ", "RT", "Metabolite", "Superclass", "Class", "HMDB", "KEGG", 'Level1', 'Level2', 'Level3', 'Pathway')]
  rownames(group_annotation) = annotation$ID

  rownames(group_data) <- group_data$ID
  group_data_total <- merge(group_annotation, group_data, by=0)

  ## 不区分 neg / pos
  draw_sig_metabolite_heatmap(group_data_total, group_out_dir, group_new_name, pca_score$class, target_sample_cols, argv$pvalue, argv$log2ratio, argv$pvalue_type, argv$vip)
  draw_sig_heatmap(group_data_total, group_out_dir, group_new_name, pca_score$class, target_sample_cols, argv$pvalue, argv$log2ratio, argv$pvalue_type, argv$vip)

  draw_volcano_sig_meta(group_data_total, group_out_dir, group_new_name, argv$pvalue, argv$log2ratio, argv$pvalue_type, argv$vip)
  draw_volcano(group_data_total, group_out_dir, group_new_name, argv$pvalue, argv$log2ratio, argv$pvalue_type, argv$vip)

  ## export all
  group_data_total$regulated <- "none"

  ## significant p < 0.05 and ratio > 2 or ratio < 2 and VIP > 1
  pvaue_cutoff_log10 <- -log10(argv$pvalue)
  FC = log2(group_data_total$ratio)
  log_Q_value = -log10(group_data_total[argv$pvalue_type])

  group_data_total[which(FC> argv$log2ratio & log_Q_value>pvaue_cutoff_log10 & group_data_total$VIP > argv$vip), 'regulated']="up"
  group_data_total[which(FC< -argv$log2ratio & log_Q_value>pvaue_cutoff_log10 & group_data_total$VIP > argv$vip), 'regulated']="down"

  # drop some columns
  group_data_total<-within(group_data_total, rm("Row.names", "ID.y"))
  colnames(group_data_total) <- c("ID", colnames(group_data_total[2:ncol(group_data_total)]))
  write.csv(group_data_total, file=paste0(group_out_dir, "/", group_new_name, ".all.csv"), row.names=F, sep=",")

  group_data_total.significant <- group_data_total[which(group_data_total$regulated != "none"), ]

  write.csv(group_data_total.significant, file=paste0(group_out_dir, "/", group_new_name, ".significant.csv"), row.names=F, sep=",")

  ## exec kegg analysis script
  print(paste("python ", kegg_analysis_script, " --input ", paste0(group_out_dir, "/", group_new_name, ".significant.csv") ," --outdir ", paste0(group_out_dir, "/KEGG"), " --compound" , argv$compound, "--mode", argv$mode, sep=" "))
  if(dim(group_data_total.significant)[1] > 0){
    system(paste("python", kegg_analysis_script, " --input ", paste0(group_out_dir, "/", group_new_name, ".significant.csv") ," --outdir ", paste0(group_out_dir, "/KEGG"), " --compound" , argv$compound, "--mode", argv$mode, sep=" "))
    system(paste(rscript_bin,  kegg_plot_script, "-l" , paste0(group_out_dir, "/KEGG/kegg.level2.csv"), "-k", paste0(group_out_dir, "/KEGG/kegg.csv"), "-o", paste0(group_out_dir, "/KEGG"), sep=" "))
    system(paste(rscript_bin,  kegg_color_script, "--kegg_result" , paste0(group_out_dir, "/KEGG/kegg.csv"), "--kegg_dir", argv$kgml, "--de_result", paste0(group_out_dir, "/", group_new_name, ".significant.csv"), "--out_dir", paste0(group_out_dir, "/KEGG/kegg_map"), sep=" "))

  }
  stats[nrow(stats) + 1,] = c(paste(group1, group2, sep="/"),
                              sum(grepl("pos-", group_data_total.significant$ID)),
                              sum(grepl("pos-", group_data_total.significant$ID) & (group_data_total.significant$regulated == "up")),
                              sum(grepl("pos-", group_data_total.significant$ID) & (group_data_total.significant$regulated == "down")),
                              sum(grepl("neg-", group_data_total.significant$ID)),
                              sum(grepl("neg-", group_data_total.significant$ID) & (group_data_total.significant$regulated == "up")),
                              sum(grepl("neg-", group_data_total.significant$ID) & (group_data_total.significant$regulated == "down")),
                              nrow(group_data_total.significant))
}
# summary stats plot
draw_stat_data(stats, outdir)
write.csv(stats, file=paste0(outdir, "/comparison.stat.csv"), row.names=F, sep=",")
