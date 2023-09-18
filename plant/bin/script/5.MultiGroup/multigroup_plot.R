suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Mfuzz))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(htmltools))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(argparser))
suppressPackageStartupMessages(library(ropls, quietly = T))

draw_stat_plot <- function(multigroup.significant.stat, outdir) {
  stat_matrix = c(
    multigroup.significant.stat$significant
    # multigroup.significant.stat$significant_MS2
  )
  # significant_group = colnames(multigroup.significant.stat)[c(3, 4)]
  significant_group = colnames(multigroup.significant.stat)[c(3)]
  stat_df = data.frame(significant_group, stat_matrix)
  colnames(stat_df) = c("name", "Number_of_Feature")
  
  
  p <-ggplot(data = stat_df) + 
    geom_bar(aes(x = significant_group,
                 y = Number_of_Feature, 
                 fill =Number_of_Feature),
             stat = "identity") +
    theme_bw() +
    geom_text(aes(x = significant_group, y = Number_of_Feature, label =
                    Number_of_Feature),
              vjust = 0,size=2.5) +
    ylab("Number of features") +
    theme(panel.grid = element_blank(),
          legend.position = "")
  ggsave(paste0(outdir, "/multigroup.significant.stat.png"), p, width = 1200,height=800,units = "px")
}

draw_anova_plsda_plot <- function(anova_plsda_file,outdir, p_cut=0.05, v_cut=1, pvalue_type="pvalue"){
  log_pvalue = -log10(anova_plsda_file[, pvalue_type])
  vip = anova_plsda_file$VIP
  regulated = rep("none",length(vip))
  regulated[which(vip>v_cut & anova_plsda_file[, pvalue_type]<p_cut)]="significant"
  anova_experss_df = data.frame(log_pvalue,vip,regulated)
  
  p <-ggplot(data=anova_experss_df, aes(x=vip, y=log_pvalue,color=regulated,shape=regulated)) + geom_point(alpha=1, size=2) +
    geom_hline(yintercept =-log10(p_cut), linetype = 3,size =1) +
    geom_vline(xintercept =  v_cut, linetype = 3,size =1) +
    annotate("text",x = 1.5*v_cut,y=-log10(p_cut),label = paste0("pvalue<",p_cut) )+
    annotate("text",x = v_cut,y=5*-log10(p_cut),label = paste0("VIP>",v_cut) )+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      legend.title = element_text(size =10),
      legend.text = element_text(size =8)
    ) + 
    xlab("VIP")+
    ylab("-log10(Anova_pvalue)")+
    scale_color_manual(values = c("#BEBEBE","#EE0000"),breaks = c("none","significant"))+
    scale_shape_manual(values = c(1,16),breaks =c("none","significant"))
  ggsave(paste0(outdir,"/anova_plsda.png"),p,width = 1800,height = 1200,units = "px")
  ggsave(paste0(outdir,"/anova_plsda.pdf"),p,width = 1800,height = 1200,units = "px")
}

draw_plsda_score_plot = function(plsda_score_file,outdir, pca){
  plsda_score_noqc = plsda_score_file[which(plsda_score_file["class"]!="QC"),]

  colourCount = length(unique(plsda_score_file$class))
  colourCount_noqc = length(unique(plsda_score_file$class))
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))

  p_2d = ggplot(data = plsda_score_file,aes(x=plsda_score_file[,2],y=plsda_score_file[,3],colour=class))+
    geom_point()+
    theme_bw()+
    scale_color_manual(values=getPalette(colourCount))+
    stat_ellipse(level = 0.9)+
    xlab(paste0("PC1", "(", pca[1], "%)"))+
    ylab(paste0("PC2", "(", pca[2], "%)"))+
    theme(legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold",size=14),
          axis.title = element_text(face = "bold"),
          panel.grid=element_blank(),
    )

  p_2d_no_qc = ggplot(data = plsda_score_noqc,aes(x=plsda_score_noqc[,2],y=plsda_score_noqc[,3],colour=class))+
    geom_point()+
    theme_bw()+
    scale_color_manual(values=getPalette(colourCount_noqc))+
    stat_ellipse(level = 0.9)+
    xlab(paste0("PC1", "(", pca[1], "%)"))+
    ylab(paste0("PC2", "(", pca[2], "%)"))+
    theme(legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold",size=14),
          axis.title = element_text(face = "bold"),
          panel.grid=element_blank(),
    )

  ggsave(filename = paste0(outdir,"/plsda.score.png"),p_2d,width = 1500,height = 1200,units = "px")
  ggsave(filename = paste0(outdir,"/plsda.score.pdf"),p_2d,width = 1500,height = 1200,units = "px")
  ggsave(filename = paste0(outdir,"/plsda.score_no_qc.png"),p_2d_no_qc,width = 1500,height = 1200,units = "px")
  ggsave(filename = paste0(outdir,"/plsda.score_no_qc.pdf"),p_2d_no_qc,width = 1500,height = 1200,units = "px")
  p_html = ggplotly(p_2d,width = 1500,height = 1200)
  save_html(p_html,paste0(outdir,"/plsda.score.2d.html"))


  # 3D plot
  p_html_3d <- plot_ly(data = plsda_score_file, x = plsda_score_file[,2], y =  plsda_score_file[,3], z = plsda_score_file[,4],
                       color = ~class, colors = getPalette(colourCount))
  p_html_3d <- p_html_3d %>% add_markers()
  p_html_3d <- p_html_3d %>% layout(scene = list(xaxis = list(title = 'P1'),
                                                 yaxis = list(title = 'P2'),
                                                 zaxis = list(title = 'P3')))
  save_html(p_html_3d,paste0(outdir,"/plsda.score.3d.html"))
}

draw_sig_eachsample_metabolite_heatmap = function(significant.eachsample.class.mfuzz_file,target_sample_cols, target_sample_class, outdir){
  all.feature_file.metabolite = significant.eachsample.class.mfuzz_file[which(significant.eachsample.class.mfuzz_file$Metabolite!=""),]

  meta_intensity = all.feature_file.metabolite[, target_sample_cols]
  cluster= all.feature_file.metabolite$cluster
  superClass = all.feature_file.metabolite$Superclass

  # make annotation col
  annotation_df_col = as.data.frame(target_sample_class)
  row.names(annotation_df_col) = colnames(meta_intensity)

  meta_intensity_numeric = apply(meta_intensity, 2, as.numeric)
  intensity_df = as.data.frame(meta_intensity_numeric)
  row.names(intensity_df) = as.character(paste0(all.feature_file.metabolite$Metabolite,"@",all.feature_file.metabolite$ID))

  # make annotation row
  annotation_df_row = data.frame(cluster=as.factor(cluster),Metaboliteclass=superClass)
  row.names(annotation_df_row) = row.names(intensity_df)

  pheatmap(
    intensity_df,
    border_color = "#666666",
    show_rownames = T,
    show_colnames = F,
    fontsize_row = 3,
    annotation_col = annotation_df_col,
    annotation_row = annotation_df_row,
    annotation_names_col = F,
    annotation_names_row = F,
    cluster_cols = F,
    cluster_rows = F,
    scale = "row",
    colorRampPalette(c("blue", "white", "red"))(50),
    na_col=0,
    filename = paste0(outdir,"/significant.eachsample.mfuzz.Metabolite.heatmap.png"),
    width = 10,
    height = 20
  )

  pheatmap(
    intensity_df,
    border_color = "#666666",
    show_rownames = T,
    show_colnames = F,
    fontsize_row = 3,
    annotation_col = annotation_df_col,
    annotation_row = annotation_df_row,
    annotation_names_col = F,
    annotation_names_row = F,
    cluster_cols = F,
    cluster_rows = F,
    scale = "row",
    colorRampPalette(c("blue", "white", "red"))(50),
    na_col=0,
    filename = paste0(outdir,"/significant.eachsample.mfuzz.Metabolite.heatmap.pdf"),
    width = 10,
    height = 20
  )
}

draw_sig_eachsample_heatmap = function(significant.eachsample.class.mfuzz_file,target_sample_cols, target_sample_class, outdir){
  all.feature_file.metabolite = significant.eachsample.class.mfuzz_file

  meta_intensity = all.feature_file.metabolite[, target_sample_cols]
  cluster= all.feature_file.metabolite$cluster
  superClass = all.feature_file.metabolite$Superclass

  # make annotation col
  annotation_df_col = as.data.frame(target_sample_class)
  row.names(annotation_df_col) = colnames(meta_intensity)

  meta_intensity_numeric = apply(meta_intensity, 2, as.numeric)
  intensity_df = as.data.frame(meta_intensity_numeric)
  row.names(intensity_df) = as.character(paste0(all.feature_file.metabolite$Metabolite, "@",all.feature_file.metabolite$ID))

  # make annotation row
  annotation_df_row = data.frame(cluster=as.factor(cluster))
  row.names(annotation_df_row) = row.names(intensity_df)

  pheatmap(
    intensity_df,
    border_color = "#666666",
    show_rownames = F,
    show_colnames = T,
    fontsize_row = 3,
    annotation_col = annotation_df_col,
    annotation_row = annotation_df_row,
    annotation_names_col = F,
    annotation_names_row = F,
    cluster_cols = T,
    cluster_rows = F,
    scale = "row",
    colorRampPalette(c("blue", "white", "red"))(50),
    na_col=0,
    filename = paste0(outdir,"/significant.eachsample.mfuzz.heatmap.png"),
    width = 12,
    height = 8
  )

  pheatmap(
    intensity_df,
    border_color = "#666666",
    show_rownames = T,
    show_colnames = F,
    fontsize_row = 3,
    annotation_col = annotation_df_col,
    annotation_row = annotation_df_row,
    annotation_names_col = F,
    annotation_names_row = F,
    cluster_cols = F,
    cluster_rows = F,
    scale = "row",
    colorRampPalette(c("blue", "white", "red"))(50),
    na_col=0,
    filename = paste0(outdir,"/significant.eachsample.mfuzz.heatmap.pdf"),
    width = 12,
    height = 8
  )
}

draw_sig_class_metabolite_heatmap = function(significant.class.mfuzz_file, sample_class, outdir){
  all.feature_file.metabolite = significant.class.mfuzz_file[which(significant.class.mfuzz_file$Metabolite!=""),]

  meta_intensity = all.feature_file.metabolite[, sample_class]
  cluster= all.feature_file.metabolite$cluster

  intensity_df = meta_intensity
  row.names(intensity_df) = as.character(paste0(all.feature_file.metabolite$Metabolite,"@",all.feature_file.metabolite$ID))

  # make annotation row
  annotation_df_row = data.frame(cluster=as.character(cluster))
  row.names(annotation_df_row) = row.names(intensity_df)

  pheatmap(
    intensity_df,
    cellwidth = 20,
    border_color = "#666666",
    show_rownames = T,
    show_colnames = F,
    fontsize_row = 2,
    annotation_row = annotation_df_row,
    annotation_names_col = F,
    annotation_names_row = F,
    cluster_cols = F,
    cluster_rows = F,
    scale = "row",
    colorRampPalette(c("blue", "white", "red"))(50),
    na_col=0,
    filename = paste0(outdir,"/significant.mfuzz.Metabolite.heatmap.png"),
    width = 8,
    height = 12
  )

  pheatmap(
    intensity_df,
    cellwidth = 20,
    border_color = "#666666",
    show_rownames = T,
    show_colnames = F,
    fontsize_row = 2,
    annotation_row = annotation_df_row,
    annotation_names_col = F,
    annotation_names_row = F,
    cluster_cols = F,
    cluster_rows = F,
    scale = "row",
    colorRampPalette(c("blue", "white", "red"))(50),
    na_col=0,
    filename = paste0(outdir,"/significant.mfuzz.Metabolite.heatmap.pdf"),
    width = 8,
    height = 12
  )
}

draw_sig_class_heatmap = function(significant.class.mfuzz_file, sample_class, outdir){

  meta_intensity = significant.class.mfuzz_file[, sample_class]
  cluster= significant.class.mfuzz_file$cluster

  intensity_df = meta_intensity

  # make annotation row
  annotation_df_row = data.frame(cluster=as.character(cluster))
  row.names(annotation_df_row) = row.names(intensity_df)

  pheatmap(
    intensity_df,
    border_color = "#666666",
    show_rownames = F,
    show_colnames = T,
    fontsize_row = 5,
    annotation_row = annotation_df_row,
    annotation_names_col = F,
    annotation_names_row = F,
    cluster_cols = T,
    cluster_rows = F,
    scale = "row",
    colorRampPalette(c("blue", "white", "red"))(50),
    na_col=0,
    filename = paste0(outdir,"/significant.class.mfuzz.heatmap.png"),
    width = 8,
    height = 8
  )

    pheatmap(
    intensity_df,
    border_color = "#666666",
    show_rownames = F,
    show_colnames = T,
    fontsize_row = 5,
    annotation_row = annotation_df_row,
    annotation_names_col = F,
    annotation_names_row = F,
    cluster_cols = T,
    cluster_rows = F,
    scale = "row",
    colorRampPalette(c("blue", "white", "red"))(50),
    na_col=0,
    filename = paste0(outdir,"/significant.class.mfuzz.heatmap.pdf"),
    width = 8,
    height = 8
  )
}

draw_sig_mfuzz_class = function(data_sig_class_mfuzz, zscore, sample_class, outdir){
  Sample = rep(sample_class,each=nrow(data_sig_class_mfuzz))
  Cluster = rep(data_sig_class_mfuzz$cluster,length(sample_class))
  Membership = rep(data_sig_class_mfuzz$membership,length(sample_class))
  Z_score = as.vector(unlist(zscore))
  Group = rep(data_sig_class_mfuzz$ID,length(sample_class))


  n_row = round(length(unique(Cluster))^0.5)

  data_group_mfuzz_df = data.frame(Sample=Sample,Cluster=Cluster,Membership=Membership,Z_score=Z_score,Group =Group)

  p = ggplot(data = data_group_mfuzz_df,aes(x=Sample,y=Z_score,color=Membership,group=Group)) +
    geom_line(alpha=0.5)+
    facet_wrap(~Cluster,nrow = n_row)+
    theme_bw()+
    ylab("Z-score")+
    scale_colour_gradientn(colours = terrain.colors(6))+
    theme(panel.grid = element_blank())
  ggsave(paste0(outdir,"/significant.class.mfuzz.png"),p,width = n_row*900,height = n_row*900,units = "px")
  ggsave(paste0(outdir,"/significant.class.mfuzz.pdf"),p,width = n_row*900,height = n_row*900,units = "px")
}

draw_sig_mfuzz_eachsample = function(significant.eachsample.class.mfuzz_file, zscore, target_sample_cols, target_sample_class, outdir){

  Sample = rep(target_sample_cols,each=nrow(significant.eachsample.class.mfuzz_file))
  Cluster = rep(significant.eachsample.class.mfuzz_file$cluster,length(target_sample_cols))
  Membership = rep(significant.eachsample.class.mfuzz_file$membership,length(target_sample_cols))
  Z_score = as.vector(unlist(zscore))
  Group = rep(significant.eachsample.class.mfuzz_file$ID, length(target_sample_cols))

  n_row = round(length(unique(Cluster))^0.5)

  data_group_mfuzz_df = data.frame(Sample=Sample,Cluster=Cluster,Membership=Membership,Z_score=Z_score,Group =Group)

  p = ggplot(data = data_group_mfuzz_df,aes(x=Sample,y=Z_score,color=Membership,group=Group)) +
    geom_line(alpha=0.5)+
    facet_wrap(~Cluster,nrow = n_row)+
    theme_bw()+
    ylab("Z-score")+
    scale_colour_gradientn(colours = terrain.colors(6))+
    theme(panel.grid = element_blank(),
            axis.text.x = element_text(angle=85,hjust=1))
  ggsave(paste0(outdir,"/significant.class.eachsample.mfuzz.png"),p,width = n_row*900,height = n_row*900,units = "px")
  ggsave(paste0(outdir,"/significant.class.eachsample.mfuzz.pdf"),p,width = n_row*900,height = n_row*900,units = "px")
}

draw_plsda_validation_plot = function(plsda_object, out_file){
  plotdat <- data.frame(results.plsda@suppLs$permMN) %>%
    #rename relevant columns
    rename(R2Y = R2Y.cum., Q2Y = Q2.cum.) %>%

    #Select columns relevant to plot and turn into a tidy/long df
    select(R2Y, Q2Y, sim)

  colnames(plotdat) <- c("R2", "Q2", "cor")
  x1 <- plotdat$cor[1]
  y1 <- plotdat$R2[1]
  y2 <- plotdat$Q2[1]
  plotdat$cor <- abs(plotdat$cor)
  plotdat <- plotdat[order(plotdat$cor),]
  png(file = paste0(out_file, ".png"))
  par(mar=c(3,3,2,1),mgp=c(1.6,0.6,0),cex.lab=1.2,cex.main=0.9)
  plot(plotdat$cor,plotdat$R2,ylim=c(min(plotdat$R2,plotdat$Q2),1),xlim=c(0, 2),pch=16,
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
  dev.off()

  pdf(file = paste0(out_file, ".pdf"))
  par(mar=c(3,3,2,1),mgp=c(1.6,0.6,0),cex.lab=1.2,cex.main=0.9)
  plot(plotdat$cor,plotdat$R2,ylim=c(min(plotdat$R2,plotdat$Q2),1),,xlim=c(0, 2), pch=16,
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
  dev.off()
}

#__________________main____________________

argv <- list(
  annotation="Output/2.MetaboliteIdentification/identification.raw.intensity.csv",
  output_dir="Output/5.MultiGroup",
  intensity="Output/metaX/neg/data/neg-norm-metaboAnalystInput.csv",
  sample_file="Output/metaX/neg/sampleList.txt",
  pvalue=0.05,
  vip=1,
  mode="all"
);

p <- arg_parser("multi-group analysis, anovar analysis/VIP")
p <- add_argument(p, "--annotation", help="input all featue.xls")
p <- add_argument(p, "--intensity", help="metaX intensity output")
p <- add_argument(p, "--sample_file", help="sample file")
p <- add_argument(p, "--compound", help="compound file")
p <- add_argument(p, "--output_dir", help="output dir")
p <- add_argument(p, "--pvalue", help="pvalue cutoff", default=0.05)
p <- add_argument(p, "--vip", help="vip cutoff", default=1)
p <- add_argument(p, "--pvalue_type", help="pvalue, default is pvalue, options are pvalue, fdr", default="pvalue")
p <- add_argument(p, "--mode", help="mode", default="all")
p <- add_argument(p, "--kgml", help="kgml path")

argv <- parse_args(p)

# print(argv)
outdir = argv$output_dir
pvalue_cutoff <- argv$pvalue
vip_cutoff <- argv$vip

sample_list <- read.csv(file=argv$sample_file, sep="\t", stringsAsFactors=FALSE)
sample_list[is.na(sample_list$class), 'class'] <- "QC"
rownames(sample_list) <- sample_list$sample

annotation <- read.csv(file=argv$annotation)
rownames(annotation) <- annotation$ID
intensity <- as.list(read.csv(argv$intensity, nrows=1, row.names=1, stringsAsFactors=FALSE, colClasses = c("character")))
target_sample_cols <- names(which(intensity != "QC"))
target_sample_class <- as.character(intensity[which(intensity != "QC")])
unique_classes = unique(target_sample_class)
group_out_dir = paste(outdir, paste(unique_classes, collapse="~"), sep="/")
dir.create(group_out_dir, showWarnings=FALSE, recursive=TRUE)

## calculate annovar p value
annotation[is.na(annotation)] <- 0
print("calculate annovar pvalue")
for(i in 1:nrow(annotation)){
  tmp <- data.frame(
    treatment = as.factor(target_sample_class),
    yield = unlist(annotation[i, target_sample_cols], use.names=FALSE)
  )
  one.way <- aov(yield ~ treatment, data = tmp)
  annotation[i, "pvalue"] <- summary(one.way)[[1]][["Pr(>F)"]][1]
}
annotation[, "fdr"] <- p.adjust(annotation[, 'pvalue'], method="fdr")

## calculate VIP
annotation.t <- t(annotation[, target_sample_cols])
# get from params
# results.plsda <- plsda(annotation.t, target_sample_class, ncomp=3,max.iter=100)
# VIP <- mixOmics::vip(results.plsda )
# VIP.value <- as.data.frame(VIP[order(VIP[,ncol(VIP)],decreasing=TRUE),ncol(VIP)])
# colnames(VIP.value) <- "VIP"
# annotation <- merge(annotation, VIP.value, by=0)
# use ropls do plsda

print("calculate VIP")
results.plsda <- opls(annotation.t, target_sample_class, fig.pdfC='none', info.txtC='none', predI = 3)
draw_plsda_validation_plot(results.plsda, paste0(group_out_dir, "/", paste(unique_classes, collapse="~"), "-plsda.permutation"))
vip <- data.frame(VIP=getVipVn(results.plsda))
annotation <- merge(x=annotation, y=vip, by=0)
annotation <- within(annotation, rm("Row.names"))
rownames(annotation) <- annotation$ID

annotation$regulated <- 'none'
annotation[annotation$VIP > vip_cutoff & annotation[, argv$pvalue_type] < pvalue_cutoff , "regulated"] = "significant"

results.plsda.summary <- results.plsda@summaryDF
results.plsda.summary$Group <- paste(unique_classes, collapse=",")
write.csv(results.plsda.summary, paste0(outdir, "/multigroup.plsda.stat.csv"), row.names=F,)


## significant VIP > 1 and p.value < 0.05
significant <- annotation[(annotation$VIP > vip_cutoff) & (annotation[, argv$pvalue_type] < pvalue_cutoff), ]
rownames(significant) <- significant$ID
## mfuzz analysis by sample
est_sample.data <- as.matrix(significant[, target_sample_cols])
est_sample <- new("ExpressionSet",exprs = est_sample.data)
est_sample.std <- standardise(est_sample)
set.seed(48)
c_sample <- 16
m_sample <- mestimate(est_sample.std)
cl_sample <- mfuzz(est_sample.std, c = c_sample, m = m_sample)
cl_sample_membership <- data.frame(cluster=cl_sample$cluster, membership=apply(cl_sample$membership, 1, max))

## mfuzz analysis by class

est_class.data <- matrix(, nrow=nrow(significant),ncol=length(unique_classes))
for(i in 1:length(unique_classes)){
  est_class.data[, i] <- as.numeric(apply(significant[, names(which(intensity == unique_classes[i]))], 1, median ,na.rm=T))
}
rownames(est_class.data) = rownames(significant)
est_class <- new("ExpressionSet",exprs = est_class.data)
est_class.std <- standardise(est_class)
set.seed(48)
c_class <- 9
m_class <- mestimate(est_class.std)
cl_class <- mfuzz(est_class.std, c = c_class, m = m_class)
cl_class_membership <- data.frame(cluster=cl_class$cluster, membership=apply(cl_class$membership, 1, max))

# plsda score data
plsda_score_data <- as.data.frame(results.plsda@scoreMN)
plsda_score_data$class <- target_sample_class
pca <- results.plsda@modelDF$R2Y * 100

## plot
draw_anova_plsda_plot(annotation, group_out_dir, pvalue_cutoff, vip_cutoff, argv$pvalue_type)
draw_plsda_score_plot(plsda_score_data, group_out_dir, pca)

significant.class.mfuzz <- as.data.frame(est_class.data)
colnames(significant.class.mfuzz) <- unique_classes
significant.class.mfuzz$ID <- rownames(significant.class.mfuzz)
significant.class.mfuzz <- merge(significant.class.mfuzz, cl_class_membership, by=0)
significant.class.mfuzz <- within(significant.class.mfuzz, rm("Row.names"))
rownames(significant.class.mfuzz) <- significant.class.mfuzz$ID
significant.class.mfuzz <-merge(significant.class.mfuzz, annotation[, c("Metabolite", "Superclass")], by=0)
significant.class.mfuzz <- within(significant.class.mfuzz, rm("Row.names"))
significant.class.mfuzz <- significant.class.mfuzz[order(significant.class.mfuzz$cluster), ]
rownames(significant.class.mfuzz) <- significant.class.mfuzz$ID

draw_sig_class_heatmap(significant.class.mfuzz, unique_classes, group_out_dir)
draw_sig_class_metabolite_heatmap(significant.class.mfuzz,unique_classes, group_out_dir)


significant.sample.mfuzz <- merge(significant, cl_sample_membership, by=0)
significant.sample.mfuzz <- within(significant.sample.mfuzz, rm("Row.names"))
rownames(significant.sample.mfuzz) <- significant.sample.mfuzz$ID


draw_sig_eachsample_heatmap(significant.sample.mfuzz, target_sample_cols, target_sample_class, group_out_dir)
draw_sig_eachsample_metabolite_heatmap(significant.sample.mfuzz, target_sample_cols, target_sample_class, group_out_dir)

mfuzz_zscore_sample <- est_sample.std@assayData$exprs
mfuzz_zscore_sample <- mfuzz_zscore_sample[rownames(significant.sample.mfuzz), ]

draw_sig_mfuzz_eachsample(significant.sample.mfuzz, mfuzz_zscore_sample, target_sample_cols, target_sample_class, group_out_dir)

mfuzz_zscore_class <- est_class.std@assayData$exprs
mfuzz_zscore_class <- mfuzz_zscore_class[rownames(significant.class.mfuzz), ]
draw_sig_mfuzz_class(significant.class.mfuzz, mfuzz_zscore_class, unique_classes, group_out_dir)

#stats <- data.frame(multigroup=paste(unique_classes, collapse="-"), all=nrow(annotation), significant=nrow(significant), significant_MS2=length(which(significant$MS2!="No MS2")))
stats <- data.frame(multigroup=paste(unique_classes, collapse="-"), all=nrow(annotation), significant=nrow(significant))
write.csv(stats, paste(outdir, "multigroup.significant.stat.csv", sep="/"), row.names=F)
draw_stat_plot(stats, outdir)


write.csv(annotation, paste(group_out_dir, "all.feature.csv", sep="/"), row.names=F)
write.csv(significant, paste(group_out_dir, "significant.feature.csv", sep="/"), row.names=F)

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

system(paste("python", kegg_analysis_script, " --input ", paste0(group_out_dir, "/significant.feature.csv") ," --outdir ", paste0(group_out_dir, "/KEGG"), " --compound" , argv$compound, "--mode", argv$mode, sep=" "))
system(paste("Rscript", kegg_plot_script, "-l" , paste0(group_out_dir, "/KEGG/kegg.level2.csv"), "-k", paste0(group_out_dir, "/KEGG/kegg.csv"), "-o", paste0(group_out_dir, "/KEGG"), sep=" "))
system(paste("Rscript", kegg_color_script, "--kegg_result" , paste0(group_out_dir, "/KEGG/kegg.csv"), "--kegg_dir", argv$kgml, "--out_dir", paste0(group_out_dir, "/KEGG/kegg_map"), sep=" "))

