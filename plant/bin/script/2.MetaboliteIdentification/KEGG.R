library(ggplot2)
library(argparser)
draw_kegg_stat_bar_plot= function(data_level2,data_level3,outdir){
  data_level2 <- data_level2[order(data_level2$Level1), ]
  data_level2$Level2  = factor(data_level2$Level2,levels = data_level2$Level2)

  deta_width_level2 = length(data_level2$Level2)*70
  if (deta_width_level2<4800) {
    deta_width_level2 = 4800
  }

  deta_width_Pathway = length(data_level3$Pathway)*70
  if (deta_width_Pathway<4800) {
    deta_width_Pathway = 4800
  }
  p=ggplot(data = data_level2,aes(x = Level2,y= NumberFeature,fill=Level1)) +
    geom_bar(stat ="identity",na.rm = TRUE)+
    theme_bw()+
    ylab("Number of Feature")+
    xlab("")+
    theme(
      axis.text.x = element_text(angle = 80,hjust = 1,size=10), # Remove x axis tick labels
      axis.title=element_text(size=14,face="bold"),
      legend.title = element_text(face = "bold"),
      panel.grid=element_blank()
    )
  ggsave(filename = paste0(outdir,"/kegg.level2.feature.png"),p,height = 3200,width = deta_width_level2,units ="px", limitsize = FALSE)
  ggsave(filename = paste0(outdir,"/kegg.level2.feature.pdf"),p,height = 3200,width = deta_width_level2,units ="px", limitsize = FALSE)

  p1 = p=ggplot(data = data_level2,aes(x = Level2,y= NumberCompound,fill=Level1)) +
    geom_bar(stat ="identity",na.rm = TRUE)+
    theme_bw()+
    ylab("Number of Compound")+
    xlab("")+
    theme(
      axis.text.x = element_text(angle = 80,hjust = 1,size=10), # Remove x axis tick labels
      axis.title=element_text(size=14,face="bold"),
      legend.title = element_text(face = "bold"),
      panel.grid=element_blank()
    )
  ggsave(filename = paste0(outdir,"/kegg.level2.png"),p,height = 3200,width = deta_width_level2,units ="px", limitsize = FALSE)
  ggsave(filename = paste0(outdir,"/kegg.level2.pdf"),p,height = 3200,width = deta_width_level2,units ="px", limitsize = FALSE)

  data_level3$Pathway  = factor(data_level3$Pathway,levels = data_level3$Pathway)
  data_level3 <- data_level3[order(data_level3$Level2), ]
  p2=ggplot(data = data_level3,aes(x =Pathway,y= NumberFeature,fill=Level2)) +
    geom_bar(stat ="identity",na.rm = TRUE)+
    theme_bw()+
    xlab("")+
    ylab("Number of Feature")+
    theme(
      axis.text.x = element_text(angle = 80,hjust = 1,size=10), # Remove x axis tick labels
      axis.title=element_text(size=14,face="bold"),
      legend.title = element_text(face = "bold"),
      panel.grid=element_blank()
    )
  ggsave(filename = paste0(outdir,"/kegg.level3.MetabolismPathway.feature.png"),p2,height = 3200,width = deta_width_Pathway,units ="px", limitsize = FALSE)
  ggsave(filename = paste0(outdir,"/kegg.level3.MetabolismPathway.feature.pdf"),p2,height = 3200,width = deta_width_Pathway,units ="px", limitsize = FALSE)

  p3=ggplot(data = data_level3,aes(x = Pathway,y= NumberCompound,fill=Level2)) +
    geom_bar(stat ="identity",na.rm = TRUE)+
    theme_bw()+
    xlab("")+
    ylab("Number of Feature")+
    theme(
      axis.text.x = element_text(angle = 80,hjust = 1,size=10), # Remove x axis tick labels
      axis.title=element_text(size=14,face="bold"),
      legend.title = element_text(face = "bold"),
      panel.grid=element_blank()
    )
  ggsave(filename = paste0(outdir,"/kegg.level3.MetabolismPathway.png"),p3,height = 3200,width = deta_width_Pathway,units ="px", limitsize = FALSE)
  ggsave(filename = paste0(outdir,"/kegg.level3.MetabolismPathway.pdf"),p3,height = 3200,width = deta_width_Pathway,units ="px", limitsize = FALSE)
}
# draw_kegg_stat_bar_plot= function(cc,data_level3,outdir){
#   data_level2$Level2  = factor(data_level2$Level2,levels = data_level2$Level2)
#   p=ggplot(data = data_level2,aes(x = Level2,y= NumberFeature,fill=Level1)) +
#     geom_bar(stat ="identity",na.rm = TRUE)+
#     theme_bw()+
#     ylab("Number of Feature")+
#     xlab("")+
#     theme(
#       axis.text.x = element_text(angle = 80,hjust = 1,size=10), # Remove x axis tick labels
#       axis.title=element_text(size=14,face="bold"),
#       legend.title = element_text(face = "bold"),
#       panel.grid=element_blank()
#     )
#   ggsave(filename = paste0(outdir,"/kegg.level2.feature.png"),p,height = 3200,width = 4800,units ="px")
#   ggsave(filename = paste0(outdir,"/kegg.level2.feature.pdf"),p,height = 3200,width = 4800,units ="px")
#
#   p1 = p=ggplot(data = data_level2,aes(x = Level2,y= NumberCompound,fill=Level1)) +
#     geom_bar(stat ="identity",na.rm = TRUE)+
#     theme_bw()+
#     ylab("Number of Compound")+
#     xlab("")+
#     theme(
#       axis.text.x = element_text(angle = 80,hjust = 1,size=10), # Remove x axis tick labels
#       axis.title=element_text(size=14,face="bold"),
#       legend.title = element_text(face = "bold"),
#       panel.grid=element_blank()
#     )
#   ggsave(filename = paste0(outdir,"/kegg.level2.png"),p,height = 3200,width = 4800,units ="px")
#   ggsave(filename = paste0(outdir,"/kegg.level2.pdf"),p,height = 3200,width = 4800,units ="px")
#
#   data_level3$Pathway  = factor(data_level3$Pathway,levels = data_level3$Pathway)
#   p2=ggplot(data = data_level3,aes(x =Pathway,y= NumberFeature,fill=Level2)) +
#     geom_bar(stat ="identity",na.rm = TRUE)+
#     theme_bw()+
#     xlab("")+
#     ylab("Number of Feature")+
#     theme(
#       axis.text.x = element_text(angle = 80,hjust = 1,size=10), # Remove x axis tick labels
#       axis.title=element_text(size=14,face="bold"),
#       legend.title = element_text(face = "bold"),
#       panel.grid=element_blank()
#     )
#   ggsave(filename = paste0(outdir,"/kegg.level3.MetabolismPathway.feature.png"),p2,height = 3200,width = 4800,units ="px")
#   ggsave(filename = paste0(outdir,"/kegg.level3.MetabolismPathway.feature.pdf"),p2,height = 3200,width = 4800,units ="px")
#
#   p3=ggplot(data = data_level3,aes(x = Pathway,y= NumberCompound,fill=Level2)) +
#     geom_bar(stat ="identity",na.rm = TRUE)+
#     theme_bw()+
#     xlab("")+
#     ylab("Number of Feature")+
#     theme(
#       axis.text.x = element_text(angle = 80,hjust = 1,size=10), # Remove x axis tick labels
#       axis.title=element_text(size=14,face="bold"),
#       legend.title = element_text(face = "bold"),
#       panel.grid=element_blank()
#     )
#   ggsave(filename = paste0(outdir,"/kegg.level3.MetabolismPathway.png"),p3,height = 3200,width = 4800,units ="px")
#   ggsave(filename = paste0(outdir,"/kegg.level3.MetabolismPathway.pdf"),p3,height = 3200,width = 4800,units ="px")
# }

draw_kegg_level3_top20_plot = function(data_level3,outdir){
  sort_idx = order(data_level3$NumberFeature,decreasing = TRUE)
  data = data_level3[sort_idx,]
  data = data[c(1:20),]
  row.names(data) = c(1:length(data[,1]))
  p = ggplot(data=data) + geom_bar(aes(x =reorder(Pathway,-NumberFeature),y=NumberFeature),fill="#4D76B3",stat = "identity")+
    geom_text(aes(x =reorder(Pathway,-NumberFeature),y=NumberFeature,label=NumberFeature),hjust=-0.3,angle =90)+
    theme_bw()+
    xlab("")+
    ylim(c(0,max(data$NumberFeature)*1.1))+
    theme(
      axis.text.x = element_text(angle = 75,hjust = 1,size = 10,vjust = 1),
      panel.grid = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
    )
  ggsave(filename = paste0(outdir,"/kegg.level3.top20.png"), p)
  ggsave(filename = paste0(outdir,"/kegg.level3.top20.pdf"), p)
}


draw_kegg_enrichment_pathway_bubble_plot = function(data_level3,outdir){
  Level1= kegg_file$Level1
  Pvalue= kegg_file$Pvalue
  Qvalue= kegg_file$FDR
  deta_height = length(Qvalue)*80
  if(deta_height<5000){
    deta_height=5000
  }
  # Encrichment_factor = kegg_file$enrichment_factor
  # Enrichment_factor = runif(length(Pvalue),min=0,max=1)
  Enrichment_factor = kegg_file$NumberCompound / kegg_file$Background
  Pathway=kegg_file$Pathway
  Number = kegg_file$NumberCompound
  df = data.frame(Pathway,Enrichment_factor,Number,Pvalue,Qvalue,Level1)

  p = ggplot(data=df,aes(x = Enrichment_factor,y=Pathway,color=Pvalue,size=Number))+
    geom_point()+
    facet_grid(vars(Level1),scales="free",space="free")+
    theme_bw()+
    scale_color_gradient(low = "#FC0012",high = "#3500F8")+
    theme(panel.grid = element_blank(),
          strip.text.y = element_text(angle = 0,size =12),
          axis.text = element_text(size = 12)
    )
  ggsave(filename = paste0(outdir,"/kegg.level3.pvalue.png"),p,width = 4000,height = deta_height,units = "px", limitsize = FALSE)
  ggsave(filename = paste0(outdir,"/kegg.level3.pvalue.pdf"),p,width = 4000,height = deta_height,units = "px", limitsize = FALSE)

  p1 = ggplot(data=df,aes(x = Enrichment_factor,y=Pathway,color=Qvalue,size=Number))+
    geom_point()+
    facet_grid(vars(Level1),scales="free",space="free")+
    theme_bw()+
    scale_color_gradient(low = "#FC0012",high = "#3500F8")+
    theme(panel.grid = element_blank(),
          strip.text.y = element_text(angle = 0,size =12),
          axis.text = element_text(size = 12)
    )
  ggsave(filename = paste0(outdir,"/kegg.level3.qvalue.png"),p1,width = 4000,height = deta_height,units = "px", limitsize = FALSE)
  ggsave(filename = paste0(outdir,"/kegg.level3.qvalue.pdf"),p1,width = 4000,height = deta_height,units = "px", limitsize = FALSE)
}
# draw_kegg_enrichment_pathway_bubble_plot = function(kegg_file,outdir){
#   Level1= kegg_file$Level1
#   Pvalue= kegg_file$Pvalue
#   Qvalue= kegg_file$FDR
#   # Encrichment_factor = kegg_file$enrichment_factor
#   Enrichment_factor = runif(length(Pvalue),min=0,max=1)
#   Pathway=kegg_file$Pathway
#   Number = kegg_file$NumberCompound
#   df = data.frame(Pathway,Enrichment_factor,Number,Pvalue,Qvalue,Level1)
#
#   p = ggplot(data=df,aes(x = Enrichment_factor,y=Pathway,color=Pvalue,size=Number))+
#     geom_point()+
#     facet_grid(vars(Level1),scales="free",space="free")+
#     theme_bw()+
#     theme(panel.grid = element_blank(),
#           strip.text.y = element_text(angle = 0,size =12),
#           axis.text = element_text(size = 12)
#     )
#   ggsave(filename = paste0(outdir,"/kegg.level3.pvalue.png"),p,width = 4000,height = 5000,units = "px")
#   ggsave(filename = paste0(outdir,"/kegg.level3.pvalue.pdf"),p,width = 4000,height = 5000,units = "px")
#
#   p1 = ggplot(data=df,aes(x = Enrichment_factor,y=Pathway,color=Qvalue,size=Number))+
#     geom_point()+
#     facet_grid(vars(Level1),scales="free",space="free")+
#     theme_bw()+
#     theme(panel.grid = element_blank(),
#           strip.text.y = element_text(angle = 0,size =12),
#           axis.text = element_text(size = 12)
#     )
#   ggsave(filename = paste0(outdir,"/kegg.level3.qvalue.png"),p1,width = 4000,height = 5000,units = "px")
#   ggsave(filename = paste0(outdir,"/kegg.level3.qvalue.pdf"),p1,width = 4000,height = 5000,units = "px")
# }
#--------- main ----------

# data prepare
p <- arg_parser("section2 KEGG plot")
p <- add_argument(p, "--level2", help="kegg level2")
p <- add_argument(p, "--kegg", help="kegg pvalue")
p <- add_argument(p, "-o", "--output_dir", help="output dir")
argv <- parse_args(p)
# outdir = "G:/Lemonx/metabolomics/参考数据/2.MetaboliteIdentification/"
kegg.level2_file = read.csv(argv$level2)
kegg_file = read.csv(argv$kegg)
if(dim(kegg.level2_file)[1] > 0){
  # draw plot
  draw_kegg_stat_bar_plot(kegg.level2_file,kegg_file,argv$o)
  draw_kegg_level3_top20_plot(data_level3 = kegg_file,argv$o)
  draw_kegg_enrichment_pathway_bubble_plot(data_level3 =kegg_file,outdir = argv$o)
}


