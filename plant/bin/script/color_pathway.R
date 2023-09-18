suppressPackageStartupMessages(library(pathview))
suppressPackageStartupMessages(library(argparser))
suppressPackageStartupMessages(library(parallel))


p <- arg_parser("color kegg pathway")
p <- add_argument(p, "--kegg_result", help="kegg enrichment result")
p <- add_argument(p, "--kegg_dir", help="kegg png/kgml dir")
p <- add_argument(p, "--de_result", help="differential file")
p <- add_argument(p, "--out_dir", help="output dir")
p <- add_argument(p, "--thread", help="num of threads, default=1", default=1)
p <- add_argument(p, "--kgml", help="kgml path")

argv <- parse_args(p)

kegg_result <- read.csv(argv$kegg_result, stringsAsFactors=F)
if(!is.na(argv$de_result)){
  de_result <- read.csv(argv$de_result)
}else{
  de_result <- data.frame()
}
kegg_dir <- normalizePath(argv$kegg_dir)

setwd(argv$out_dir)

mclapply(1:nrow(kegg_result), function(i){
  tryCatch({
    pathway <- kegg_result[i, "KEGG"]
    pathway.id <- substr(pathway, 4, nchar(pathway))
    entries <- strsplit(kegg_result[i, "Compound"], ";")[[1]]
    cpd_data <- rep(0, length(entries))
    names(cpd_data) <- entries
    for(entry in entries){
      if(nrow(de_result[(de_result$KEGG == entry) & (de_result$regulated == "up"), ]) > 0){
        cpd_data[entry] = 1
      }else if(nrow(de_result[(de_result$KEGG == entry) & (de_result$regulated == "down"), ]) > 0){
        cpd_data[entry] = -1
      }
    }
    if(file.exists(paste0(kegg_dir, "/ko", pathway.id, ".png")) && file.size(paste0(kegg_dir, "/ko", pathway.id, ".png")) < 200000){
        pv.out <- pathview(pathway.id=pathway.id, kegg.dir=kegg_dir, cpd.data=cpd_data, new.signature=F, plot.col.key=F,
                       low = "blue", mid = "yellow", high="red", out.suffix="pathview", species="ko")
        file.rename(paste0("ko", pathway.id, ".pathview.png"), paste0(pathway, ".png"))
    }
  }, warning = function(w){
  }, error = function(e){
  },finally = {

  })


}, mc.cores=argv$thread)

# Intermediate file
intermediate_files <- Sys.glob("ko*.pathview.png")
file.remove(intermediate_files)
