conn_res_cell <- function(conn_res){
  new <- as.data.frame(t(sapply(1:nrow(conn_res), function(i) unlist(strsplit(as.character(conn_res$set[i]), "_")))), stringsAsFactors=FALSE)
  colnames(new) = c("drug", "cell")
  conn_res <- cbind(new, conn_res)
  list <- lapply(unique(conn_res$cell), function(i) conn_res[conn_res$cell == i, ! colnames(conn_res) %in% c("set")])
  names(list) <- unique(conn_res$cell)
  return(list)
}

drug_targets_enrich <- function(drugs=NULL, c_cmap, cell="MCF7", trend, top=20, type="GO", ont="MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 2, maxGSSize = 500){
  if(is.null(drugs)){
    if(class(c_cmap)=="CMAPResults"){
      c_cmap <- cmapTable(c_cmap)
    }
    if(!is.data.frame(c_cmap)){
      stop('cmap_res should be either "CMAPResults" class from gCMAP package or data.frame')
    }
    c_cmap_list <- conn_res_cell(c_cmap)
    c_cmap <- c_cmap_list[[cell]]
    
    if(trend=="up"){
      drugs <- as.character(c_cmap$drug[1:top])
    }
    if(trend=="down"){
      drugs <- as.character(rev(c_cmap$drug)[1:top])
    }
  }
  drugs <- unique(drugs)
  tmpfile <- tempfile(fileext=".csv")
  download.file("http://biocluster.ucr.edu/~yduan004/drugbankR/drug_links_5.0.7.csv", tmpfile, quiet = TRUE)
  db_drug_links <- read.delim(tmpfile, sep=",")
  unlink(tmpfile)
  
  dbid <- as.character(db_drug_links[tolower(db_drug_links$Name) %in% tolower(drugs),"DrugBank.ID"])
  library(drugbankR)
  targets <- queryDB(ids=dbid, type = "getTargets")
  gnset <- na.omit(unlist(lapply(targets$t_gn_sym, function(i) unlist(strsplit(i, split = ";")))))
  
  tmpfile <- tempfile(fileext=".txt")
  download.file("http://biocluster.ucr.edu/~yduan004/drugbankR/universe_targets_in_drugbank.txt", tmpfile, quiet = TRUE)
  universe <- readLines(tmpfile)
  unlink(tmpfile)
  
  suppressPackageStartupMessages(library(clusterProfiler2))
  library(org.Hs.eg.db)
  if(type=="GO"){
    ego <- clusterProfiler2::enrichGO(gene = gnset, universe = universe, OrgDb = org.Hs.eg.db, keytype = 'SYMBOL', ont = ont, pAdjustMethod = pAdjustMethod, pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, minGSSize = minGSSize, maxGSSize = maxGSSize)
    return(ego)
  }
  
  if(type=="KEGG"){
    gnset_entrez <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = gnset, keytype = "SYMBOL", columns = "ENTREZID")$ENTREZID)
    univ_entrez <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = universe, keytype = "SYMBOL", columns = "ENTREZID")$ENTREZID)
    kk <- clusterProfiler2::enrichKEGG(gene=gnset_entrez, organism = "hsa", keyType = "kegg", pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, universe=univ_entrez, minGSSize = minGSSize, maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff)
    return(kk)
  }
}
