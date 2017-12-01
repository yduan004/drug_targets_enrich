# drug_targets_enrich
Do drug targets GO/KEGG enrichment analysis by using weighted hypergeometric test

# Install packages

Download, extract `tar.gz`, and install `clusterProfiler2`, `DOSE2`, `drugbankR` package from source
```{r ins, eval=FALSE}
install.packages("clusterProfiler2", repos = NULL, type="source")
install.packages("DOSE2", repos = NULL, type="source")
if(! file.exists("./drugbankR_1.3.tar.gz")){
  download.file("http://biocluster.ucr.edu/~yduan004/drugbankR/drugbankR_1.3.tar.gz", "./drugbankR_1.3.tar.gz")
}
install.packages("drugbankR", repos = NULL, type="source")
```

# Do `drug_targets_enrich` analysis by providing connectivity result from `gCMAP` package.

`c_cmap`: connectivity result from `connectivity_score` function in `gCMAP` package

`cell`: select the connected drugs in specific cell

`trend`: setting as "up" means select the drugs that have similar effect as query signature

`top`: number, select top connected drugs to do targets enrichment analysis

`type`: "GO" or "KEGG"

`ont`: when type is "GO", select the ontology of enriched GO terms
```{r examp, eval=TRUE}
library(gCMAP); library(utils)
if(! file.exists("./degList.rds")){
  download.file("http://biocluster.ucr.edu/~yduan004/cmap02/degList.rds", "./degList.rds")
}

logMA <- readRDS("degList.rds")$logFC  # 12437 X 3478
cmap <- induceCMAPCollection(logMA, element="members", higher=1, lower=-1)
# members(cmap)[1:10,1:10]
query <- logMA[, "sirolimus_MCF7", drop=FALSE]
query <- CMAPCollection(query, signed=rep(TRUE, ncol(query)))
c_cmap <- connectivity_score(query, cmap, element="members")
cmapTable(c_cmap)[1:20,]

source("drug_targets_enrich.R")
c_ego <- drug_targets_enrich(drugs=NULL, c_cmap, cell="MCF7", trend="up", top=20, type="GO", ont="MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 2, maxGSSize = 500)
head(c_ego[,c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "geneID")],10)

c_kegg <- drug_targets_enrich(drugs=NULL, c_cmap, cell="MCF7", trend="up", top=20, type="KEGG", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 2, maxGSSize = 500)
head(c_kegg[,c("ID", "Description", "GeneRatio", "BgRatio", "pvalue")],10)
```

# Do `drug_targets_enrich` analysis by providing drug names.
```{r exp2, eval=TRUE}
test_drug <- c("chlorpromazine","thioridazine","fluphenazine","trifluoperazine","prochlorperazine")

test_ego <- drug_targets_enrich(drugs=test_drug, type="GO", ont="MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 2, maxGSSize = 500)
head(test_ego[,c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "geneID")],10)

test_kegg <- drug_targets_enrich(drugs=test_drug, type="KEGG", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 2, maxGSSize = 500)
head(test_kegg[,c("ID", "Description", "GeneRatio", "BgRatio", "pvalue")],10)
```
