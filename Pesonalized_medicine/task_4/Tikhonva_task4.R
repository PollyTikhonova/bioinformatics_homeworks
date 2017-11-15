library(genefilter)
library("hgu133plus2.db")
library("hgu133a.db")
load('datasets.Rdata')
library("clusterProfiler")
library("EnrichmentBrowser")
gs <- get.kegg.genesets("hsa")
cancer_pathways = c(gs[291:302],gs[304:306])
names(cancer_pathways)

find_cancers <- function(DATA, libre) {
  DATA.filtered <- nsFilter(DATA)[[1]]
  sampleinfo = pData(DATA)
  abnormalities <- factor(pData(DATA)[[1]], levels=c('c', 'd'), labels=c(0,1))
  diff_genes = rowttests(exprs(DATA.filtered), abnormalities)$p.value
  
  zonds = rownames(DATA.filtered)
  if (sum(p.adjust(diff_genes, method = "BH")<=0.1) == 0) {
    chosen_genes = zonds
  }
  else {
    chosen_genes = zonds[p.adjust(diff_genes, method = "BH")<=0.1]
  }
  
  x <- eval(parse(text=paste0(libre,'ENTREZID')))
  mapped_probes <- mappedkeys(x)
  gene_ENTREZID <- as.environment(as.list(x[mapped_probes]))
  DAT.ENTREZID <- sapply(chosen_genes, function(x){mget(x,gene_ENTREZID)})
  
  ncg <- enrichKEGG(DAT.ENTREZID,pvalueCutoff=0.2)

  sn=1
  cancers = c()
  for (ncgi in ncg$ID){
    if ((as.numeric(strsplit(ncgi, 'hsa')[[1]][2]) >= 5203) & (as.numeric(strsplit(ncgi, 'hsa')[[1]][2]) <= 5225)) {
      for (canpath in names(cancer_pathways)) {
        if (strsplit(canpath, '_')[[1]][1] == ncgi) {
          cancers = c(cancers, canpath)
          sn = sn+ 1
          break
        }
      }
    }
    
    if (sn == 2) {
      break
    }
  }
  return(cancers)
   
}

for (i in DAT) {
  print(c(find_cancers(i, i@annotation)))
}
