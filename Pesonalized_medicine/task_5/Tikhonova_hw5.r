
library(fossil)
library(iClusterPlus)

load('brca.dat')

data.full = data
d.1 = data.full$mRNA
d.2 = data.full$miRNA
d.3 = data.full$proteomics
elements = list(list(d.1, d.2, d.3), 
                list(d.1, d.2, NULL), 
                list(d.2, d.3, NULL), 
                list(d.1, d.3, NULL), 
                list(d.1, NULL, NULL), 
                list(d.2, NULL, NULL), 
                list(d.3, NULL, NULL))

icluster = sapply(elements, function(x){
    fit.all=iClusterPlus(
                  dt1 = t(x[[1]]),
                  dt2 <- if (!is.null(x[[2]])) t(x[[2]]) else NULL,
                  dt3 <- if (!is.null(x[[3]])) t(x[[3]]) else NULL,
                  type=c("gaussian","gaussian","gaussian"),
                  lambda=c(0.04,0.61,0.90),K=2,maxiter=10)
    return(rand.index(fit.all$clusters, as.numeric(factor(Y, labels = c(1,2,3)))))
})

kmeans = sapply(elements, function(x){ 
    big.data = rbind(scale(x[[1]]), 
                     if (!is.null(x[[2]])) scale(x[[2]]) else NULL,
                     if (!is.null(x[[3]])) scale(x[[3]]) else NULL)
    big.data = scale(big.data)
    set.seed(105)
    kmeans.all.fit = kmeans(t(big.data),3)
    return(rand.index(kmeans.all.fit$cluster,as.numeric(factor(Y, labels = c(1,2,3)))))
})

hclust = sapply(elements, function(x){
    big.data = rbind(scale(x[[1]]), 
                     if (!is.null(x[[2]])) scale(x[[2]]) else NULL,
                     if (!is.null(x[[3]])) scale(x[[3]]) else NULL)
    big.data = scale(big.data)
    d <- dist(t(big.data), method = "euclidean") # distance matrix
    fit <- hclust(d, method="average") 
    groups <- cutree(fit, k=3)
    return(rand.index(groups,as.numeric(factor(Y, labels = c(1,2,3)))))
})

table = rbind(icluster, kmeans, hclust)
colnames(table) = c('All', 'mRNA&miRNA', 'miRNA&prot', 'mRNA&prot', 'mRNA', 'miRNA', 'prot')
.=c('','','')
best_result = sapply(list(icluster, kmeans, hclust), max)
best_type =  sapply(list(icluster, kmeans, hclust), function(x){colnames(table)[which(x==max(x))]})
cbind(table,., best_result, best_type)
cat('Замечание: Icluster и Kmeans с каждым запуском могут давать разные результаты.')

#install.packages("mixOmics")
library(mixOmics)
df = list(mRNA = t(d.1), miRNA = t(d.2), proteomics = t(d.3))
mixOm = mixOmics(df, Y=Y, ncomp=3)
plot(mixOm)
