# It includes the pre-processing steps for raw data.
#------------------------------------
library(mvtnorm)
library(ks)
library(Matrix)
library(stringi)
#-----------------------------------------------------

# Unify the gene_id and retain 15 digits for all genes
data.all$gene_id = stri_sub(data.all$gene_id,1,15)

# the names of tissues
ALL.tissue.name = data.all$tissue
ALL.tissue.name = ALL.tissue.name[!duplicated(ALL.tissue.name)]

# Extract the genes we focus on
MD_137 = read.csv("gene_id.csv",header = F)

# Retain only the variable data for the MD_137 gene
F.tissue = list()
gene.name = list()
for (l in 1:length(ALL.tissue.name)) {
  
  F.tissue[[l]] = subset(data.all,tissue == ALL.tissue.name[l])
  F.tissue[[l]] = subset(F.tissue[[l]],gene_id %in% MD_137$V1)
  
  # remove the missing values
  F.tissue[[l]] = data.frame(t(na.omit(t(F.tissue[[l]]))))
  
  # extract the gene names for each tissues
  gene.name[[l]] = F.tissue[[l]]$gene_id
}

# Obtain the intersect genes
intersect = Reduce(intersect,gene.name)
for (l in 1:length(ALL.tissue.name)) {
  
  F.tissue[[l]] = subset(F.tissue[[l]],gene_id %in% intersect)
}

#--------------------------------------------------


# the location of response gene
chr_21 = subset(F.tissue[[1]],X.chr == "chr21")
response.candidate = chr_21$gene_id
r = 2
response.gene = response.candidate[r]
Y.id = which(intersect == response.gene)
Num.response = c()
for (l in 1:length(ALL.tissue.name)) {
  
  row = subset(F.tissue[[l]],gene_id == response.gene)
  Num.response[l] = length(which(row[1,] != "NA")) - 5
}

# get the remained tissue data after processing
delete.tissue = which(Num.response < 150)
remain.tissue = F.tissue[-delete.tissue]
remain.tissue.name = ALL.tissue.name[-delete.tissue]

# remove the irrelavant variables
for (l in 1:length(remain.tissue.name)) {
  
  remain.tissue[[l]] = remain.tissue[[l]][,-c(1:5)]
  remain.tissue[[l]] = t(remain.tissue[[l]])
  remain.tissue[[l]] = apply(remain.tissue[[l]],2,as.numeric)
  colnames(remain.tissue[[l]]) = intersect
}

# Extract the target data and souce data, respectively.
p = length(intersect) - 1
target.list = remain.tissue[7:15]
target.tissue.name = remain.tissue.name[7:15]     
source.list = remain.tissue[-c(7:15)]
source.tissue.name = remain.tissue.name[-c(7:15)]



