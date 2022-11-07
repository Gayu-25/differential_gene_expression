#Reading an input cancer data file
data<-read.csv("D:/MSc Bioinformatics/Sem3/Cancer genomics/liver cancer.csv",row.names = 1)
print(data)
view(data)
#Create a count per matrix
cpmatrix=data
for(i in 1:ncol(data)){
  cpmatrix[,i]=(data[,i]/sum(data[,i]))*1000000
}

#Calculate a log of cpm
logcpm=log2(cpmatrix+1)

summary(logcpm)

#Calculate a z score
library(matrixStats)
z_score = (logcpm - rowMeans(logcpm))/rowSds(as.matrix(logcpm))[row(logcpm)]
z_score[is.na(z_score)]=0
#Calculate variance using log 
variance = apply(logcpm, 1, var)
variance = sort(variance,decreasing = T)
top50 = variance[1:50]
pmat = z_score[names(top50),]

#Create a heatmap
library(ComplexHeatmap)
Heatmap(pmat)


#To identify genes which are differential in tumor vs control samples

mat=matrix(NA,ncol=4,nrow = nrow(logcpm))
rownames(mat)=rownames(logcpm)
colnames(mat)=c('meanTumor','meanControl','pvalue','log2FC')

for(i in 1:nrow(logcpm)){
  vector1 = as.numeric(logcpm[i, 1:5])
  
  vector2 = as.numeric(logcpm[i, 6:10 ])
  
  res=t.test(vector1, vector2, paired = F, alternative = "two.sided")
  mat[i,1]=res$estimate[[1]]
  mat[i,2]=res$estimate[[2]]
  mat[i,3]=res$p.value
  mat[i,4]=mat[i,1]-mat[i,2]
  
}

mat=as.data.frame(mat)
num=which(is.nan(mat$pvalue))
mat[num,'pvalue']=1

library(EnhancedVolcano)
EnhancedVolcano(mat,lab = rownames(mat),x = 'log2FC' ,y ='pvalue')                
                