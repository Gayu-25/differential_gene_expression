Differential Gene analysis: 

Differential expression analysis means taking the normalised read count data and performing statistical analysis to discover quantitative changes in expression levels between experimental groups.
Here we take normalized data logCPM matrix value.

vec1=gene[i,control]

vec2=gene[i,tumor]

log2cpm=mean(vec1) - mean(vec2)
