
s.w=lapply(w,function(x){(x-min(x))/(max(x)-min(x))})%>%as.data.frame()
s.w[is.na(s.w)] <- 0
for (i in seq(1, ncol(s.w))) {
  if(length(which((s.w[,i]==0)))>=20){
    print(i)
    del1 <- append(del1, -i)
  }
}
s.w <- s.w[,del1]

BIEC=matrix(nrow=n,ncol=n)
for (i in 1:n) {
  for(j in 1:n){
    BIEC[i,j]=sum(t(div[,i])*div[,j])*(sqrt(pre[i])*sqrt(pre[j]))^(-1)
  }
}
rownames(BIEC)=colnames(y)
colnames(BIEC)=colnames(y)

saver<- matrix(nrow=22538,ncol=500)
savenames<- matrix(nrow=22538,ncol=500)
genenames=dimnames(BIEC)[[2]]
n=15
test=BIEC
pt<- matrix(nrow=1,ncol=22538)
k=0
for (i in 1:22538) {
  for (j in i:22538){ 
    if (test[i,j]>0 && test[i,j]<1){
      pt[j]=(1-pt(test[i,j]*sqrt((n-2)/(1-test[i,j]^2)),n-2))*2
    }
    else if(test[i,j]<0){
      pt[j]=pt(test[i,j]*sqrt((n-2)/(1-test[i,j]^2)),n-2)*2
    }else{
      pt[j]=1e-10
    }
    
  }
  add=c()
  for(m in 1:22538){
    if(pt[m]<0.0001 && test[i,m]>=0.99){
      add=append(add,m)
    }
  }
  for(h in 1:length(add)){
    saver[i,h]= test[i,add[h]]
    savenames[i,h]=genenames[add[h]]
  }
  k=k+1
  print(k)
  rm(add) 
}

del2=c()
for (i in seq(1, nrow(saver))) {
  if(sum(saver[i,])==1){
    #print(i)
    del2 <- append(del2, -i)
  }
}
saver <- saver[del2,]
savenames= savenames[del2,]

sum=0
p<- matrix(nrow=1,ncol=dim(savenames)[1])
regulatory=matrix(nrow=22538*100,ncol=3)
for(i in 1:dim(savenames)[1]){
  p[i]=length(which(saver[i,]>0))
  if(i>=2){
    for(j in (sum+1):(sum+p[i])){
      regulatory[j,1]=savenames[i,1]
      regulatory[j,2]=savenames[i,j-sum+1] 
      regulatory[j,3]=saver[i,j-sum]
    }
  }else{
    for(j in 1:p[i]){
      regulatory[j,1]=savenames[i,1]
      regulatory[j,2]=savenames[i,j+1] 
      regulatory[j,3]=saver[i,j]
    }
  }
  sum=sum+p[i]
}



for(i in 1:dim(p)[2]){
  pv<- matrix(nrow=15,ncol=p[i])
  for(k in 1:p[i]){
    pv[,k]=s.w[,savenames[i,k]]
  }
  colnames(pv)=c(savenames[i,which(saver[i,]>0)])
  tune.model = e1071::tune.nnet(pv[,1] ~ ., data = pv, size = 5:30)
  model = tune.model$best.parameters
  res <- neuralnet(pv[,1]~.,data=pv,hidden=model,linear.output=T)
  fact=res$generalized.weights[[1]]
  
  factd=decostand(fact,"range",1)
  facta=apply(factd,2,sum)
  factsort=order(facta,decreasing=T)
  selectnumber=ceiling(p[i]/10)
  factuse=factsort[1:selectnumber]
  for(j in 1:length(factuse)){
    regulatory1[i,j]=savenames[i,factuse[j]+1]
  }
  rm(pv)
  rm(res)
}
rownames(regulatory1)=rownames(savenames)


num.genes <- dim(weight.matrix1)[1]
genes <- colnames(weight.matrix1)
matrix.length <- length(weight.matrix)
relatory=matrix(nrow=num.genes*100,ncol=2)
colnames(relatory)=c("from.gene","to.gene")

k=0

for(i in 1:num.genes){
  genesort=order(weight.matrix1[i,],decreasing=T)
  for(j in 1:num.genes){
    if(weight.matrix1[i,genesort[j]]>0.012){
      relatory[k+j,1]=genes[i]
      relatory[k+j,2]=genes[genesort[j]]
    }
    else{
      break
    }  
  } 
  k=k+j-1
  print(i)
}

relatory =na.omit(relatory)

init.igraph<-function(data,dir=F,rem.multi=T)
{
  labels<-union(unique(data[,1]),unique(data[,2]))
  ids<-1:length(labels);names(ids)<-labels
  from<-as.character(data[,1]);to<-as.character(data[,2])
  edges<-matrix(c(ids[from],ids[to]),nc=2)
  g<-graph.empty(directed = dir)
  g<-add.vertices(g,length(labels))
  V(g)$label=labels
  g<-add.edges(g,t(edges))
  if (rem.multi)
  {
    E(g)$weight<-count.multiple(g)
    g<-simplify(g,remove.multiple = TRUE,
                remove.loops = TRUE,edge.attr.comb = "mean")
  }
  g
}
g <- init.igraph(relatory,dir=F,rem.multi=T)
kdegree=degree(g,mode="total")
sweight=matrix(nrow = length(y),ncol=1)
for(i in 1:length(y)){
  sweight[i]=0
  generowsort=order(weight.matrix1[which(rownames(weight.matrix1)==y[i],arr.ind = TRUE),],decreasing=T)
  genecolsort=order(weight.matrix1[,which(rownames(weight.matrix1)==y[i],arr.ind = TRUE)],decreasing=T)
  for (j in 1:num.genes){
    if(weight.matrix1[which(rownames(weight.matrix1)==y[i],arr.ind = TRUE),generowsort[j]]>0.012){
      sweight[i]=sweight[i]+weight.matrix1[which(rownames(weight.matrix1)==y[i],arr.ind = TRUE),generowsort[j]]
    }
    else{
      break
    }  
  }
  for (j in 1:num.genes){
    if(weight.matrix1[genecolsort[j],which(colnames(weight.matrix1)==y[i],arr.ind = TRUE)]>0.012){
      sweight[i]=sweight[i]+weight.matrix1[genecolsort[j],which(rownames(weight.matrix1)==y[i],arr.ind = TRUE)]
    }
    else{
      break
    }  
  }
  print(i)
}
alpha=1.3
sim=matrix(nrow = length(y),ncol=1)
for (i in 1:length(y)) {
  sim[i]=sweight[i]^alpha*kdegree[i]^(1-alpha)
  
}

simsort=order(sim,decreasing=T)
y1=matrix(nrow = length(y),ncol=1)
for (i in 1:length(y)) {
  y1[i]=y[simsort[i]]
}


y=matrix(nrow=2*dim(relatory)[1],ncol=1)
for(i in 1:dim(relatory)[1]){
  y[i]=relatory[i,1]
  y[i+dim(relatory)[1]]=relatory[i,2]
}
y=y[!duplicated(y)]

rt = read.table('diff_gene_list.txt',header = T,sep = '\t')
y = rt$gene_id


eg <- bitr(y, 
           fromType="SYMBOL", 
           toType=c("ENTREZID","ENSEMBL",'SYMBOL'),
           OrgDb="org.Hs.eg.db")

# Run GO enrichment analysis 
go <- enrichGO(eg$ENTREZID, 
               OrgDb = org.Hs.eg.db, 
               ont='ALL',
               pAdjustMethod = 'BH',
               pvalueCutoff = 0.05, 
               qvalueCutoff = 0.2,
               keyType = 'ENTREZID',
               readable      = TRUE)
dim(go)


dim(go[go$ONTOLOGY=='BP',])
dim(go[go$ONTOLOGY=='CC',])
dim(go[go$ONTOLOGY=='MF',])


barplot(go,showCategory=20,drop=T)
dotplot(go,showCategory=20)


kegg <- enrichKEGG(eg$ENTREZID, 
                   organism = 'hsa',  ## hsa为人的简写，bta是牛的简写 
                   keyType = 'kegg', 
                   pvalueCutoff = 0.05,
                   pAdjustMethod = 'BH', 
                   minGSSize = 3,
                   maxGSSize = 500,
                   qvalueCutoff = 0.2,
                   use_internal_data = FALSE
                   readable      = TRUE)
head(kegg)


dotplot(kegg, showCategory=20) 
barplot(kegg,showCategory=20,drop=T) 
cnetplot(kegg, foldChange=geneList) 
heatplot(kegg) 