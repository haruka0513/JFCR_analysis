#JFCR39

#ヒートマップ作成

#サンプルの中のどれかで2倍以上発現しているIDを格納

data <- read.table("JFCR39_mas5.txt",header=T,sep="\t")

row_samp <- 2:130
nsamp <- length(row_samp)
gid <- as.character(data[,1])
gf_gid_list <- c()
hy_gid_list <- c()
th.rat <- 2
th.val <- 300
sname <- substr(colnames(data), 9, nchar(colnames(data))-21)

seq1 <- seq(2, nsamp+1, by = 3)
seq2 <- seq(3, nsamp+1, by = 3)
seq3 <- seq(4, nsamp+1, by = 3)

#gfについて3倍以上/1/3以下発現しているサンプルが一つでもあればgid_listに追加

for(i in 1:43){
	
	gf_rat <- data[,seq1[i]]/data[,seq3[i]]	
	gf_up_sign <- gf_rat >= th.rat & data[,seq1[i]] >= th.val
	gf_down_sign <- gf_rat <= 1/th.rat & data[,seq3[i]] >= th.val
	gf_gid_list <- c(gf_gid_list,gid[gf_up_sign],gid[gf_down_sign])
		
}

#hypoについて3倍以上発現しているサンプルが一つでもあればgid_listに追加

for(i in 1:43){
	
	hy_rat <- data[,seq2[i]]/data[,seq3[i]]		
	hy_up_sign <- hy_rat >=th.rat & data[,seq2[i]] >= th.val
	hy_down_sign <- hy_rat <= 1/th.rat & data[,seq3[i]] >= th.val
	hy_gid_list <- c(hy_gid_list,gid[hy_up_sign],gid[hy_down_sign])
		
}

gf_gid_list_uniq <- unique(gf_gid_list)
hy_gid_list_uniq <- unique(hy_gid_list)

#元データから、サンプルの中のどれかで3倍以上発現しているIDの行のみを抽出

data_useful <- read.table("JFCR39_mas5.txt",header=T,row.names=1,sep="\t")

gf_data <- data_useful[gf_gid_list_uniq,]
hy_data <- data_useful[hy_gid_list_uniq,]


#ratio計算

gf_ngene <- nrow(gf_data)
hy_ngene <- nrow(hy_data)
row_samp <- 1:129
nsamp <- length(row_samp)
gf_ratio_mat <- matrix(0,gf_ngene,nsamp/3)
hy_ratio_mat <- matrix(0,hy_ngene,nsamp/3)

seq1 <- seq(1, nsamp, by = 3)
seq2 <- seq(2, nsamp, by = 3)
seq3 <- seq(3, nsamp, by = 3)

for(i in 1:43){
	
	ct <- pmax(gf_data[,seq3[i]],50)
	x <- pmax(gf_data[,seq1[i]],50)
	gf_ratio_mat[,i] <- log(x/ct)
	
}

for(i in 1:43){
	
	ct <- pmax(hy_data[,seq3[i]],50)
	x <- pmax(hy_data[,seq2[i]],50)
	hy_ratio_mat[,i] <- log(x/ct)
	
}

##seq1を2に変更してみた↑

colnames(gf_ratio_mat) <- sname[seq2]
rownames(gf_ratio_mat) <- gf_gid_list_uniq
colnames(hy_ratio_mat) <- sname[seq3]
rownames(hy_ratio_mat) <- hy_gid_list_uniq

#write.table(sname[seq2],"sname.txt",append=T,row.names=F,col.names=F,quote=F,sep="\t")

#gf_ratio_mat_1 <- matrix(0,100,nsamp/3)
#gf_ratio_mat_1 <- gf_ratio_mat_sort[1:70,]
#gf_ratio_mat_2 <- matrix(0,100,nsamp/3)
#gf_ratio_mat_2 <- gf_ratio_mat_sort[71:140,]
#gf_ratio_mat_3 <- matrix(0,100,nsamp/3)
#gf_ratio_mat_3 <- gf_ratio_mat_sort[141:210,]
#gf_ratio_mat_4 <- matrix(0,100,nsamp/3)
#gf_ratio_mat_4 <- gf_ratio_mat_sort[211:290,]
#gf_ratio_mat_5 <- matrix(0,100,nsamp/3)
#gf_ratio_mat_5 <- gf_ratio_mat_sort[291:370,]


#sort　遺伝子もクラスタリングするなら使用しない
#order関数は、並べ替えて元の列番号を返す関数
#gf_od <- order(gf_ratio_mat[,1],decreasing=T)
#gf_ratio_mat_sort <- gf_ratio_mat[gf_od,]
#hy_od <- order(hy_ratio_mat[,1],decreasing=T)
#hy_ratio_mat_sort <- hy_ratio_mat[hy_od,]


# ヒートマップ表示
install.packages("gplots")
#ターミナルで確認できる→/var/folders/2b/yfm8lh4s4pxdz1jfy_nj7yzm0000gq/T//RtmpMDSlV9/downloaded_packages 
library(gplots)
＃pはピアソンの相関の定義
dist.p <- function(x) as.dist((1-cor(t(x)))/2)
dist.s <- function(x) as.dist((1-cor(t(x), method="spearman"))/2)
#hclust.ward <- function(x) hclust(x, "ward")
hclust.ward <- function(x) hclust(x, "ward.D")
col1 <- heat.colors(60)

png(file="test-heatmap.png", width=600, height=600)

#Rowv = FALSE  dendrogramにboth,rowを指定した時にTRUEにする必要あり

hv <- heatmap.2(hy_ratio_mat, Rowv=T, Colv=T, distfun=dist.p, hclustfun=hclust.ward, dendrogram="both", margin=c(9.5,10), breaks=seq(-3,3,0.1), col=col1, scale="none", key=TRUE, keysize=1, symkey=FALSE, density.info="none", trace="none", cexCol=1)


#保存する        
dev.off()     











#Cmap
#Cmap用に、GF,Hypoそれぞれのratio_2の結果に対して、Cmapに含まれるプローブidのみを抜き出したい

#data <- read.table("test.txt",header=T,row.names=1)
#cmap_gene_data <- read.table("test_geneid.txt")
#texteditよりCoteditで編集した方がいい

data <- read.table("0417_JFCR39_mas5_gene_GF_up_array.txt",header=T,row.names=1,sep="\t")
cmap_gene_data <- read.table("U133A-psid.txt")

cmap_gene_list <- c()
cmap_gene_list <- as.character(cmap_gene_data[,1])

pick_up_gene <- rownames(data[cmap_gene_list,])
NA_id <- grep("NA",pick_up_gene)
#pick_up_geneベクトルの中でNAとなっている要素番号
pick_up_gene_2 <- pick_up_gene[-NA_id]
#要素の削除
write.table(t(pick_up_gene_2),"test_pick_up_gene.txt",append=T,sep="\n",row.names=F,col.names=F,quote=F)










#DAVID
#U251で発現が2倍以上・1/2以下の遺伝子名を出力する

data <- read.table("JFCR39_mas5.txt",header=T,sep="\t")

row_samp <- 2:130
nsamp <- length(row_samp)
gid <- as.character(data[,1])

th.rat <- 2
th.val <- 300

seq1 <- seq(2,nsamp+1,by=3)
seq2 <- seq(3,nsamp+1,by=3)
seq3 <- seq(4,nsamp+1,by=3)

gf_rat <- data[,seq1[39]]/data[,seq3[39]]
gf_up_sign <- gf_rat >= th.rat & data[,seq1[39]] >= th.val
gf_down_sign <- gf_rat <= 1/th.rat & data[,seq3[39]] >= th.val
#FALSE/TRUE

write.table(gid[gf_up_sign],"u251_up.txt",append=T,sep="\t",row.names=F,col.names=F,quote=F)
write.table(gid[gf_down_sign],"u251_down.txt",append=T,sep="\t",row.names=F,col.names=F,quote=F)

#SKOV3で発現が2倍以上・1/2以下の遺伝子名を出力する

s_gf_rat <- data[,seq1[42]]/data[,seq3[42]]
s_gf_up_sign <- gf_rat >= th.rat & data[,seq1[42]] >= th.val
s_gf_down_sign <- gf_rat <= 1/th.rat & data[,seq3[42]] >= th.val
#FALSE/TRUE

write.table(gid[s_gf_up_sign],"skov3_up.txt",append=T,sep="\t",row.names=F,col.names=F,quote=F)
write.table(gid[s_gf_down_sign],"skov3_down.txt",append=T,sep="\t",row.names=F,col.names=F,quote=F)


#A群全てで発現が2倍以上・1/2以下の遺伝子名を出力する

data <- read.table("JFCR39_mas5.txt",header=T,sep="\t")

row_samp <- 2:130
nsamp <- length(row_samp)
gid <- as.character(data[,1])
n_gid <- length(gid)
gid_list <- c()

th.rat <- 2
th.val <- 300

seq1 <- seq(2,nsamp+1,by=3)
seq2 <- seq(3,nsamp+1,by=3)
seq3 <- seq(4,nsamp+1,by=3)

#A群
a <- c(41,27,32,5,6,1,38,35,39)
n_a <- length(a)
count <- 0

for(j in 1:n_gid){
for(i in 1:n_a){

gf_rat <- data[j,seq1[a[i]]]/data[j,seq3[a[i]]]

if(gf_rat >= th.rat & data[j,seq1[a[i]]] >= th.val){
count = count + 1
}

}
if(count==n_a){
gid_list <- c(gid_list,gid[j])
count <- 0
}
}

#A群のいずれかで発現が2倍以上・1/2以下の遺伝子名を出力する

gf_gid_list <- c()
gf_gid_list_2 <- c()

for(i in 1:n_a){

gf_rat <- data[,seq1[a[i]]]/data[,seq3[a[i]]]
gf_up_sign <- gf_rat >= th.rat & data[,seq1[a[i]]] >= th.val
gf_down_sign <- gf_rat <= 1/th.rat & data[,seq3[a[i]]] >= th.val
gf_gid_list <- c(gf_gid_list,gid[gf_up_sign])
gf_gid_list_2 <- c(gf_gid_list_2,gid[gf_down_sign])

}

gf_gid_list_uniq <- unique(gf_gid_list)
gf_gid_list_2_uniq <- unique(gf_gid_list_2)

write.table(gf_gid_list_uniq,"a_group_up.txt",append=T,sep="\t",row.names=F,col.names=F,quote=F)
write.table(gf_gid_list_2_uniq,"a_group_down.txt",append=T,sep="\t",row.names=F,col.names=F,quote=F)


#Cmap for group_a
#group_aで2倍以上・1/2以下発現しているプローブについてCmapに含まれるプローブidのみを抜き出したい

data <- read.table("a_group_up.txt",row.names=1)
cmap_gene_data <- read.table("U133A-psid.txt")

cmap_gene_list <- c()
cmap_gene_list <- as.character(cmap_gene_data[,1])

pick_up_gene <- rownames(data[cmap_gene_list,])
NA_id <- grep("NA",pick_up_gene)
#pick_up_geneベクトルの中でNAとなっている要素番号
pick_up_gene_2 <- pick_up_gene[-NA_id]
#要素の削除
write.table(t(pick_up_gene_2),"test_pick_up_gene.txt",append=T,sep="\n",row.names=F,col.names=F,quote=F)











#主成分分析
#GFのレシオのみを出力する
data <- read.table("JFCR39_mas5.txt",header=T,sep="\t",row.names=1)

ngene <- nrow(data)
nsamp <- 43*3
seq1 <- seq(1,nsamp,by=3)
seq2 <- seq(2,nsamp,by=3)
seq3 <- seq(3,nsamp,by=3)

output <- matrix(0,ngene,nsamp/3)

for(i in 1:43){
  
  rat <- data[,seq1[i]]/data[,seq3[i]]
  output[,i] <- rat
  
}

write.table(output,"jfcr_ratio_for_pc.txt",append=T,sep="\t",row.names=F,col.names=F,quote=F)


#分析
data <- read.table("jfcr_ratio_for_pc.txt",sep="\t")
res <- prcomp(t(data))

#主成分分析の結果として標準偏差を返す
summary(res)

#PC1,PC2でプロット
#主成分得点は$xに列単位で記録される
res$x

#プロットするときのためにサンプルめいを格納しておく
data <- read.table("JFCR39_mas5.txt",header=T,sep="\t")
sname <- substr(colnames(data), 9, nchar(colnames(data))-21)
row_samp <- 2:130
nsamp <- length(row_samp)
seq1 <- seq(2, nsamp+1, by = 3)

plot(res$x[,1],res$x[,2],type="n")
text(res$x[,1],res$x[,2],sname[seq1],cex=0.5)

plot(res$x[,1],res$x[,2],xlim=c(-100,100),ylim=c(-100,100),type="n")
text(res$x[,1],res$x[,2],sname[seq1],cex=0.5)

#各合成変数における主成分（係数）を表示（43×54000）
#主成分は$rotationに列を単位に記録される
res$rotation
#第一主成分のみ表示したい場合
res$rotation[,1]

#プロットするときのために遺伝子名を格納しておく
data <- read.table("0414_JFCR39_mas5_gene.txt",header=T,sep="\t",quote="\"")
gene <- as.character(data[,2])

#///で区切られて遺伝子が2つ以上記載されているプローブがあってプロットするとみにくいから消す
gene2 <- gsub("///","",gene)

plot(res$rotation[,1],res$rotation[,2],type="n")
text(res$rotation[,1],res$rotation[,2],gene2,cex=0.5)

#結合変数の組み合わせをかえてみる
plot(res$x[,1],res$x[,3],type="n")
text(res$x[,1],res$x[,3],sname[seq1],cex=0.5)

plot(res$x[,1],res$x[,4],type="n")
text(res$x[,1],res$x[,4],sname[seq1],cex=0.5)

plot(res$x[,1],res$x[,5],type="n")
text(res$x[,1],res$x[,5],sname[seq1],cex=0.5)

plot(res$x[,2],res$x[,3],type="n")
text(res$x[,2],res$x[,3],sname[seq1],cex=0.5)

plot(res$x[,2],res$x[,4],type="n")
text(res$x[,2],res$x[,4],sname[seq1],cex=0.5)

plot(res$x[,2],res$x[,5],type="n")
text(res$x[,2],res$x[,5],sname[seq1],cex=0.5)

plot(res$x[,3],res$x[,4],type="n")
text(res$x[,3],res$x[,4],sname[seq1],cex=0.5)

plot(res$x[,3],res$x[,5],type="n")
text(res$x[,3],res$x[,5],sname[seq1],cex=0.5)



#GFのレシオ、4サンプル以上で発現2倍以上/1/2以下のプローブのみを出力する
data <- read.table("JFCR39_mas5.txt",header=T,sep="\t")

ngene <- nrow(data)
nsamp <- ncol(data)-1
nsamp3 <- nsamp/3

seq1 <- seq(2,nsamp+1,by=3)
seq2 <- seq(3,nsamp+1,by=3)
seq3 <- seq(4,nsamp+1,by=3)

up_count <- 0
down_count <- 0
th.rat <- 2
th.count <- 20

output <- c()

for(i in 1:ngene){
	
	for(j in 1:nsamp3){
		
		rat <- pmax(data[i,seq1[j]],50)/pmax(data[i,seq3[j]],50)
			
			if(rat >= th.rat){		
				up_count = up_count + 1						
			}
			if(rat <= 1/th.rat){				
				down_count = down_count + 1						
			}
	}	

	if(up_count >= th.count | down_count >= th.count){
		
		output <- c(output,as.character(data[i,1]))
		
	}
	
	up_count <- 0
	down_count <- 0
			
}

write.table(output,"jfcr_gene_for_pc_20.txt",append=T,sep="\t",row.names=F,col.names=F,quote=F)

#出力したプローブリストより、その行のレシオを出力
data <- read.table("JFCR39_mas5.txt",header=T,sep="\t",row.names=1)
gene <- read.table("jfcr_gene_for_pc_20.txt",header=F,sep="\t")

gene_list <- as.character(gene[,1])

nsamp <- ncol(data)
nsamp3 <- nsamp/3
ngene <- length(gene_list)

seq1 <- seq(1,nsamp,by=3)
seq2 <- seq(2,nsamp,by=3)
seq3 <- seq(3,nsamp,by=3)

output <- matrix(0,ngene,nsamp/3)

for(i in 1:nsamp3){

rat <- log(data[gene_list,seq1[i]]/data[gene_list,seq3[i]])
output[,i] <- rat

}

write.table(output,"jfcr_ratio_for_pc_20_log.txt",append=T,sep="\t",row.names=F,col.names=F,quote=F)

#分析

data <- read.table("jfcr_ratio_for_pc_15_log.txt",sep="\t")
res <- prcomp(t(data))

#サンプル名格納
data <- read.table("JFCR39_mas5.txt",header=T,sep="\t")
sname <- substr(colnames(data), 9, nchar(colnames(data))-21)
row_samp <- 2:130
nsamp <- length(row_samp)
seq1 <- seq(2, nsamp+1, by = 3)

plot(res$x[,1],res$x[,2],type="n")
text(res$x[,1],res$x[,2],sname[seq1],cex=0.5)

#主成分分析の結果として標準偏差を返す
summary(res)

#プロットするときのために遺伝子名を格納しておく
data <- read.table("0414_JFCR39_mas5_gene.txt",header=T,sep="\t",quote="\"",row.names=1)
gene <- as.character(data[gene_list,1])

#///で区切られて遺伝子が2つ以上記載されているプローブがあってプロットするとみにくいから消す
gene2 <- gsub("///","",gene)

plot(res$rotation[,1],res$rotation[,2],type="n")
text(res$rotation[,1],res$rotation[,2],gene2,cex=0.5)









#GFのレシオとHypoのレシオ、Xサンプル以上で発現2倍以上/1/2以下のプローブのみを出力する
data <- read.table("JFCR39_mas5.txt",header=T,sep="\t")

ngene <- nrow(data)
nsamp <- ncol(data)-1
nsamp3 <- nsamp/3

seq1 <- seq(2,nsamp+1,by=3)
seq2 <- seq(3,nsamp+1,by=3)
seq3 <- seq(4,nsamp+1,by=3)

up_count <- 0
down_count <- 0
th.rat <- 2
th.count <- 10

output <- c()

for(i in 1:ngene){
	
	for(j in 1:nsamp3){
		
		rat_gf <- pmax(data[i,seq1[j]],50)/pmax(data[i,seq3[j]],50)
		rat_hy <- pmax(data[i,seq2[j]],50)/pmax(data[i,seq3[j]],50)
					
			if(rat_gf >= th.rat | rat_hy >= th.rat){	
			if(rat_gf >= th.rat & rat_hy >= th.rat){		
				up_count = up_count + 2					
			}else{						
				up_count = up_count + 1										
			}
			}
			if(rat_gf <= 1/th.rat | rat_hy <=1/th.rat){	
			if(rat_gf <= 1/th.rat & rat_hy <=1/th.rat){				
				#down_count = down_count + 2						
			}else{							
				#down_count = down_count + 1						
			}
			}		
	}	

	if(up_count >= th.count | down_count >= th.count){
		
		output <- c(output,as.character(data[i,1]))
		
	}
	
	up_count <- 0
	down_count <- 0
			
}

write.table(output,"jfcr_gene_for_pc_10_gf_hy_uponly.txt",append=T,sep="\t",row.names=F,col.names=F,quote=F)

#出力したプローブリストより、その行のレシオを出力
data <- read.table("JFCR39_mas5.txt",header=T,sep="\t",row.names=1)
gene <- read.table("jfcr_gene_for_pc_10_gf_hy_uponly.txt",header=F,sep="\t")

gene_list <- as.character(gene[,1])

nsamp <- ncol(data)
nsamp3 <- nsamp/3
ngene <- length(gene_list)

seq1 <- seq(1,nsamp,by=3)
seq2 <- seq(2,nsamp,by=3)
seq3 <- seq(3,nsamp,by=3)

output <- matrix(0,ngene,(nsamp3)*2)

for(i in 1:nsamp3){

rat_gf <- log(data[gene_list,seq1[i]]/data[gene_list,seq3[i]])
rat_hy <- log(data[gene_list,seq2[i]]/data[gene_list,seq3[i]])

output[,i] <- rat_gf
output[,nsamp3+i] <- rat_hy

}

write.table(output,"jfcr_ratio_for_pc_10_log_gf_hy_uponly.txt",append=T,sep="\t",row.names=F,col.names=F,quote=F)

#分析

data <- read.table("jfcr_ratio_for_pc_5_39_50_log.txt",sep="\t")
res <- prcomp(t(data))

#サンプル名格納
data <- read.table("JFCR39_mas5.txt",header=T,sep="\t")
sname <- substr(colnames(data), 9, nchar(colnames(data))-21)
row_samp <- 2:130
nsamp <- length(row_samp)
nsamp3 <- nsamp/3
seq1 <- seq(2, nsamp+1, by = 3)
seq2 <- seq(3, nsamp+1, by = 3)

plot(res$x[,1],res$x[,2],type="n")
text(res$x[1:nsamp3,1],res$x[1:nsamp3,2],sname[seq1],cex=0.5)
text(res$x[(nsamp3+1):(nsamp3*2),1],res$x[(nsamp3+1):(nsamp3*2),2],sname[seq2],cex=0.5,col=2)


#主成分分析の結果として標準偏差を返す
summary(res)

#プロットするときのために遺伝子名を格納しておく
data <- read.table("0414_JFCR39_mas5_gene.txt",header=T,sep="\t",quote="\"",row.names=1)
gene <- as.character(data[gene_list,1])

#///で区切られて遺伝子が2つ以上記載されているプローブがあってプロットするとみにくいから消す
gene2 <- gsub("///","",gene)

plot(res$rotation[,1],res$rotation[,2],type="n")
text(res$rotation[,1],res$rotation[,2],gene2,cex=0.5)


#hypo,gf両方で10/86 or 20/86 サンプル2倍以上1/2以下変動している遺伝子群についてのヒートマップ0723

gene_list <- read.table("jfcr_gene_for_pc_10_gf_hy.txt")
data_useful <- read.table("JFCR39_mas5.txt",header=T,row.names=1,sep="\t")
gene <- as.character(gene_list[,1])

data <- data_useful[gene,]

#ratio計算

ngene <- nrow(data)
row_samp <- 1:129
nsamp <- length(row_samp)
gf_ratio_mat <- matrix(0,ngene,nsamp/3)
hy_ratio_mat <- matrix(0,ngene,nsamp/3)

seq1 <- seq(1, nsamp, by = 3)
seq2 <- seq(2, nsamp, by = 3)
seq3 <- seq(3, nsamp, by = 3)

for(i in 1:43){
	
	ct <- pmax(data[,seq3[i]],50)
	gf <- pmax(data[,seq1[i]],50)
	hy <- pmax(data[,seq2[i]],50)
	
	gf_ratio_mat[,i] <- log(gf/ct)
	hy_ratio_mat[,i] <- log(hy/ct)	
}


colnames(gf_ratio_mat) <- sname[seq2]
rownames(gf_ratio_mat) <- gene
colnames(hy_ratio_mat) <- sname[seq3]
rownames(hy_ratio_mat) <- gene

#write.table(gf_ratio_mat,"test.txt",append=T,sep="\t",row.names=F,col.names=F,quote=F)

#ヒートマップ表示
install.packages("gplots")
#ターミナルで確認できる→/var/folders/2b/yfm8lh4s4pxdz1jfy_nj7yzm0000gq/T//RtmpMDSlV9/downloaded_packages 
library(gplots)
＃pはピアソンの相関の定義
dist.p <- function(x) as.dist((1-cor(t(x)))/2)
dist.s <- function(x) as.dist((1-cor(t(x), method="spearman"))/2)
#hclust.ward <- function(x) hclust(x, "ward")
hclust.ward <- function(x) hclust(x, "ward.D")
col1 <- heat.colors(60)

png(file="0723_heatmap_agree_hy.png", width=600, height=600)

#Rowv = FALSE  dendrogramにboth,rowを指定した時にTRUEにする必要あり

hv <- heatmap.2(hy_ratio_mat, Rowv=T, Colv=T, distfun=dist.p, hclustfun=hclust.ward, dendrogram="both", margin=c(9.5,10), breaks=seq(-3,3,0.1), col=col1, scale="none", key=TRUE, keysize=1, symkey=FALSE, density.info="none", trace="none", cexCol=1)


#保存する        
dev.off()     








#GFのレシオ、5~10/43サンプルで発現2倍以上/1/2以下のプローブのみを出力する0729
	data <- read.table("JFCR39_mas5_39.txt",header=T,sep="\t")
	
	ngene <- nrow(data)
	nsamp <- ncol(data)-1
	nsamp3 <- nsamp/3
	
	seq1 <- seq(2,nsamp+1,by=3)
	seq2 <- seq(3,nsamp+1,by=3)
	seq3 <- seq(4,nsamp+1,by=3)
	
	up_count <- 0
	down_count <- 0
	th.rat <- 2
	th.count <- 11
	
	output <- c()
	
	for(i in 1:ngene){
		
		for(j in 1:nsamp3){
			
			rat <- pmax(data[i,seq1[j]],50)/pmax(data[i,seq3[j]],50)
			#rat <- data[i,seq1[j]]/data[i,seq3[j]]
						
				if(rat >= th.rat){		
					up_count = up_count + 1						
				}
				if(rat <= 1/th.rat){				
					down_count = down_count + 1						
				}
		}	
	
		if(up_count >= th.count | down_count >= th.count){
			
			output <- c(output,as.character(data[i,1]))
			
		}
		
		up_count <- 0
		down_count <- 0
				
	}
	
	write.table(output,"jfcr_gene_for_pc_11_39_50.txt",append=T,sep="\t",row.names=F,col.names=F,quote=F)
	
	#出力したプローブリストより、その行のレシオを出力
	data <- read.table("JFCR39_mas5_39.txt",header=T,sep="\t",row.names=1)
	gene <- read.table("jfcr_gene_for_pc_11_39_50.txt",header=F,sep="\t")
	
	gene_list <- as.character(gene[,1])
	
	nsamp <- ncol(data)
	nsamp3 <- nsamp/3
	ngene <- length(gene_list)
	
	seq1 <- seq(1,nsamp,by=3)
	seq2 <- seq(2,nsamp,by=3)
	seq3 <- seq(3,nsamp,by=3)
	
	output <- matrix(0,ngene,nsamp/3)
	
	for(i in 1:nsamp3){
	
	rat <- log(pmax(data[gene_list,seq1[i]],50)/pmax(data[gene_list,seq3[i]],50))
	#rat <- log(data[gene_list,seq1[i]]/data[gene_list,seq3[i]])
	
	output[,i] <- rat
	
	}
	
	write.table(output,"jfcr_ratio_for_pc_11_39_50_log.txt",append=T,sep="\t",row.names=F,col.names=F,quote=F)
	
	
	
	
#分析
	
data <- read.table("jfcr_ratio_for_pc_5_39_50_log.txt",sep="\t")
res <- prcomp(t(data))
	
#サンプル名格納
data <- read.table("JFCR39_mas5_39.txt",header=T,sep="\t")
	
nsamp <- 117
	
seq1 <- seq(2,nsamp+1,by=3)
seq2 <- seq(3,nsamp+1,by=3)
seq3 <- seq(4,nsamp+1,by=3)
	
colname_gf <- colnames(data[seq1])
colname_n <- colnames(data[seq3])
	
sname <- substr(colnames(data), 9, nchar(colnames(data))-21)
row_samp <- 2:118
nsamp <- length(row_samp)
seq1 <- seq(2, nsamp+1, by = 3)
	
plot(res$x[,1],res$x[,2],type="n")
text(res$x[,1],res$x[,2],sname[seq1],cex=0.5)
	
#主成分分析の結果として標準偏差を返す
summary(res)

#プロットするときのために遺伝子名を格納しておく
data <- read.table("0414_JFCR39_mas5_gene.txt",header=T,sep="\t",quote="\"",row.names=1)
gene <- read.table("jfcr_gene_for_pc_5_39_50.txt",header=F,sep="\t")	
gene_list <- as.character(gene[,1])
gene <- as.character(data[gene_list,1])

#///で区切られて遺伝子が2つ以上記載されているプローブがあってプロットするとみにくいから消す
gene2 <- gsub("///","",gene)

plot(res$rotation[,1],res$rotation[,2],type="n")
text(res$rotation[,1],res$rotation[,2],gene2,cex=0.5)




#3group 分類 for GSEA

x1 <- res$rotation[,1]
x2 <- res$rotation[,2]

ngene <- length(gene_list)

group1_gene <- c()
group2_gene <- c()
group3_gene <- c()


for(i in 1:ngene){
	
	if(x1[i] < 0){
		
		group1_gene <- c(group1_gene,gene_list[i])
		
	}
	
	if(x1[i] >= 0 && x2[i] >= 0){
		
		group2_gene <- c(group2_gene,gene_list[i])
		
	}
	
	if(x1[i] >= 0 && x2[i] <= 0){
		
		group3_gene <- c(group3_gene,gene_list[i])
		
	}		
	
}

x1_samp <- res$x[,1]
x2_samp <- res$x[,2]

nsamp <- length(colname_gf)

group1_samp <- c()
group2_samp <- c()
group3_samp <- c()

for(i in 1:nsamp){
	
	if(x1_samp[i] < 0){
		
		group1_samp <- c(group1_samp,colname_gf[i])
		group1_samp <- c(group1_samp,colname_n[i])
		
	}
	if(x1_samp[i] >= 0 && x2_samp[i] >= 0){
		
		group2_samp <- c(group2_samp,colname_gf[i])
		group2_samp <- c(group2_samp,colname_n[i])
				
	}
	if(x1_samp[i] >= 0 && x2_samp[i] <= 0){
		
		group3_samp <- c(group3_samp,colname_gf[i])
		group3_samp <- c(group3_samp,colname_n[i])
				
	}
			
}

data <- read.table("JFCR39_mas5_39.txt",header=T,row.names=1,sep="\t")

#data_group1 <- as.list(NA)で初期化できます
data_group1_c <- data[group1_gene,group1_samp]
colnames(data_group1) <- group1_samp
rownames(data_group1) <- group1_gene
write.table(data_group1,"JFCR39_group1.txt",append=T,sep="\t",row.names=T,col.names=T,quote=F)

data_group2_c <- data[group2_gene,group2_samp]
colnames(data_group2) <- group2_samp
rownames(data_group2) <- group2_gene
write.table(data_group2,"JFCR39_group2.txt",append=T,sep="\t",row.names=T,col.names=T,quote=F)

data_group3_c <- data[group3_gene,group3_samp]
colnames(data_group3) <- group3_samp
rownames(data_group3) <- group3_gene
write.table(data_group3,"JFCR39_group3.txt",append=T,sep="\t",row.names=T,col.names=T,quote=F)



#ヒートマップ作成

data <- read.table("jfcr_ratio_for_pc_5_39_50_log.txt",header=F,sep="\t")
nrowdata <- nrow(data)
ncoldata <- ncol(data)
output <- matrix(0,nrowdata,ncoldata)

for(i in 1:ncoldata){
	
	output[,i] <- data[,i]
	
}

#mode(output)でoutputが何型かみれる

colnames(output) <- sname[seq1]
rownames(output) <- as.character(gene2)

# ヒートマップ表示
#install.packages("gplots")
#ターミナルで確認できる→/var/folders/2b/yfm8lh4s4pxdz1jfy_nj7yzm0000gq/T//RtmpMDSlV9/downloaded_packages 
library(gplots)
＃pはピアソンの相関の定義
dist.p <- function(x) as.dist((1-cor(t(x)))/2)
dist.s <- function(x) as.dist((1-cor(t(x), method="spearman"))/2)
#hclust.ward <- function(x) hclust(x, "ward")
hclust.ward <- function(x) hclust(x, "ward.D")
col1 <- redgreen(60)

png(file="jfcr_heatmap_gf_5_39_50_log.png", width=600, height=600)

#Rowv = FALSE  dendrogramにboth,rowを指定した時にTRUEにする必要あり

hv <- heatmap.2(output, Rowv=T, Colv=T, distfun=dist.p, hclustfun=hclust.ward, dendrogram="both", margin=c(9.5,10), breaks=seq(-3,3,0.1), col=col1, scale="none", key=TRUE, keysize=1, symkey=FALSE, density.info="none", trace="none", cexCol=1)


#保存する        
dev.off()     





#0808相関比	
data <- read.table("jfcr_ratio_for_pc_20_39_50_log.txt",sep="\t")
res <- prcomp(t(data))

#######いらないかもしれない
#data <- read.table("JFCR39_mas5_39.txt",header=T,sep="\t")
#sname <- substr(colnames(data), 9, nchar(colnames(data))-21)
#row_samp <- 2:118
#nsamp <- length(row_samp)
#seq1 <- seq(2, nsamp+1, by = 3)

#plot(res$x[,1],res$x[,2],type="n")
#text(res$x[,1],res$x[,2],sname[seq1],cex=0.5)
########

x1_samp <- res$x[,1]
x2_samp <- res$x[,2]

nsamp <- length(res$x[,1])

group1_samp_i <- c()
group2_samp_i <- c()
group3_samp_i <- c()

for(i in 1:nsamp){
	
	if(x1_samp[i] < 0){
		
		group1_samp_i <- c(group1_samp_i,i)
		
	}
	if(x1_samp[i] >= 0 && x2_samp[i] >= 0){
		
		group2_samp_i <- c(group2_samp_i,i)
				
	}
	if(x1_samp[i] >= 0 && x2_samp[i] <= 0){
		
		group3_samp_i <- c(group3_samp_i,i)
				
	}
			
}

#グループ内平均、全体平均を第一主成分第二主成分それぞれでだす

group1_1_mean <- mean(res$x[group1_samp_i,1])
group1_2_mean <- mean(res$x[group1_samp_i,2])
group2_1_mean <- mean(res$x[group2_samp_i,1])
group2_2_mean <- mean(res$x[group2_samp_i,2])
group3_1_mean <- mean(res$x[group3_samp_i,1])
group3_2_mean <- mean(res$x[group3_samp_i,2])
group_all_1_mean <- mean(res$x[,1])
group_all_2_mean <- mean(res$x[,1])

#plot(res$x[,1],res$x[,2],type="n")
#text(res$x[group1_samp_i,1],res$x[group1_samp_i,2],sname[seq1[group1_samp_i]],cex=0.5)

#群内変動
ngroup1 <- length(group1_samp_i)
ngroup2 <- length(group2_samp_i)
ngroup3 <- length(group3_samp_i)

group_all_ch <- 0
group1_ch <- 0
group2_ch <- 0
group3_ch <- 0

for(i in 1:ngroup1){
group1_ch <-group1_ch + (res$x[group1_samp_i[i],1]-group1_1_mean)^2 + (res$x[group1_samp_i[i],2]-group1_2_mean)^2
}
for(i in 1:ngroup2){
group2_ch <-group2_ch + (res$x[group2_samp_i[i],1]-group2_1_mean)^2 + (res$x[group2_samp_i[i],2]-group2_2_mean)^2
}
for(i in 1:ngroup3){
group3_ch <-group3_ch + (res$x[group3_samp_i[i],1]-group3_1_mean)^2 + (res$x[group3_samp_i[i],2]-group3_2_mean)^2
}

group_all_ch <- group1_ch + group2_ch + group3_ch

#群間変動
bet_group_ch <- ((group1_1_mean-group_all_1_mean)^2 + (group1_2_mean-group_all_2_mean)^2)*ngroup1 
+ ((group2_1_mean-group_all_1_mean)^2 + (group2_2_mean-group_all_2_mean)^2)*ngroup2
+ ((group3_1_mean-group_all_1_mean)^2 + (group3_2_mean-group_all_2_mean)^2)*ngroup3

#相関比
rel_rat <- bet_group_ch / (bet_group_ch + group_all_ch)

write.table(rel_rat,"relation_ratio.txt",append=T,sep="\t",row.names=F,col.names=F,quote=F)



#paste関数の使い方

for(i in 1:5){
	
	test <- paste("aaa",i,"bbb")
	filename <- paste(i,".txt",sep="")
	#sepがないと空白が間に入る
	write.table(test,filename,append=T,sep="\t",row.names=F,col.names=F,quote=F)
		
}



#0827相関比 4group

for(j in 1:20){
	
filename <- paste("jfcr_ratio_for_pc_",j,"_39_50_log.txt",sep="") 

data <- read.table(filename,sep="\t")
res <- prcomp(t(data))

#######いらないかもしれない
#data <- read.table("JFCR39_mas5_39.txt",header=T,sep="\t")
#sname <- substr(colnames(data), 9, nchar(colnames(data))-21)
#row_samp <- 2:118
#nsamp <- length(row_samp)
#seq1 <- seq(2, nsamp+1, by = 3)

#plot(res$x[,1],res$x[,2],type="n")
#text(res$x[,1],res$x[,2],sname[seq1],cex=0.5)
########

x1_samp <- res$x[,1]
x2_samp <- res$x[,2]

nsamp <- length(res$x[,1])

group1_samp_i <- c()
group2_samp_i <- c()
group3_samp_i <- c()
group4_samp_i <- c()

for(i in 1:nsamp){
	
	if(x1_samp[i] < 0 && x2_samp[i] >= 0){
		
		group1_samp_i <- c(group1_samp_i,i)
		
	}
		if(x1_samp[i] < 0 && x2_samp[i] <= 0){
		
		group2_samp_i <- c(group1_samp_i,i)
		
	}
		if(x1_samp[i] >= 0 && x2_samp[i] >= 0){
		
		group3_samp_i <- c(group3_samp_i,i)
				
	}
	if(x1_samp[i] >= 0 && x2_samp[i] <= 0){
		
		group4_samp_i <- c(group4_samp_i,i)
				
	}
			
}

#グループ内平均、全体平均を第一主成分第二主成分それぞれでだす

group1_1_mean <- mean(res$x[group1_samp_i,1])
group1_2_mean <- mean(res$x[group1_samp_i,2])
group2_1_mean <- mean(res$x[group2_samp_i,1])
group2_2_mean <- mean(res$x[group2_samp_i,2])
group3_1_mean <- mean(res$x[group3_samp_i,1])
group3_2_mean <- mean(res$x[group3_samp_i,2])
group4_1_mean <- mean(res$x[group4_samp_i,1])
group4_2_mean <- mean(res$x[group4_samp_i,2])
group_all_1_mean <- mean(res$x[,1])
group_all_2_mean <- mean(res$x[,1])

#plot(res$x[,1],res$x[,2],type="n")
#text(res$x[group1_samp_i,1],res$x[group1_samp_i,2],sname[seq1[group1_samp_i]],cex=0.5)

#群内変動
ngroup1 <- length(group1_samp_i)
ngroup2 <- length(group2_samp_i)
ngroup3 <- length(group3_samp_i)
ngroup4 <- length(group4_samp_i)

group_all_ch <- 0
group1_ch <- 0
group2_ch <- 0
group3_ch <- 0
group4_ch <- 0

for(i in 1:ngroup1){
group1_ch <-group1_ch + (res$x[group1_samp_i[i],1]-group1_1_mean)^2 + (res$x[group1_samp_i[i],2]-group1_2_mean)^2
}
for(i in 1:ngroup2){
group2_ch <-group2_ch + (res$x[group2_samp_i[i],1]-group2_1_mean)^2 + (res$x[group2_samp_i[i],2]-group2_2_mean)^2
}
for(i in 1:ngroup3){
group3_ch <-group3_ch + (res$x[group3_samp_i[i],1]-group3_1_mean)^2 + (res$x[group3_samp_i[i],2]-group3_2_mean)^2
}
for(i in 1:ngroup4){
group4_ch <-group4_ch + (res$x[group4_samp_i[i],1]-group4_1_mean)^2 + (res$x[group4_samp_i[i],2]-group4_2_mean)^2
}

group_all_ch <- group1_ch + group2_ch + group3_ch + group4_ch

#群間変動
bet_group_ch <- ((group1_1_mean-group_all_1_mean)^2 + (group1_2_mean-group_all_2_mean)^2)*ngroup1 
+ ((group2_1_mean-group_all_1_mean)^2 + (group2_2_mean-group_all_2_mean)^2)*ngroup2
+ ((group3_1_mean-group_all_1_mean)^2 + (group3_2_mean-group_all_2_mean)^2)*ngroup3
+ ((group4_1_mean-group_all_1_mean)^2 + (group4_2_mean-group_all_2_mean)^2)*ngroup4
#相関比
rel_rat <- bet_group_ch / (bet_group_ch + group_all_ch)

write.table(rel_rat,"relation_ratio_4.txt",append=T,sep="\t",row.names=F,col.names=F,quote=F)

}





#0902主成分がん種類　5/39,7/39,11/39

num <- c(5,7,11)
for(z in num){
	
filename <- paste("jfcr_ratio_for_pc_",z,"_39_50_log.txt",sep="")

data <- read.table(filename,sep="\t")
res <- prcomp(t(data))

#res$xはサンプル39個の第○主成分（新軸）の主成分得点が書いてある
#res$rotationは第○主成分（新軸）に対してそれぞれの遺伝子がどんな係数になってるか書いてある
	
#サンプル名格納
#リストの使い方
#ベクトルでは要素一つ一つが同じ方で同じ重みでなくてはいけないので今回はリストを用いた
#リストの中の要素を指定する場合は二重カッコになる点に注意
data_sname <- read.table("sname_english.txt",sep="\t")

sname <- as.character(data_sname[,1])
type <- as.character(data_sname[,2])
nsamp <- length(type)

type_snum <- list(c(),c(),c(),c(),c(),c(),c(),c())
cancer <- c("lung","prostate","breast","colon","stomach","ovary","brain","other")
ncancer <- length(cancer)

for(i in 1:nsamp){
	
	for(j in 1:ncancer){
		if(as.character(type[i]) == cancer[j]){		
			type_snum[[j]] <- c(type_snum[[j]],i)
		}
	}
	
}

outputname <- paste("jfcr39_cancer_type_plot_",z,".pdf",sep="")
pdf(outputname)
plot(res$x[,1],res$x[,2],type="n")
for(i in 1:6){	
	
	text(res$x[type_snum[[i]],1],res$x[type_snum[[i]],2],sname[type_snum[[i]]],cex=0.8,col=i)		
}
#col8,9がみにくかったから
text(res$x[type_snum[[7]],1],res$x[type_snum[[7]],2],sname[type_snum[[7]]],cex=0.8,col="darkorange")
text(res$x[type_snum[[8]],1],res$x[type_snum[[8]],2],sname[type_snum[[8]]],cex=0.8,col="darkviolet")

dev.off()
}



#141003
#プログレス用ヒートマップ作成

data <- read.table("jfcr_ratio_for_pc_5_39_50_log.txt",header=F,sep="\t")
nrowdata <- nrow(data)
ncoldata <- ncol(data)
output <- matrix(0,nrowdata,ncoldata)
gene_list <- as.character(data[,1])

for(i in 1:ncoldata){
	
	output[,i] <- data[,i]
	
}

data_samp <- read.table("JFCR39_mas5_39.txt",header=T,sep="\t",row.names=1)
seq_gf <- seq(1, nsamp, by = 3)
sname <- substr(colnames(data_samp), 9, nchar(colnames(data_samp))-21)
colname <- sname[seq_gf]
data_gene <- read.table("0414_JFCR39_mas5_gene.txt",header=T,sep="\t",quote="\"",row.names=1)	
gene <- as.character(data_gene[gene_list,1])
#///で区切られて遺伝子が2つ以上記載されているプローブがあってプロットするとみにくいから消す
gene2 <- gsub("///","",gene)

colnames(output) <- colname
rownames(output) <- as.character(gene2)

# ヒートマップ表示
#install.packages("gplots")
#ターミナルで確認できる→/var/folders/2b/yfm8lh4s4pxdz1jfy_nj7yzm0000gq/T//RtmpMDSlV9/downloaded_packages 
library(gplots)
＃pはピアソンの相関の定義
dist.p <- function(x) as.dist((1-cor(t(x)))/2)
dist.s <- function(x) as.dist((1-cor(t(x), method="spearman"))/2)
#hclust.ward <- function(x) hclust(x, "ward")
hclust.ward <- function(x) hclust(x, "ward.D")
col1 <- redgreen(60)

png(file="jfcr_heatmap_gf_5_39_50_log.png", width=600, height=600)

#Rowv = FALSE  dendrogramにboth,rowを指定した時にTRUEにする必要あり

hv <- heatmap.2(output, Rowv=T, Colv=T, distfun=dist.p, hclustfun=hclust.ward, dendrogram="both", margin=c(9.5,10), breaks=seq(-3,3,0.1), col=col1, scale="none", key=TRUE, keysize=1, symkey=FALSE, density.info="none", trace="none", cexCol=1)


#保存する        
dev.off()  


#hypo

#hypoのレシオ、5サンプル以上で発現2倍以上/1/2以下のプローブのみを出力する
data <- read.table("JFCR39_mas5_39.txt",header=T,sep="\t")

ngene <- nrow(data)
nsamp <- ncol(data)-1
nsamp3 <- nsamp/3

seq1 <- seq(2,nsamp+1,by=3)
seq2 <- seq(3,nsamp+1,by=3)
seq3 <- seq(4,nsamp+1,by=3)

up_count <- 0
down_count <- 0
th.rat <- 2
th.count <- 20

output <- c()

for(i in 1:ngene){
	
	for(j in 1:nsamp3){
		
		rat <- pmax(data[i,seq1[j]],50)/pmax(data[i,seq3[j]],50)
			
			if(rat >= th.rat){		
				up_count = up_count + 1						
			}
			if(rat <= 1/th.rat){				
				down_count = down_count + 1						
			}
	}	

	if(up_count >= th.count | down_count >= th.count){
		
		output <- c(output,as.character(data[i,1]))
		
	}
	
	up_count <- 0
	down_count <- 0
			
}

write.table(output,"jfcr_gene_for_pc_20.txt",append=T,sep="\t",row.names=F,col.names=F,quote=F)

#出力したプローブリストより、その行のレシオを出力
data <- read.table("JFCR39_mas5.txt",header=T,sep="\t",row.names=1)
gene <- read.table("jfcr_gene_for_pc_20.txt",header=F,sep="\t")

gene_list <- as.character(gene[,1])

nsamp <- ncol(data)
nsamp3 <- nsamp/3
ngene <- length(gene_list)

seq1 <- seq(1,nsamp,by=3)
seq2 <- seq(2,nsamp,by=3)
seq3 <- seq(3,nsamp,by=3)

output <- matrix(0,ngene,nsamp/3)

for(i in 1:nsamp3){

rat <- log(data[gene_list,seq1[i]]/data[gene_list,seq3[i]])
output[,i] <- rat

}

write.table(output,"jfcr_ratio_for_pc_20_log.txt",append=T,sep="\t",row.names=F,col.names=F,quote=F)









#141014_主成分分析_係数が大きい上位10個のプローブを出力
	
data <- read.table("jfcr_ratio_for_pc_5_39_50_log.txt",sep="\t")
res <- prcomp(t(data))
	
#サンプル名格納
data <- read.table("JFCR39_mas5_39.txt",header=T,sep="\t")
	
nsamp <- 117
	
seq1 <- seq(2,nsamp+1,by=3)
seq2 <- seq(3,nsamp+1,by=3)
seq3 <- seq(4,nsamp+1,by=3)
	
colname_gf <- colnames(data[seq1])
colname_n <- colnames(data[seq3])
	
sname <- substr(colnames(data), 9, nchar(colnames(data))-21)
row_samp <- 2:118
nsamp <- length(row_samp)
seq1 <- seq(2, nsamp+1, by = 3)
	
plot(res$x[,1],res$x[,2],type="n")
text(res$x[,1],res$x[,2],sname[seq1],cex=0.5)
	
#主成分分析の結果として標準偏差を返す
summary(res)

#プロットするときのために遺伝子名を格納しておく
data <- read.table("0414_JFCR39_mas5_gene.txt",header=T,sep="\t",quote="\"",row.names=1)
gene <- read.table("jfcr_gene_for_pc_5_39_50.txt",header=F,sep="\t")	
gene_list <- as.character(gene[,1])
gene <- as.character(data[gene_list,1])

#///で区切られて遺伝子が2つ以上記載されているプローブがあってプロットするとみにくいから消す
gene2 <- gsub("///","",gene)

plot(res$rotation[,1],res$rotation[,2],type="n")
text(res$rotation[,1],res$rotation[,2],gene2,cex=0.5)

#ここから新しい

pc1_gene <- matrix(0,length(gene2),2)
pc1_gene[,1] <- abs(res$rotation[,1])
pc1_gene[,2] <- gene2
pc2_gene <- matrix(0,length(gene2),2)
pc2_gene[,1] <- abs(res$rotation[,2])
pc2_gene[,2] <- gene2

#sort 
od <- order(pc1_gene[,1],decreasing=T)
pc1_gene <- pc1_gene[od,]
od <- order(pc2_gene[,1],decreasing=T)
pc2_gene <- pc2_gene[od,]

write.table(pc1_gene,"jfcr_pc1_gene_5_39_50_log.txt",append=T,sep="\t",row.names=F,col.names=F,quote=F)
write.table(pc2_gene,"jfcr_pc2_gene_5_39_50_log.txt",append=T,sep="\t",row.names=F,col.names=F,quote=F)



