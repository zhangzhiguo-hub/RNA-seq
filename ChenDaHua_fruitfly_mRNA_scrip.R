###ballgown包处理数据,输出bg_filt_trans
{
  #加载ballgown包
  library(ballgown)
  #查看当前路径
  getwd()
  #设置存储路径
  setwd("E:/scientific_research/Transcriptome_analysis/fruitfly_maternal")
  #重新查看更改过后的当前路径
  getwd()
  #读取定量分析的表型数据（批量操作型）
  bg = ballgown(dataDir = "./ballgown", samplePattern = "C", meas = "all")
  #过滤掉表达量低的基因,即每一行所代表的基因，其后的各值都为0的
  bg_filt=as.data.frame(gexpr(bg))
  bg_filt_row=as.data.frame(bg@expr$trans)
  bg_filt_trans = subset(bg_filt,apply(bg_filt,1,max)>20,genomesubset=TRUE)
  colnames(bg_filt_trans)
  
  #feature:在那个水平计算，可选“gene”,"transcript","exon","intron"
  #协变量covariate是对因变量有影响的变量，它不是研究者研究的自变量，不为实验者操作，但对试验有影响。
  #主变量adjustvars
  #getFC可以指定输出结果显示组间表达量的foldchange(差异倍数，FC=Exp x1/x2)
  #meas:measurement采用哪种计算方式
  #mod0=model.matrix(~?, pData(bg_filt)$treat)
  #mod=model.matrix(~?, pData(bg_filt)$treat)
  #results_tanscripts=stattest(bg_filt,feature = "transcript",covariat = ?,adjustvars = ?,mod0 = mod0,mod = mod,meas = "FPKM")本例中无重复，此法不可用来做差异化分析。
}
#读入列表bg_filt_row，重命名并重排序，读出原始列表bg_filt_row_1
{
  library(dplyr)
  bg_filt_row_1= select(bg_filt_row,c(t_name,gene_name,FPKM.C29_L4_340340,FPKM.C29_L4_340340,FPKM.C30_L4_341341,FPKM.C31_L4_342342,FPKM.C32_L4_343343,FPKM.C33_L4_345345,FPKM.C34_L4_346346,FPKM.C35_L4_347347,FPKM.C36_L1_349349,FPKM.C37_L4_350350,FPKM.C38_L4_351351,FPKM.C39_L4_352352,FPKM.C40_L4_353353))
  #重命名
  {
    library(reshape)
    #使用该函数时变量名经常出问题。使用dplyr::rename重命名了变量，以便名称不再取决于变量，而取决于字符。后面如果再出现问题，跟换变量前后的位置即可
    bg_filt_row_1 = dplyr::rename(bg_filt_row_1,c("WT0.5-1h_1" = FPKM.C29_L4_340340,
                                                  "WT0.5-1h_2" = FPKM.C30_L4_341341,
                                                  "MT0.5-1h_1" = FPKM.C31_L4_342342,
                                                  "MT0.5-1h_2" = FPKM.C32_L4_343343,
                                                  "WT2-3h_1" = FPKM.C33_L4_345345,
                                                  "WT2-3h_2" = FPKM.C34_L4_346346,
                                                  "MT2-3h_1" = FPKM.C35_L4_347347,
                                                  "MT2-3h_2" = FPKM.C36_L1_349349,
                                                  "WT5.5-6h_1" = FPKM.C37_L4_350350,
                                                  "WT5.5-6h_2" = FPKM.C38_L4_351351,
                                                  "MT5.5-6h_1" = FPKM.C39_L4_352352,
                                                  "MT5.5-6h_2" = FPKM.C40_L4_353353))
    head(bg_filt_row_1,1)#查看首行
  }
  #给列重排序
  {
    #dplyr包，该宏包功能强大。是数据科学tidyverse集合的核心部件之一。学习链接https://bookdown.org/wangminjie/R4DS/colwise.html
    library(dplyr,warn.conflicts = F)
    bg_filt_row_1 = bg_filt_row_1 %>% select("t_name","gene_name","WT0.5-1h_1","WT0.5-1h_2","WT2-3h_1","WT2-3h_2","WT5.5-6h_1","WT5.5-6h_2","MT0.5-1h_1","MT0.5-1h_2","MT2-3h_1","MT2-3h_2","MT5.5-6h_1","MT5.5-6h_2")
  }
}
####读出并保存数据
{
  write.table(bg_filt_trans,file = "fruitfly_maternal_bg_filt_trans.csv",sep = ",",row.names = T) #列表保存为.csv格式
  write.table(bg_filt_trans,file = "fruitfly_maternal_bg_filt_trans.txt",sep = " ",row.names = T) #列表保存为.txt格式
}
getwd()
#设置存储路径
setwd("E:/scientific_research/Transcriptome_analysis/fruitfly_maternal")
#重新查看更改过后的当前路径
getwd()
#读入数据为列表格式，以第一列为行名称
{
  data = read.table("fruitfly_maternal_bg_filt_trans.txt",
                    header = T,
                    row.names = 1)
}
#去掉列表中中表达量异常的基因
{
  data[data==0] = NA
  data = na.omit(data)
}
#给列表重命名并重排序
{
  #重命名
  {
    library(reshape)
    #使用该函数时变量名经常出问题。使用dplyr::rename重命名了变量，以便名称不再取决于变量，而取决于字符。后面如果再出现问题，跟换变量前后的位置即可
    data = dplyr::rename(data,c("WT0.5-1h_1" = FPKM.C29_L4_340340,
                                "WT0.5-1h_2" = FPKM.C30_L4_341341,
                                "MT0.5-1h_1" = FPKM.C31_L4_342342,
                                "MT0.5-1h_2" = FPKM.C32_L4_343343,
                                "WT2-3h_1" = FPKM.C33_L4_345345,
                                "WT2-3h_2" = FPKM.C34_L4_346346,
                                "MT2-3h_1" = FPKM.C35_L4_347347,
                                "MT2-3h_2" = FPKM.C36_L1_349349,
                                "WT5.5-6h_1" = FPKM.C37_L4_350350,
                                "WT5.5-6h_2" = FPKM.C38_L4_351351,
                                "MT5.5-6h_1" = FPKM.C39_L4_352352,
                                "MT5.5-6h_2" = FPKM.C40_L4_353353))
    head(data,1)#查看首行
  }
  #归类
  {
    WT = c("WT0.5-1h_1","WT0.5-1h_2","WT2-3h_1","WT2-3h_2","WT5.5-6h_1","WT5.5-6h_2")
    MT = c("MT0.5-1h_1","MT0.5-1h_2","MT2-3h_1","MT2-3h_2","MT5.5-6h_1","MT5.5-6h_2")
  }
  #给列重排序
  {
    #dplyr包，该宏包功能强大。是数据科学tidyverse集合的核心部件之一。学习链接https://bookdown.org/wangminjie/R4DS/colwise.html
    library(dplyr,warn.conflicts = F)
    data = data %>% select("WT0.5-1h_1","WT0.5-1h_2","WT2-3h_1","WT2-3h_2","WT5.5-6h_1","WT5.5-6h_2","MT0.5-1h_1","MT0.5-1h_2","MT2-3h_1","MT2-3h_2","MT5.5-6h_1","MT5.5-6h_2")
  }
}
#取均值，输出newdata
#成功#！要比对重复的时勿执行
{
  i = 1
  newdata = data
  while (i < length(colnames(data))){
    newdata$A <- rowMeans(data[,c(i,i+1)])
    a = c(substr(colnames(data[i]),1,nchar(colnames(data[i]))-2))
    colnames(newdata)[length(newdata)] = a
    i = i+2
  }
  newdata = newdata[,-which(colnames(newdata) %in% colnames(data))]
}
#WT原始聚类,输出data_wt_cluster，并取log2，输入为data(newdata)
data = newdata#根据均值计算
WT = colnames(newdata)[c(1,2,3)]
{
  data_wt = data[,which(colnames(data) %in% WT)]
  #data_wt_cluster
  i = 1
  data_wt_cluster = data.frame()
  while(i <= length(data_wt))
  {
    a =  subset(data_wt,data_wt[,i] > apply(data_wt[,colnames(data_wt)[-i]],1,max))
    #给行按照colnames(data_wt)[i]的值升序排序
    a = a[order(-a[i]),]
    data_wt_cluster = rbind(data_wt_cluster,a)
    i <- i + 1
  }
  #求取log2值
  data_wt_cluster = log2(data_wt_cluster[,1:base::ncol(data_wt_cluster)])
}
#WT和MT共聚类，筛选输出data_cluster，并取log2
data = newdata#根据均值计算
{
  data_cluster = data[rownames(data_wt_cluster),]
  #取log2
  data_cluster = log2(data_cluster)
}

#筛选数据中 突变型“最大值列与下一列的差”小于 野生型“最大值列与下一列的差”的list
#Maternal_mRNA :WT(0.5~1h)-WT(5.5~6h) > MUT(0.5~1h)-MUT(5.5~6h)##输出data_cluster1,结果取log2
{
  data_cluster1 = subset(data_wt,data_wt[,1]>20)
  data_cluster1 = subset(data_wt,data_wt[,1] > data_wt[,2] & data_wt[,2] >= data_wt[,3])
  data_cluster1 = newdata[rownames(data_cluster1),]
  data_cluster1 = subset(data_cluster1,(data_cluster1[,1] - data_cluster1[,3]) > (data_cluster1[,4] - data_cluster1[,6]),genomesubset=TRUE)
  data_cluster1 = log2(data_cluster1)
  a = tibble::as_tibble(bg_filt_row_1[which(bg_filt_row_1$gene_name %in% rownames(data_cluster1)),])
  #去掉列表中中表达量异常的基因
  {
    a[a==0] = NA
    a = na.omit(a)
  }
  a$gene_name = factor(a$gene_name)
  a = a %>% group_by(gene_name) %>% top_n(1,`WT0.5-1h_1`)
  write.csv(a,file = paste(getwd(),c("Maternal_RNA.csv"),sep = "/"))
}
#Zygote_mRNA ##未知##输出data_cluster3
{
  data_cluster3 = subset(data_wt,apply(data_wt[,colnames(data_wt)[-1]],1,max)>20)
  data_cluster3 = subset(data_cluster3,apply(data_cluster3[,colnames(data_cluster3)[-1]],1,min)/data_cluster3[,1]>1.5)
  data_cluster3 = newdata[rownames(data_cluster3),]
  a = tibble::as_tibble(bg_filt_row_1[which(bg_filt_row_1$gene_name %in% rownames(data_cluster3)),])
  #去掉列表中中表达量异常的基因
  {
    a[a==0] = NA
    a = na.omit(a)
  }
  a$gene_name = factor(a$gene_name)
  a = a %>% group_by(gene_name) %>% top_n(1,`WT5.5-6h_1`)
  write.csv(a,file = paste(getwd(),c("Zygote.csv"),sep = "/"))
}

#根据筛选的基因，返回上游的bedgraph文件，查看不同转录本在3'5'端的分布
{
  filename = list.files("C:/Users/Administrator/Desktop/ChengDaHua/fly_mrna/")
  i = 1
  a = data.frame()
  rm(sum)
  for (i in filename) 
  {
    A = read.table(file = paste("C:/Users/Administrator/Desktop/ChengDaHua/fly_mrna/",i,sep = ""))
    reference = read.table(file = "C:/Users/Administrator/Desktop/ChengDaHua/FF.txt")
    reference = reference[which(reference$V1 %in% unique(unlist(A$V1))),]
    A$V1 = factor(A$V1)
    #去掉列表中中表达量为0的转录本
    {
      A[A==0] = NA
      A = na.omit(A)
    }
    A = data.frame(tapply(A$V3,A$V1, median)) %>% (function(x){cbind(x,V1 = rownames(x))}) %>% (function(x){merge(A,x)})
    colnames(A)[length(A)]="Median"
    A$V3 = (A$V3/A$Median)
    A = left_join(A,reference,by = "V1",copy = T)
    #分bin
    A$site = (A$V2.x/A$V2.y)*100
    A$bin = factor(ceiling(A$site))
    #counts数除以转录本长度,为了在bin中进行标准化
    names(A)[3]=c("The number of reads")
    A$`The number of reads`=A$`The number of reads`/A$V2.y
    A = data.frame(aggregate(A$`The number of reads`~A$bin, A,sum)) %>% (function(x){cbind(x,site=rownames(x))}) 
    A$group = substr(i,1,nchar(i)-3)
    A$Type1=ifelse(unlist(strsplit(i,"_"))[2]=="Zygote","Zygote","Maternal")
    A$Type2=ifelse(substr(i,1,2)=="WT","WT","MT")
    a = rbind(a,A)
  }
  #画图
  library(ggplot2)
    A$A.bin=factor(A$A.bin)
    title="mRNA"
    ggplot(data=a,mapping=aes(x=A.bin,y=A..The.number.of.reads.,group=group,colour = group,linetype=Type1,shape=Type2))+
    geom_line(size = 0.5)+
    geom_point(size = 1.24)+
    labs(x = "Bins",y = "The number of reads",title = "mRNA")+
    theme(panel.grid = element_blank())+#去掉网格线
    scale_x_discrete(breaks=seq(0,100,10))+#修改x轴的刻度间隔,使用discrete时，行必须是因子
    scale_y_discrete(breaks=seq(0,15,3))
}   

#新大型聚类,输出cluster_data,需要运行720*6678次。。。。失败
data = newdata
{
  a = c(1:6)
  #生成1:6的所有随机组合并由小到大排列,num向量
  #order()函数的排序原则，（最小，次小，。。最大）
  {
    num = c() 
    b=c()
    a = (1:6)
    for (i in a) {
      b[1]=a[i]
      for (j in a[-i]) {
        b[2]=a[j]
        for (k in a[c(-i,-j)]) {
          b[3]=a[k]
          for (l in a[c(-i,-j,-k)]) {
            b[4]=a[l]
            for (m in a[c(-i,-j,-k,-l)]) {
              b[5]=a[m]
              for (n in a[c(-i,-j,-k,-l,-m)]) {
                b[6]=a[n]
                num = c(num,paste(b,collapse = ""))
              }
            }
          }
        }
      }
    }
    num = num[c(grep(".....1",num),grep(".....2",num),grep(".....3",num),grep(".....4",num),grep(".....5",num),grep(".....6",num))]
  }
  cluster_data = data.frame()
  i=1
  library(stringr)
  while(i<=length(num)) {
    a = data.frame()
    j=1
    while(j<=nrow(data)) {
      if (str_c(as.character(order((data[j,]))),collapse = "") == num[i]) {
        a = rbind(a,data[j,])
      }
      j = j + 1
    }
    if (nrow(a)!=0) {
      a = a[order(-a[which(a[1,] >= max(a[1,]))]),]
    }
    cluster_data = rbind(cluster_data,a)
    i = i + 1
  }
  #取log2
  cluster_data = log2(cluster_data[,1:base::ncol(cluster_data)])
}

#未完成
#GO分析———在线，先做总表的GO分析，此处为MAXrow为1以上的,文件名为RNA in nemateode mutants.csv#DAVID无法识别超过3000的ID量
{
  write.table(data.frame(geneid = row.names(data)),file = "geneid1.txt",row.names = F,col.names = F,quote = F)#将基因名作为表格输出，并保存为.txt格式，不要行名，不要列名,不要引号
  #-----
  #本数据基因id来源于Wormbase数据库,格式为“例：WBGene00009237”
  #软件名：DAVID  https://david.ncifcrf.gov/,复制粘贴提交，Identifier选择WORMBASE_GENE_ID.，这里选用其中的 Gene ID Conversion Tool，转到OFFICAL_GENE_SYMBOL（国际通用名称）下载即可，命名为all.txt
  #读入并重命名
  {
    GO_gene_name = read.table("all.txt",header = T,sep ="\t",#Tap分隔
                              fill=TRUE,#该参数可忽略数据不规则的错误，直接读取
                              stringsAsFactors = FALSE,quote = ""#不要引号
    )
    cluster = log2(data)
    cluster$Wormbase_id = c(Wormbase_id = row.names(cluster))
    cluster$MAX_MUT = apply(cluster[2:4],1,max)
    cluster$MIN_MUT = apply(cluster[2:4],1,min)
    library(reshape)
    GO_gene_name = dplyr::rename(GO_gene_name,c("Wormbase_id" = From,
                                                "Gene Name" = To,
                                                "Discription" = Gene.Name
    ))
  }
  #删除缺失值
  {
    GO_gene_name[GO_gene_name==""] = NA
    GO_gene_name = na.omit(GO_gene_name)
  }
  #合并
  {
    cluster_go = merge(cluster,GO_gene_name,sort = F,all = T)
    rownames(cluster_go) = cluster_go[,"Wormbase_id"]
    #排序
    #列排序
    library(dplyr,warn.conflicts = F)
    colnames(cluster_go)
    cluster_go = cluster_go %>% select("Wormbase_id","Gene Name","Species","Discription","CTRL","MUT_qx666","MUT_epg.5","MUT_dpy.7","MAX_MUT","MIN_MUT")
    
  }
  #读出列表及注释
  {
    write.csv(cluster_go,file = "RNA in nemateode mutants.csv")
  }
  #-----
}

#分开GO
{
  #data_cluster1,GO分析———在线，hclust聚类，输出列表为：cluster1_go,输出文件名为WT(0.5~1h)_MAX of fruitfly.csv
  {
    write.table(data.frame(geneid = row.names(data_cluster1)),file = "data_cluster1.txt",row.names = F,col.names = F,quote = F)#将基因名作为表格输出，并保存为.txt格式，不要行名，不要列名,不要引号
    #-----
    #本数据基因id为sympol
    #软件名：DAVID  https://david.ncifcrf.gov/,复制粘贴提交，Identifier选择WORMBASE_GENE_ID.，这里选用其中的 Gene ID Conversion Tool，转到OFFICAL_GENE_SYMBOL（国际通用名称）下载即可，命名为all.txt
    #读入并重命名
    {
      GO_gene_name = read.table("data_cluster1.txt",header = T,sep ="\t",#Tap分隔
                                fill=TRUE,#该参数可忽略数据不规则的错误，直接读取
                                stringsAsFactors = FALSE,quote = ""#不要引号
      )
      GO_gene_name = GO_gene_name[which(GO_gene_name[,"Species"]=="Drosophila melanogaster"),]
      cluster1 = data_cluster1#log2已经取过
      #聚类
      dist_method = "euclidean"
      #计算类间距离的方法：类平均法：average;重心法：centroid;中间距离法：median;最长距离法：complete(默认);最短距离法：single;离差平方和法：ward;密度估计法：density
      cluster_method = "complete"
      b = dist(cluster1,method = dist_method)
      b = hclust(b,method = cluster_method)
      order = b$order
      cluster1 = cluster1[order,]
      cluster1$SYMBOL_ID = c(SYMBOL_ID  = row.names(cluster1))
      #改列名
      library(reshape)
      GO_gene_name = dplyr::rename(GO_gene_name,c("SYMBOL_ID" = From,
                                                  "NCBI_ID" = To,
                                                  "Discription" = Gene.Name
      ))
    }
    #删除缺失值
    {
      GO_gene_name[GO_gene_name==""] = NA
      GO_gene_name = na.omit(GO_gene_name)
    }
    #合并
    {
      cluster1_go = merge(cluster1,GO_gene_name,sort = F,all = T)
      #rownames(cluster1_go) = cluster1_go[,"SYMBOL_ID"]
      #排序
      #列排序
      library(dplyr,warn.conflicts = F)
      colnames(cluster1_go)
      cluster1_go = cluster1_go %>% select("SYMBOL_ID","NCBI_ID","Species","Discription","WT0.5-1h","WT2-3h","WT5.5-6h","MT0.5-1h","MT2-3h","MT5.5-6h")
    }
    #读出列表及注释
    {
      write.csv(cluster1_go,file = "WT(0.5~1h)_MAX of fruitfly.csv")
    }
    #-----
  }
  #data_cluster2,GO分析———在线，hclust聚类，输出列表为：cluster2_go,输出文件名为WT(2~3h)_MAX of fruitfly.csv
  {
    write.table(data.frame(geneid = row.names(data_cluster2)),file = "data_cluster2.txt",row.names = F,col.names = F,quote = F)#将基因名作为表格输出，并保存为.txt格式，不要行名，不要列名,不要引号
    #-----
    #本数据基因id为sympol
    #软件名：DAVID  https://david.ncifcrf.gov/,复制粘贴提交，Identifier选择WORMBASE_GENE_ID.，这里选用其中的 Gene ID Conversion Tool，转到OFFICAL_GENE_SYMBOL（国际通用名称）下载即可，命名为all.txt
    #读入并重命名
    {
      GO_gene_name = read.table("data_cluster2.txt",header = T,sep ="\t",#Tap分隔
                                fill=TRUE,#该参数可忽略数据不规则的错误，直接读取
                                stringsAsFactors = FALSE,quote = ""#不要引号
      )
      GO_gene_name = GO_gene_name[which(GO_gene_name[,"Species"]=="Drosophila melanogaster"),]
      cluster2 = data_cluster2#log2已经取过
      #聚类
      dist_method = "euclidean"
      #计算类间距离的方法：类平均法：average;重心法：centroid;中间距离法：median;最长距离法：complete(默认);最短距离法：single;离差平方和法：ward;密度估计法：density
      cluster_method = "complete"
      b = dist(cluster2,method = dist_method)
      b = hclust(b,method = cluster_method)
      order = b$order
      cluster2 = cluster2[order,]
      cluster2$SYMBOL_ID = c(SYMBOL_ID  = row.names(cluster2))
      #改列名
      library(reshape)
      GO_gene_name = dplyr::rename(GO_gene_name,c("SYMBOL_ID" = From,
                                                  "NCBI_ID" = To,
                                                  "Discription" = Gene.Name
      ))
    }
    #删除缺失值
    {
      GO_gene_name[GO_gene_name==""] = NA
      GO_gene_name = na.omit(GO_gene_name)
    }
    #合并
    {
      cluster2_go = merge(cluster2,GO_gene_name,sort = F,all = T)
      #rownames(cluster1_go) = cluster1_go[,"SYMBOL_ID"]
      #排序
      #列排序
      library(dplyr,warn.conflicts = F)
      colnames(cluster2_go)
      cluster2_go = cluster2_go %>% select("SYMBOL_ID","NCBI_ID","Species","Discription","WT0.5-1h","WT2-3h","WT5.5-6h","MT0.5-1h","MT2-3h","MT5.5-6h")
    }
    #读出列表及注释
    {
      write.csv(cluster2_go,file = "WT(2~3h)_MAX of fruitfly.csv")
    }
    #-----
  }
}


#绘大热图
{
  library(pheatmap)
  gaps_row = c()
  sum = 0
  Cluster = c()
  length_Cluster = c()
  i = 1
  while (i <= length(colnames(data_wt_cluster))) #绘制不同的图记得改名儿
  {
    a = nrow(subset(data_wt_cluster,data_wt_cluster[,i] > apply(data_wt_cluster[,colnames(data_wt_cluster)[-i]],1,max)))
    b = paste ("Max",
               colnames(data_wt_cluster)[i],
               sep = "_")
    b = paste (b,
               a,
               sep = " ")
    b = paste (b,
               round(a/nrow(data_wt_cluster)*100),
               sep = "(")
    Cluster[i] = paste (b,
                        c("%)"),
                        sep = ""
    )
    length_Cluster[i] = a
    sum = sum+a
    gaps_row[i] = sum
    i = i + 1
  }
  #画data_cluster聚类的图
  {
    #添加注释
    annotation_col= data.frame(Time = c(0.75,2.5,5.75,0.75,2.5,5.75))
    annotation_row = data.frame(Cluster = rep(Cluster,length_Cluster))
    color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")#调色
    #运行没反应时，运行dev.off()
    {
      rownames(annotation_col)=colnames(data_cluster)
      rownames(annotation_row)=rownames(data_cluster)
      #data_cluster
      {
        pheatmap(data_cluster,
                 color = colorRampPalette(color.key)(50),
                 border_color = NA,
                 cluster_cols = F,
                 cluster_rows = F,
                 #分隔
                 gaps_row = gaps_row,
                 annotation_col = annotation_col,
                 annotation_row = annotation_row,
                 main = "RNA of WT and MT",
                 #scale = "row",
                 show_rownames = F,show_colnames = T,
                 cellwidth = 50, cellheight = 0.1,
                 fontsize = 20
        )
      }
    }
  }
  
}  
#绘制小热图
data_cluster=data_cluster1
main = "WT(0.5~1h)_MAX of fruitfly:
WT(0.5~1h)-WT(2~3h) > MUT(0.5~1h)-MUT(2~3h)"
data_cluster=data_cluster2
main = "WT(2~3h)_MAX of fruitfly:
WT(2~3h)-WT(5.5~6h) > MUT(2~3h)-MUT(5.5~6h)"
{
  color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")#调色
  #运行没反应时，运行dev.off()
  {
    #data_cluster
    {
      library(pheatmap)
      pheatmap(data_cluster,
               color = colorRampPalette(color.key)(50),
               border_color = NA,
               cluster_cols = F,
               cluster_rows = T,
               #分隔
               #gaps_row = gaps_row,
               #annotation_col = annotation_col,
               #annotation_row = annotation_row,
               main = main,
               #scale = "row",
               show_rownames = F,show_colnames = T,
               cellwidth = 50, cellheight = 0.5,
               fontsize = 20
      )
    }
  }
}



##画热图
{
  library(pheatmap)
  #将data转化为matrix格式
  #分别画AB两组的图
  {
    #添加注释
    annotation_col = data.frame(Grop= rep(c("Control","Treat"),c(1,5)),Time = c(0,1:5))
    color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")#调色
    #运行没反应时，运行dev.off()
    
    #画A组的基因热图
    {
      annotation_row_a = data.frame(Cluster = rep(c("Max_CTRL_173(11%)","Max_A15_775(51%)","Max_A30_165(11%)","Max_A60_114(8%)","Max_A90_115(8%)","Max_A120_169(11%)"),c(173,775,165,114,115,169)))
      rownames(annotation_col)=colnames(data_filt_a_cluster)
      rownames(annotation_row_a)=rownames(data_filt_a_cluster)
      #data_filt_a_cluster
      {
        a = pheatmap(data_filt_a_cluster,
                     color = colorRampPalette(color.key)(50),
                     border_color = NA,
                     cluster_cols = F,
                     cluster_rows = F,
                     #分隔
                     gaps_row = 
                       c(
                         base::nrow(data_filt_a_cluster1), 
                         base::nrow(data_filt_a_cluster1) + base::nrow(data_filt_a_cluster2),
                         base::nrow(data_filt_a_cluster1) + base::nrow(data_filt_a_cluster2) + base::nrow(data_filt_a_cluster3),
                         base::nrow(data_filt_a_cluster1) + base::nrow(data_filt_a_cluster2) + base::nrow(data_filt_a_cluster3) + base::nrow(data_filt_a_cluster4),
                         base::nrow(data_filt_a_cluster1) + base::nrow(data_filt_a_cluster2) + base::nrow(data_filt_a_cluster3) + base::nrow(data_filt_a_cluster4) + base::nrow(data_filt_a_cluster5)
                       ),
                     annotation_col = annotation_col,
                     annotation_row = annotation_row_a,
                     main = "Cluster about A",
                     #scale = "row",
                     show_rownames = F,show_colnames = T,
                     cellwidth = 100, cellheight = 0.5
        )
      }
    }
    
    #画B组的基因热图
    {
      annotation_row_b = data.frame(Cluster = rep(c("Max_CTRL_106(16%)","Max_B15_98(15%)","Max_B30_107(17%)","Max_B60_96(15%)","Max_B90_102(16%)","Max_B120_132(21%)"),c(106,98,107,96,102,132)))
      rownames(annotation_col)=colnames(data_filt_b_cluster)
      rownames(annotation_row_b)=rownames(data_filt_b_cluster)
      #data_filt_a_cluster
      {
        b = pheatmap(data_filt_b_cluster,
                     color = colorRampPalette(color.key)(50),
                     border_color = NA,
                     cluster_cols = F,
                     cluster_rows = F,
                     #分隔
                     gaps_row = 
                       c(
                         base::nrow(data_filt_b_cluster1), 
                         base::nrow(data_filt_b_cluster1) + base::nrow(data_filt_b_cluster2),
                         base::nrow(data_filt_b_cluster1) + base::nrow(data_filt_b_cluster2) + base::nrow(data_filt_b_cluster3),
                         base::nrow(data_filt_b_cluster1) + base::nrow(data_filt_b_cluster2) + base::nrow(data_filt_b_cluster3) + base::nrow(data_filt_b_cluster4),
                         base::nrow(data_filt_b_cluster1) + base::nrow(data_filt_b_cluster2) + base::nrow(data_filt_b_cluster3) + base::nrow(data_filt_b_cluster4) + base::nrow(data_filt_b_cluster5)
                       ),
                     annotation_col = annotation_col,
                     annotation_row = annotation_row_b,
                     main = "Cluster about B",
                     #scale = "row",
                     show_rownames = F,show_colnames = T,
                     cellwidth = 50, cellheight = 0.5
        )
      }
    }
  }
  #画AB两组共有的图
  {
    #添加注释
    annotation_col_ab = data.frame(Grop= rep(c("Control","TreatA","TreatB"),c(1,5,5)),Dose = c(0,1:5,1:5))
    annotation_row_ab = data.frame(Cluster = rep(c("Max_CTRL_33(7%)","Max_A15_259(53%)","Max_A30_23(5%)","Max_A60_17(3%)","Max_A90_17(3%)","Max_A120_52(10%)","Max_B15_18(4%)","Max_B30_19(4%)","Max_B60_22(5%)","Max_B90_17(3%)","Max_B120_16(3%)"),c(33,259,23,17,17,52,18,19,22,17,16)))
    color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")#调色
    #运行没反应时，运行dev.off()
    #画ab组共有的基因热图
    {
      rownames(annotation_col_ab)=colnames(data_filt_ab_cluster)
      rownames(annotation_row_ab)=rownames(data_filt_ab_cluster)
      #data_filt_a_cluster
      {
        ab = pheatmap(data_filt_ab_cluster,
                      color = colorRampPalette(color.key)(50),
                      border_color = NA,
                      cluster_cols = F,
                      cluster_rows = F,
                      #分隔
                      gaps_row =c(33,
                                  33+259,
                                  33+259+23,
                                  33+259+23+17,
                                  33+259+23+17+17,
                                  33+259+23+17+17+52,
                                  33+259+23+17+17+52+18,
                                  33+259+23+17+17+52+18+19,
                                  33+259+23+17+17+52+18+19+22,
                                  33+259+23+17+17+52+18+19+22+17,
                                  33+259+23+17+17+52+18+19+22+17+16
                      ),
                      annotation_col = annotation_col_ab,
                      annotation_row = annotation_row_ab,
                      main = "Cluster about commonality of A and B",
                      #scale = "row",
                      show_rownames = F,show_colnames = T,
                      cellwidth = 100, cellheight = 0.5
        )
      }
    }
  }
  
}
##画venn图
{
  # Check R packages required here
  if(!requireNamespace("VennDiagram"))
  {
    BiocManager::install()
    BiocManager::install(c("cowplot","showtext"))
  }
  library(VennDiagram)
  l = list(A = rownames(data_filt_a_cluster),
           B = rownames(data_filt_b_cluster))
  v = venn.diagram(l,
                   #height = 4000 ,
                   #width = 4000 ,
                   #resolution=500,
                   #alpha = 0.1,#设置透明度
                   filename = NULL,#imagetype="tiff",
                   category = c("A Treat", "B Treat"),
                   lwd=1 ,#设置圆弧宽度
                   #lty=2 ,#设置圆弧线性
                   col=c('cornflowerblue','yellow'),
                   fill = c("cornflowerblue", "yellow"),
                   cat.col = c("black", "black"),#集合名称的显示颜色
                   label.col="black",
                   cex =1.3,#区域内部数字的字体大小 
                   cat.cex=1.3,#分类名称的字体大小
                   #cat.dist = 0.09,   #分类名称距离边的距离 实际调整 
                   #cat.just = list(c(-1, -1), c(1, 1)),  #分类名称的位置  ，圈内或者圈外
                   #ext.pos = 30,  #线的角度 默认是正上方12点位置 
                   #ext.dist = -0.05,   #外部线的距离  跟根据圆圈的大小适当调整
                   #ext.length = 0.85,  #外部线长度 
                   ext.line.lwd = 0.1,  #外部线的宽度 
                   main = "The commonality of Treat A and Treat B",
                   main.cex = 2
  )
  dev.off()
  grid.draw(v)
}
##UpSetR包绘图
{
  #1，表格形式，在R语言里就是数据框了。行表示元素，列表示数据集分配和额外信息。
  #2，元素名的集合(没见过，不知道。。)fromList
  #3，venneuler包引入的用于描述集合交集的向量fromExpression
  library(UpSetR)
  #定义函数
  fromExpression <- function(input)
  {
    input <- list(input)
    intersections <- lapply(input, function(x) strsplit(names(unlist(x)), "&"))
    intersections <- lapply(intersections[[1]], function(x) unlist(as.list(x)))
    sets <- unique(unlist(intersections))
    data <- na.omit(data.frame(matrix(NA, ncol = length(sets))))
    names(data) <- sets
    counts <- lapply(input, function(x) unlist(x))
    names(counts[[1]]) <- NULL
    counts[[1]] <- as.numeric(counts[[1]])
    
    for(i in seq(intersections)) {
      cols <- match(names(data), intersections[[i]])
      cols[!is.na(cols)] <- 1
      cols[is.na(cols)] <- 0
      cols <- rep(cols, times = counts[[1]][i])
      cols <- matrix(cols, ncol = length(sets), byrow = T)
      cols <- data.frame(cols)
      names(cols) <- sets
      data <- rbind(data, cols)
    }
    return(data)
  }
  #输入序列
  {
    input = c(
      "MAX_A15" = length(rownames(data_filt_a_cluster2)),
      "MAX_A30" = length(rownames(data_filt_a_cluster3)),
      "MAX_A60" = length(rownames(data_filt_a_cluster4)),
      "MAX_A90" = length(rownames(data_filt_a_cluster5)),
      "MAX_A120" = length(rownames(data_filt_a_cluster6)),
      "MAX_B15" = length(rownames(data_filt_b_cluster2)),
      "MAX_B30" = length(rownames(data_filt_b_cluster3)),
      "MAX_B60" = length(rownames(data_filt_b_cluster4)),
      "MAX_B90" = length(rownames(data_filt_b_cluster5)),
      "MAX_B120" = length(rownames(data_filt_b_cluster6)),
      "MAX_A15&MAX_B15" = length(intersect(rownames(data_filt_a_cluster2),rownames(data_filt_b_cluster2))),
      "MAX_A15&MAX_B30" = length(intersect(rownames(data_filt_a_cluster2),rownames(data_filt_b_cluster3))),
      "MAX_A15&MAX_B60" = length(intersect(rownames(data_filt_a_cluster2),rownames(data_filt_b_cluster4))),
      "MAX_A15&MAX_B90" = length(intersect(rownames(data_filt_a_cluster2),rownames(data_filt_b_cluster5))),
      "MAX_A15&MAX_B120" = length(intersect(rownames(data_filt_a_cluster2),rownames(data_filt_b_cluster6))),
      
      "MAX_A30&MAX_B15" = length(intersect(rownames(data_filt_a_cluster3),rownames(data_filt_b_cluster2))),
      "MAX_A30&MAX_B30" = length(intersect(rownames(data_filt_a_cluster3),rownames(data_filt_b_cluster3))),
      "MAX_A30&MAX_B60" = length(intersect(rownames(data_filt_a_cluster3),rownames(data_filt_b_cluster4))),
      "MAX_A30&MAX_B90" = length(intersect(rownames(data_filt_a_cluster3),rownames(data_filt_b_cluster5))),
      "MAX_A30&MAX_B120" = length(intersect(rownames(data_filt_a_cluster3),rownames(data_filt_b_cluster6))),
      
      "MAX_A60&MAX_B15" = length(intersect(rownames(data_filt_a_cluster4),rownames(data_filt_b_cluster2))),
      "MAX_A60&MAX_B30" = length(intersect(rownames(data_filt_a_cluster4),rownames(data_filt_b_cluster3))),
      "MAX_A60&MAX_B60" = length(intersect(rownames(data_filt_a_cluster4),rownames(data_filt_b_cluster4))),
      "MAX_A60&MAX_B90" = length(intersect(rownames(data_filt_a_cluster4),rownames(data_filt_b_cluster5))),
      "MAX_A60&MAX_B120" = length(intersect(rownames(data_filt_a_cluster4),rownames(data_filt_b_cluster6))),
      
      "MAX_A90&MAX_B15" = length(intersect(rownames(data_filt_a_cluster5),rownames(data_filt_b_cluster2))),
      "MAX_A90&MAX_B30" = length(intersect(rownames(data_filt_a_cluster5),rownames(data_filt_b_cluster3))),
      "MAX_A90&MAX_B60" = length(intersect(rownames(data_filt_a_cluster5),rownames(data_filt_b_cluster4))),
      "MAX_A90&MAX_B90" = length(intersect(rownames(data_filt_a_cluster5),rownames(data_filt_b_cluster5))),
      "MAX_A90&MAX_B120" = length(intersect(rownames(data_filt_a_cluster5),rownames(data_filt_b_cluster6))),
      
      "MAX_A120&MAX_B15" = length(intersect(rownames(data_filt_a_cluster6),rownames(data_filt_b_cluster2))),
      "MAX_A120&MAX_B30" = length(intersect(rownames(data_filt_a_cluster6),rownames(data_filt_b_cluster3))),
      "MAX_A120&MAX_B60" = length(intersect(rownames(data_filt_a_cluster6),rownames(data_filt_b_cluster4))),
      "MAX_A120&MAX_B90" = length(intersect(rownames(data_filt_a_cluster6),rownames(data_filt_b_cluster5))),
      "MAX_A120&MAX_B120" = length(intersect(rownames(data_filt_a_cluster6),rownames(data_filt_b_cluster6)))
    )
  }
  #转化
  data_2 <- fromExpression(input)
  upset(data_2,nsets = 10,nintersects = 26,
        text.scale = c(1.3, 1.3, 1, 1, 1.5, 1)) #六个数字，分别控制c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars))
  
  
}

#拼接
if(FALSE)#全段注释
{
  #安装并加载包
  {
    #BiocManager::install()
    #BiocManager::install(c("cowplot","showtext"))
    library(cowplot)
    library(sysfonts)
    library(showtextdb)
    library(showtext)
    library(grid)
    library(ggplotify)
  }
  #转化热图为ggplot格式
  {
    library(ggplotify)
    a = as.ggplot(a)
    b = as.ggplot(b)
    ab = as.ggplot(ab)
    v = as.ggplot(v)
  }
  #建画布
  {
    grid.newpage() # 新建画布
    layout_1 <- grid.layout(nrow = 2, ncol = 2,widths = c(10.7,80), 
                            heights = c(base::nrow(data_filt_a_cluster1) * 0.5,
                                        base::nrow(data_filt_a_cluster2) * 0.5,
                                        base::nrow(data_filt_a_cluster3) * 0.5,
                                        base::nrow(data_filt_a_cluster4) * 0.5,
                                        base::nrow(data_filt_a_cluster5) * 0.5,
                                        base::nrow(data_filt_a_cluster6) * 0.5
                            )
    ) # 分成4行1列共4个版块
    pushViewport(viewport(layout = layout_1)) # 推出分为4个版块的视窗
  }
  #打印图
  {
    print(a1, vp = viewport(layout.pos.row = 1,
                            layout.pos.col = 2
    )) # 将a1输出到第一行
    print(a2, vp = viewport(layout.pos.row = 2,
                            layout.pos.col = c(1,2)
    )) # 将a2输出到第二行
    print(a3, vp = viewport(layout.pos.row = 3,
                            layout.pos.col = c(1,2),
                            just = bottom
    )) # 将a3输出到第三行
    print(a4, vp = viewport(layout.pos.row = 4,
                            layout.pos.col = c(1,2)
    )) # 将a4输出到第四行
    print(a5, vp = viewport(layout.pos.row = 5,
                            layout.pos.col = c(1,2)
    )) # 将a5输出到第五行
    print(a6, vp = viewport(layout.pos.row = 6,
                            layout.pos.col = c(1,2)
    )) # 将a6输出到第六行
    # 添加主标题和分标题
  }
}
#加标签
if(FALSE)#全段注释
{
  #转化热图为ggblot格式
  略
  grid.text(label="        Clucster1 78(17%)
          Low to high", x = 0.39,y = 0.85,gp=gpar(col="black",fontsize = 12,draw=TRUE))
  grid.text(label="        Clucster2 211(46%)
          Low to high to lower", x = 0.39,y =0.6,gp=gpar(col="black",fontsize = 12,draw=TRUE,just = c("left", "top")))
  grid.text(label="        Clucster3 44(10%)
          high to low", x = 0.39,y =0.34,gp=gpar(col="black",fontsize = 12,draw=TRUE))
  grid.text(label="        Clucster4 122(27%)
          high to low to higher", x = 0.39,y =0.16,gp=gpar(col="black",fontsize = 12,draw=TRUE)) 
}

#dev.off()#关闭当前设备


####GO注释
#读入数据
data_cluster = read.csv("WT(0.5~1h)_MAX of fruitfly.csv",header = T,row.names = 3)
{
  #BiocManager::install("clusterProfiler")
  #BiocManager::install("topGO")
  #BiocManager::install("Rgraphviz")
  #BiocManager::install("pathview")
  #BiocManager::install("org.Dm.eg.db"),果蝇的包
  ##下载绘图包
  #BiocManager::install("Hmisc")
  #BiocManager::install("ggplot2")
  
  library(clusterProfiler)
  library(topGO)
  library(Rgraphviz)
  library(pathview)
  library(org.Dm.eg.db)
  library(ggplot2)
  library(Hmisc)
  genename = c(rownames(data_cluster))
  #查看所有的支持和可转化类型
  keytypes(org.Dm.eg.db)
  #将基因名的SYMBOL格式转化为ENTREZID格式，其中org.Dm.eg.db为fly数据库
  DEG.entrez_id = data.frame()
  DEG.entrez_id = bitr(genename,fromType = "ENSEMBLTRANS",toType = "FLYBASE",OrgDb = "org.Dm.eg.db")
  ##GO分析
  #分子生物学功能（Moleular Function,MF，该基因在分子层面的功能是什么，它催化哪些功能）,生物学过程（Biological Process,BP，该基因参与了哪些生物学过程）,细胞学组分（Cellular Components,CC,基因存在于哪里）
  erich.go.ALL=enrichGO(gene = DEG.entrez_id$ENTREZID,
                        OrgDb = org.Hs.eg.db,
                        ont = "ALL",
                        pvalueCutoff = 1,
                        qvalueCutoff = 1,
                        readable = T)
  #分别选子集
  erich.go.MF = erich.go.ALL[which(erich.go.ALL@result$ONTOLOGY =="MF",)]#which返回下标
  erich.go.BP = erich.go.ALL[which(erich.go.ALL@result$ONTOLOGY =="BP",)]
  erich.go.CC = erich.go.ALL[which(erich.go.ALL@result$ONTOLOGY =="CC",)]
  #化GeneRatio为小数
  erich.go.MF$GeneRatio = apply(data.frame(name=erich.go.MF$GeneRatio),1,function(x) eval(parse(text=x)))
  erich.go.BP$GeneRatio = apply(data.frame(name=erich.go.BP$GeneRatio),1,function(x) eval(parse(text=x)))
  erich.go.CC$GeneRatio = apply(data.frame(name=erich.go.CC$GeneRatio),1,function(x) eval(parse(text=x)))
  ###
  # Functions to draw plots
  DrawGOBubblePlot <- function(dat, category = "BP", top.number = 10, col="blue"){
    # Draw bubble plot using DAVID function enrichment results
    
    category = toupper(category)
    if (category == "BP"){
      main.title = "Biological Process"
    } else if (category == "CC"){
      main.title = "Cellular Components"
    } else if (category == "MF"){
      main.title = "Molecular Function"
    } #else if (category == "KEGG"){
    #main.title = "KEGG"}
    else {
      return("错了! 整个正确的过来 [category].")
    }
    dat1 = dat[c(1:top.number),c("ID","Description","GeneRatio","pvalue","geneID","Count")]#输入数据的第一行到第top.number行，"ID","Description","GeneRadio","pvalue","geneID"列
    #if(category == 'KEGG'){
    #dat1$Term = substr(dat1$Term,10,200)
    #}
    #else{
    #dat1$Description = substr(dat1$Description,12,200)
    #}
    dat1$Description = capitalize(dat1$Description)
    dat1$Description = factor(dat1$Description,levels=dat1$Description[length(dat1$Description):1])
    dat1$pvalue = -log10(dat1$pvalue)
    
    p = ggplot(dat1,aes(GeneRatio,Description)) +
      geom_point(aes(size=Count,colour=pvalue)) +
      scale_colour_gradient(low=col,high="red") + 
      labs(colour=expression(-log[10]("pvalue")),size="Gene counts",  
           x="Gene Ratio",y="",title=main.title) +
      theme_bw() +
      scale_x_continuous(limits = c(0,max(dat1$GeneRatio) * 1.2)) 
    
    return(p)
  }
  
}
# Read in data and generate the plots
# Biological Process
DrawGOBubblePlot(erich.go.BP,"BP",10,"blue")
DrawGOBubblePlot(erich.go.CC,"CC",10,"blue")
DrawGOBubblePlot(erich.go.MF,"MF",10,"blue")


