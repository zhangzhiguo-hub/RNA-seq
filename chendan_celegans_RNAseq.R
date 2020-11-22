#清空环境
rm(list=ls())


###二选一
###ballgown包处理数据，maxrow > 20;输出bg_filt_trans
{
  #加载ballgown包
  library(ballgown)
  library(stringr)
  #查看当前路径
  getwd()
  #设置存储路径
  setwd("E:/scientific_research/Transcriptome_analysis/chendan_celegans_RNAseq/mut")
  #重新查看更改过后的当前路径
  getwd()
  #读取定量分析的表型数据（批量操作型）
  list.files("./ballgown")#查看文件名
  #字符串连接
  a = paste(rep("ballgown",length(list.files("./ballgown"))),#第一个向量
            list.files("./ballgown"),#第二个向量
            sep = "/" #连接符
  )
  bg = ballgown(dataDir = "./ballgown",samples = a )
  #过滤掉表达量低的基因,即每一行所代表的基因，其后的各值都为0的
  bg_filt=as.data.frame(gexpr(bg))
  bg_filt_trans = subset(bg_filt,apply(bg_filt,1,max)>20,genomesubset=TRUE)
  #gexpr()选取的结果为基因,texpr()选取的结果为转录组
  #feature:在那个水平计算，可选“gene”,"transcript","exon","intron"
  #协变量covariate是对因变量有影响的变量，它不是研究者研究的自变量，不为实验者操作，但对试验有影响。
  #主变量adjustvars
  #getFC可以指定输出结果显示组间表达量的foldchange(差异倍数，FC=Exp x1/x2)
  #meas:measurement采用哪种计算方式
  #mod0=model.matrix(~?, pData(bg_filt)$treat)
  #mod=model.matrix(~?, pData(bg_filt)$treat)
  #results_tanscripts=stattest(bg_filt,feature = "transcript",covariat = ?,adjustvars = ?,mod0 = mod0,mod = mod,meas = "FPKM")本例中无重复，此法不可用来做差异化分析。
}

###ballgown包处理数据，maxrow > 1;输出bg_filt_trans
{
  #加载ballgown包
  library(ballgown)
  library(stringr)
  #查看当前路径
  getwd()
  #设置存储路径
  setwd("E:/scientific_research/Transcriptome_analysis/chendan_celegans_RNAseq/mut")
  #重新查看更改过后的当前路径
  getwd()
  #读取定量分析的表型数据（批量操作型）
  list.files("./ballgown")#查看文件名
  #字符串连接
  a = paste(rep("ballgown",length(list.files("./ballgown"))),#第一个向量
            list.files("./ballgown"),#第二个向量
            sep = "/" #连接符
  )
  bg = ballgown(dataDir = "./ballgown",samples = a )
  #过滤掉表达量低的基因,即每一行所代表的基因，其后的各值都为0的
  bg_filt=as.data.frame(gexpr(bg))
  bg_filt_trans = subset(bg_filt,apply(bg_filt,1,max)>1,genomesubset=TRUE)
  #gexpr()选取的结果为基因,texpr()选取的结果为转录组
  #feature:在那个水平计算，可选“gene”,"transcript","exon","intron"
  #协变量covariate是对因变量有影响的变量，它不是研究者研究的自变量，不为实验者操作，但对试验有影响。
  #主变量adjustvars
  #getFC可以指定输出结果显示组间表达量的foldchange(差异倍数，FC=Exp x1/x2)
  #meas:measurement采用哪种计算方式
  #mod0=model.matrix(~?, pData(bg_filt)$treat)
  #mod=model.matrix(~?, pData(bg_filt)$treat)
  #results_tanscripts=stattest(bg_filt,feature = "transcript",covariat = ?,adjustvars = ?,mod0 = mod0,mod = mod,meas = "FPKM")本例中无重复，此法不可用来做差异化分析。
}

####读出并保存数据mut.txt
{
  write.table(bg_filt_trans,file = "mut.csv",sep = ",",row.names = T) #列表保存为.csv格式
  write.table(bg_filt_trans,file = "mut.txt",sep = " ",row.names = T) #列表保存为.txt格式
}

#读入数据为列表格式，以第一列为行名称
{
  getwd()
  #设置存储路径
  setwd("E:/scientific_research/Transcriptome_analysis/chendan_celegans_RNAseq/mut")
  #重新查看更改过后的当前路径
  getwd()
  data = read.table("mut.txt",
                    header = T,
                    row.names = 1,
                    check.names = TRUE#检查变量名（列）是否符合规则，不符合的话自动修改
  )
}

#去掉列表中中表达量异常的基因
{
  data[data==0] = NA
  data = na.omit(data)
}

#给列表重命名并取
{
  #重命名
  {
    library(reshape)
    #使用该函数时变量名经常出问题。使用dplyr::rename重命名了变量，以便名称不再取决于变量，而取决于字符。后面如果再出现问题，跟换变量前后的位置即可
    data = dplyr::rename(data,c("CTRL" = FPKM.4_N2_12_h_Ctr,
                                "MUT_qx666" = FPKM.6_qx666_Day_1,
                                "MUT_epg.5" = FPKM.7_epg.5_Day_1,
                                "MUT_dpy.7" = FPKM.8_dpy.7_Day_1
    ))
    head(data,1)#查看首行
  }
}

#GO分析———在线，先做总表的GO分析，此处为MAXrow为1以上的,文件名为RNA in nemateode mutants.csv
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
    cluster_go[cluster_go == NA] = c("Not")
    #排序
    #列排序
    library(dplyr,warn.conflicts = F)
    colnames(cluster_go)
    cluster_go = cluster_go %>% select("Wormbase_id","Gene Name","Species","Discription","CTRL","MUT_qx666","MUT_epg.5","MUT_dpy.7","MAX_MUT","MIN_MUT")
    
  }
  #读出列表及注释
  {
    write.csv(cluster_go,file = "RNA in nemateode mutants.csv")
    #a = data.frame(number = c("1：————","2：————","Value：————","Screening Criteria：———— ","Big Cluster(four)：————","Small Cluster(twenty four)：————"),'————Specification————' = c("This specification is for #Differentially expressed RNA in nemateode mutants.zip#","The data are not raw data,but were screened for genes that were considered the most likely to be differentially expressed","log2FPKM(FPKM: Fragments Per Kilobase of exon model per Million mapped fragments)","MAX per row > 20,MUT(max)/CTRL>2 or CTRL/MUT(min)>2","Sort by log2FPKM from large to small ,and cluster according to the maximum value of each group","4321 4231 3421 3241 2431 2341 4312 4132 3412 3142 1432 1342 4213 4123 2413 2143 1423 1243 3214 3124 2314 2134 1324 1234,(1 for CTRL,2 for MUT_qx666,3 for MUT_epg.5,4 for MUT_dpy.7,and The four positions represent the size of the value,From left to right means bigger and bigger"),row.names = 1,)
    #write.table(a,file = "Specification.txt",quote = F)
  }
  #-----
}

#差异表达分析
#需要原始数据时勿操作此步
{
  #求出处理组每个基因的表达的最小值和最大值
  {
    #求出处理组每个基因的表达的最小值
    treat = c(colnames(data[-1]))
    ctrl = c(colnames(data[1]))
    data$treat_min = as.matrix(apply(data[,treat],1,min))
    #求出处理组每个基因的表达的最大值
    data$treat_max = as.matrix(apply(data[,treat],1,max))
  }
  #筛选差异表达基因data_filt
  {
    data$log2FC1 = log2((data$treat_max)/(data$CTRL))
    data$log2FC2 = log2((data$CTRL)/(data$treat_min))
    data_filt = subset(data,data$log2FC1>1|data$log2FC2>1)
    data = data_filt[,c(ctrl,treat)]
  }
}
 
data = log2(data)
 
#聚类I 二十四个小聚类
{
  #生成1:4的所有随机组合并由小到大排列,num向量
  #order()函数的排序原则，（最小，次小，。。最大）
  {
    num = c()
    n = 1
    for(i in c(1,2,3,4)){
      for (j in c(1,2,3,4)) {
        for (k in c(1,2,3,4)) {
          for (l in c(1,2,3,4)) {
            if (length(unique(c(i,j,k,l))) == 4)
              num = c(num,paste(i,j,k,l,sep = ""))
            num = as.numeric(num)
            num = num[order(-num)]
            num = num[c(grep("...1",num),grep("...2",num),grep("...3",num),grep("...4",num))]
            n=n+1
          }
        }
      }
    }
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

#聚类II按最大值最优聚类 kmeansa法输出结果cluster_data
#分开聚类
data_x = data[,c(1,2)]
data_x = data[,c(1,3)]
data_x = data[,c(1,4)]
{
  {
    #BiocManager::install()
    #BiocManager::install(c("factoextra","cluster"))
    library ("factoextra")
    library ("cluster")
    #按最大值最优聚类
    {
      i = 1
      cluster_data = data.frame()
      while(i <= length(data_x))
      {
        a =  subset(data_x,data_x[,i] > data_x[,colnames(data_x)[-i]])
        #选择聚类方法
        #在聚类中求两点的距离的方法包括(用于定义“距离”的统计量)，绝对距离：manhattan;欧氏距离：euclidean默认;闵可夫斯基距离：minkowski;切比雪夫距离：chebyshev;马氏距离：mahalanobis;蓝氏距离：canberra
        dist_method = "euclidean"
        #计算类间距离的方法：
        cluster_method = "complete"
        b = dist(a,method = dist_method)
        b = hclust(b,method = cluster_method)
        order = b$order
        a = a[order,]
        cluster_data = rbind(cluster_data,a)
        i <- i + 1
      }
    }
  }
}

#集合聚类
data_x = data
{
  #BiocManager::install()
  #BiocManager::install(c("factoextra","cluster"))
  library ("factoextra")
  library ("cluster")
  #按最大值最优聚类
  {
    {
      #BiocManager::install()
      #BiocManager::install(c("factoextra","cluster"))
      library ("factoextra")
      library ("cluster")
      #按最大值最优聚类
      {
        i = 1
        cluster_data = data.frame()
        while(i <= length(data_x))
        {
          a =  subset(data_x,data_x[,i] > apply(data_x[,colnames(data_x)[-i]],1,max))
          #选择聚类方法
          #在聚类中求两点的距离的方法包括(用于定义“距离”的统计量)，绝对距离：manhattan;欧氏距离：euclidean默认;闵可夫斯基距离：minkowski;切比雪夫距离：chebyshev;马氏距离：mahalanobis;蓝氏距离：canberra
          dist_method = "euclidean"
          #计算类间距离的方法：类平均法：average;重心法：centroid;中间距离法：median;最长距离法：complete(默认);最短距离法：single;离差平方和法：ward;密度估计法：density
          cluster_method = "complete"
          b = dist(a,method = dist_method)
          b = hclust(b,method = cluster_method)
          order = b$order
          a = a[order,]
          cluster_data = rbind(cluster_data,a)
          i <- i + 1
        }
      }
    }
  }
}

#读出列表的GO格式
{
  a = read.csv("RNA in nemateode mutants.csv",
               header = T,
               row.names = 1,
               check.names = F#不检查变量名（列）是否符合规则
  )
  a = a[row.names(cluster_data),]
  write.csv(a,file = "Differentially expressed RNA in nemateode mutants(kmeans cluster).csv",row.names = F)
}

####GO注释
if (FALSE)
{
  #BiocManager::install("clusterProfiler")
  #BiocManager::install("topGO")
  #BiocManager::install("Rgraphviz")
  #BiocManager::install("pathview")
  #BiocManager::install("org.Mm.eg.db")
  ##下载绘图包
  #BiocManager::install("Hmisc")
  #BiocManager::install("ggplot2")
  
  library(clusterProfiler)
  library(topGO)
  library(Rgraphviz)
  library(pathview)
  library(org.Mm.eg.db)#小鼠的#方便根据EntreZ进行基因的注释分析。不同的五种选取不同的DroBD包,参考https://www.omicsclass.com/article/262或者http://www.bioconductor.org/packages/release/BiocViews.html#___AnnotationData
  library(ggplot2)
  library(Hmisc)
  genename = c(rownames(cluster_data))
  #查看所有的支持和可转化类型
  keytypes(org.Hs.eg.db)
  #将基因名的SYMBOL格式转化为ENTREZID格式，其中org.Hs.eg.db为人类数据库
  DEG.entrez_id = bitr(genename,fromType = "SYMBOL",toType = "GOALL",OrgDb = "org.Hs.eg.db")
  ##GO分析
  #分子生物学功能（Moleular Function,MF，该基因在分子层面的功能是什么，它催化哪些功能）,生物学过程（Biological Process,BP，该基因参与了哪些生物学过程）,细胞学组分（Cellular Components,CC,基因存在于哪里）
  erich.go.ALL=enrichGO(gene = DEG.entrez_id$ENTREZID,
                        OrgDb = org.Mm.eg.db,
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
  # Read in data and generate the plots
  # Biological Process
  DrawGOBubblePlot(erich.go.BP,"BP",10,"blue")
  DrawGOBubblePlot(erich.go.CC,"CC",10,"blue")
  DrawGOBubblePlot(erich.go.MF,"MF",10,"blue")
  DrawGOBubblePlot = function(erich.go.ALL,category = "Mf",)
    
    dotplot(erich.go.ALL)
  plotGOgraph(erich.go.ALL)
}

#ALL数据读出.csv
{
  write.csv(cluster_data,file = "all high-quality RNA of mut.csv")
}
#Discrepant数据读出.csv
{
  write.csv(cluster_data,file = "Discrepant of mut.csv")
}


#绘小热图
main = c("Differentially expressed RNA in nemateode mutants of qx666 ")
main = c("Differentially expressed RNA in nemateode mutants of epg.5 ")
main = c("Differentially expressed RNA in nemateode mutants of dpy.7 ")
{
  library(pheatmap)
  gaps_row = c()
  sum = 0
  Cluster = c()
  length_Cluster = c()
  i = 1
  while (i <= length(colnames(cluster_data))) #绘制不同的图记得改名儿
  {
    a = nrow(subset(cluster_data,cluster_data[,i] >= cluster_data[,colnames(cluster_data)[-i]]))
    b = paste("max",
              colnames(cluster_data)[i],
              sep = "_")
    b = paste (b,
               a,
               sep = " ")
    b = paste (b,
               round(a/nrow(cluster_data)*100),
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
  #画聚类的图
  {
    #添加注释
    annotation_col= data.frame(Group = colnames(cluster_data))
    annotation_row = data.frame(Cluster = rep(Cluster,length_Cluster))
    color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")#调色
    #运行没反应时，运行dev.off()
    {
      rownames(annotation_col)=colnames(cluster_data)
      rownames(annotation_row)=rownames(cluster_data)
      #cluster_data
      {
         pheatmap(cluster_data,
                 color = colorRampPalette(color.key)(50),
                 border_color = NA,
                 cluster_cols = F,
                 cluster_rows = F,
                 #分隔
                 gaps_row = gaps_row,
                 #annotation_col = annotation_col,
                 annotation_row = annotation_row,
                 main = main,
                 #scale = "row",
                 show_rownames= F,#show_colnames = T,
                 cellwidth = 80, cellheight = 0.5,
                 fontsize = 20,
                 
        )
      }
    }
  }
  
}

#绘大热图
main = c("Differentially expressed RNA in all nemateode mutants")
{
  library(pheatmap)
  gaps_row = c()
  sum = 0
  Cluster = c()
  length_Cluster = c()
  i = 1
  while (i <= length(colnames(cluster_data))) #绘制不同的图记得改名儿
  {
    a = nrow(subset(cluster_data,cluster_data[,i] > apply(cluster_data[,colnames(cluster_data)[-i]],1,max)))
    b = paste("max",
              colnames(cluster_data)[i],
              sep = "_")
    b = paste (b,
               a,
               sep = " ")
    b = paste (b,
               round(a/nrow(cluster_data)*100),
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
  #画聚类的图
  {
    #添加注释
    #annotation_col= data.frame(Group = colnames(cluster_data))
    annotation_row = data.frame(Cluster = rep(Cluster,length_Cluster))
    color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")#调色
    #运行没反应时，运行dev.off()
    {
      #rownames(annotation_col)=colnames(cluster_data)
      rownames(annotation_row)=rownames(cluster_data)
      #cluster_data
      {
        pheatmap(cluster_data,
                 color = colorRampPalette(color.key)(50),
                 border_color = NA,
                 cluster_cols = F,
                 cluster_rows = F,
                 #分隔
                 gaps_row = gaps_row,
                 #annotation_col = annotation_col,
                 annotation_row = annotation_row,
                 main = main,
                 #scale = "row",
                 show_rownames = F,show_colnames = T,
                 cellwidth = 80, cellheight = 0.8,
                 fontsize = 20
        )
      }
    }
  }
  
}





