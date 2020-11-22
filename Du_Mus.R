#清空环境
rm(list=ls())

#读入数据为列表格式，以第一列为行名称,读入为row_data输出为data
{
  getwd()
  #设置存储路径
  setwd("E:/scientific_research/du_human/du_human/du_human/X101SC19051602-Z01-J010-B10-16-result/0.SupFile")
  #重新查看更改过后的当前路径
  getwd()
  row_data = read.csv("all_compare.csv",
                    header = T,
                    sep=",",
                    row.names = 1,
                    fill=TRUE,
                    check.names = T,
                  stringsAsFactors = FALSE#检查变量名（列）是否符合规则，不符合的话自动修改
  )
  a = c("A1_fpkm","A2_fpkm","A3_fpkm","B1_fpkm","B2_fpkm","B3_fpkm","C1_fpkm","C2_fpkm","C3_fpkm","D1_fpkm","D2_fpkm","D3_fpkm")
  data = row_data[,a]
  colnames(data)
}

#去掉列表中中表达量异常的基因输入输出都为data
{
  data[data==0] = NA
  data = na.omit(data)
}

#给列表重命名输入输出都为data
{
  #重命名
  {
    library(reshape)
    #使用该函数时变量名经常出问题。使用dplyr::rename重命名了变量，以便名称不再取决于变量，而取决于字符。后面如果再出现问题，跟换变量前后的位置即可
    data = dplyr::rename(data,c("A1" = A1_fpkm,
                                "A2" = A2_fpkm,
                                "A3" = A3_fpkm,
                                "B1" = B1_fpkm,
                                "B2" = B2_fpkm,
                                "B3" = B3_fpkm,
                                "C1" = C1_fpkm,
                                "C2" = C2_fpkm,
                                "C3" = C3_fpkm,
                                "D1" = D1_fpkm,
                                "D2" = D2_fpkm,
                                "D3" = D3_fpkm
    ))
    head(data,1)#查看首行
  }
}

# DEseq2进行筛选差异表达,要求读入为rowcount,不接受处理过的FPKM/TPM等；输入mycounts ；需要两个表格，分别是mycounts和colData;输出A_vs_control,B_vs_control,C_vs_control,D_vs_control
{
  if(!requireNamespace("tidyverse"))
  {
    BiocManager::install()
    BiocManager::install("tidyverse")
  }
  if(!requireNamespace("DESeq2"))
  {
    BiocManager::install()
    BiocManager::install("DESeq2")
  }
  library(tidyverse)
  library(DESeq2)
  #读入数据为列表格式，以第一列为行名称,读入为mycount_row输出为mycounts
  {
    getwd()
    #设置存储路径
    setwd("E:/scientific_research/du_human/du_human/du_human/X101SC19051602-Z01-J010-B10-16-result/0.SupFile")
    #重新查看更改过后的当前路径
    mycount_row = read.csv("all_compare.csv",
                           header = T,
                           sep=",",
                           row.names = 1,
                           fill=TRUE,
                           check.names = T,
                           stringsAsFactors = FALSE#检查变量名（列）是否符合规则，不符合的话自动修改
    )
    a = c("A1_count","A2_count","A3_count","B1_count","B2_count","B3_count","C1_count","C2_count","C3_count","D1_count","D2_count","D3_count")
    mycount_row = mycount_row[,a]
    colnames(mycount_row)
  }
  #去掉列表中中表达量异常的基因输入输出都为mycount_row
  {
    mycount_row[mycount_row==0] = NA
    mycount_row = na.omit(mycount_row)
  }
  #给列表重命名输入输出都为mycounts
  {
    #重命名
    {
      library(reshape)
      #使用该函数时变量名经常出问题。使用dplyr::rename重命名了变量，以便名称不再取决于变量，而取决于字符。后面如果再出现问题，跟换变量前后的位置即可
      mycount_row = dplyr::rename(mycount_row,c("A1" = A1_count,
                                                "A2" = A2_count,
                                                "A3" = A3_count,
                                                "B1" = B1_count,
                                                "B2" = B2_count,
                                                "B3" = B3_count,
                                                "C1" = C1_count,
                                                "C2" = C2_count,
                                                "C3" = C3_count,
                                                "D1" = D1_count,
                                                "D2" = D2_count,
                                                "D3" = D3_count
      ))
      head(mycount_row,1)#查看首行
    }
  }
  #输入mycount_row，输出mycounts制作对照组,小鼠数据有对照组和处理组，各三个重复,这里把它们的均值作为对照组
  mycounts = mycount_row
  mycounts$mean_1 = round(apply(mycounts[,1:12],1,mean)) 
  mycounts$mean_2 = round(apply(mycounts[,1:12],1,mean))
  mycounts$mean_3 = round(apply(mycounts[,1:12],1,mean))
  #给列重排序,默认情况下，R会根据字母表顺序排列因子型变量，排在最前面的因子作为对照,输入输出都是mycounts
  {
    library(dplyr,warn.conflicts = F)
    mycounts = mycounts %>% dplyr::select(colnames(mycounts[13:15]),colnames(mycounts[1:12]))
  }
  # 制作colData，这一步很关键，要明白condition这里是因子，不是样本名称；样本信息colData，每一行对应一个样本，行名与countData的样本顺序一一对应，列为各种分组信息。
  condition <- factor(c(rep("control",3),rep("treat_A",3),rep("treat_B",3),rep("treat_C",3),rep("treat_D",3)), levels = c("control","treat_A","treat_B","treat_C","treat_D"))
  colData <- data.frame(row.names=colnames(mycounts), condition)
  #构建dds对象,开始DESeq流程，dds=DESeqDataSet Object; design参数设置公式，默认情况下，此包中的函数将使用公式中的最后一个变量来构建结果表和绘图；默认情况下，R会根据字母表顺序排列因子型变量，排在最前面的因子作为对照。
  dds = DESeqDataSetFromMatrix(countData = mycounts, colData = colData, design= ~ condition)
  #筛选行的和大于5的基因
  dds <- dds[rowSums(counts(dds)) > 50, ] 
  #校正数据，评估离散程度等：
  dds = DESeq(dds)
  #提取结果
  #res <- results(dds)默认使用样本信息的最后一个因子与第一个因子进行比
  resultsNames = resultsNames(dds)[-1]
  
  #####两两比对
  #x_vs_control
  {
    i = 1
    library(dplyr,warn.conflicts = F)
    for (i in resultsNames) 
    {
      filename = paste("C:/Users/Administrator/Desktop/Gene list of DU/DESeq2/",substr(i,17,nchar(i)),sep = "")
      filename = paste(filename,".csv",sep = "")
      assign(substr(i,17,nchar(i)),assign(substr(i,17,nchar(i)),results(dds,name = i)) %>% 
            (function(x){x = x[order(x$pvalue),]}) %>% 
            (function(x){subset(x, padj < 0.01 & abs(log2FoldChange) > 1)}) %>%
            (function(x){write.csv(x,file = filename)}) 
             )
      A = assign(substr(i,17,nchar(i)),read.csv(file = filename,row.names = 1))
      B = row_data[row.names(A),]
      A = cbind(A,B[,"gene_name"])
      row.names(A) = A$`B[, "gene_name"]`
      #A_vs_control = A_vs_control[,-length(A_vs_control)]
      A$Group = "not-significant"
      A$Group[which((A$padj<0.01)&(A$log2FoldChange > 1.5))] = "up-regulated"
      A$Group[which((A$padj<0.01)&(A$log2FoldChange < -1.5))] = "down-regulated"
      A = A[order(A$padj),]
      #高表达的基因中，选择padj最小的10个
      up_genes = head(A$`B[, "gene_name"]`[which(A$Group == "up-regulated")], 10)
      #低表达的基因中，选择padj最小的10个
      down_genes = head(A$`B[, "gene_name"]`[which(A$Group == "down-regulated")], 10)
      #将up_genes和down. genes合并，并加入到Label中
      deg.top10.genes = c(as.character(up_genes),as.character(down_genes))
      A$Lable = ""
      A$Lable[as.numeric(match(deg.top10.genes,A$`B[, "gene_name"]`))]= deg.top10.genes
      A$logP = -log10(A$padj)
      colnames(A)[7] = "gene_name"
      write.csv(A,file = filename)
      assign(substr(i,17,nchar(i)),A)
    }
  }
}
#火山图展示富集效果
input = D_vs_control
{
  if(!requireNamespace("ggpubr"))
  {
    BiocManager::install()
    BiocManager::install("ggpubr")
  }
  if(!requireNamespace("ggthemes"))
  {
    BiocManager::install()
    BiocManager::install("ggthemes")
  }
  library(ggpubr)
  library(ggthemes)
  ggscatter(input,x = "log2FoldChange",y = "logP",
            #title = input,
            color = "Group",
            palette = c("#2f5688","#BBBBBB","#CC0000"),
            size = 1,
            label = input$Lable,
            font.label = 10,#字体大小
            repel = T,
            xlab = "log2FoldChange",
            ylab = "-log10(Adjust P-value)",)+theme_base()+
    geom_hline(yintercept = 1.30,linetype="dashed")+
    geom_vline(xintercept = c(-1.5,1.5),linetype="dashed")
  
}
#ggsave()

###GSEA（基因集富集分析），GSEA无需先做差异分析，会保留更多更全面的关键信息。可以帮助我们找到那些——差异不是很明显但基因差异趋势很一致——的功能基因集
{
  #读取基因名，读取文件夹"C:/Users/Administrator/Desktop/Gene list of DU/csv"中的每一个文件，取其基因名，产一个包含所有名称的l列表
  {
    #读入数据为列表格式，以第一列为行名称,读入为row_data输出为data
    {
      getwd()
      #设置存储路径
      setwd("E:/scientific_research/du_human/du_human/du_human/X101SC19051602-Z01-J010-B10-16-result/0.SupFile")
      #重新查看更改过后的当前路径
      getwd()
      row_data = read.csv("all_compare.csv",
                          header = T,
                          sep=",",
                          row.names = 1,
                          fill=TRUE,
                          check.names = T,
                          stringsAsFactors = FALSE#检查变量名（列）是否符合规则，不符合的话自动修改
      )
      a = c("A1_fpkm","A2_fpkm","A3_fpkm","B1_fpkm","B2_fpkm","B3_fpkm","C1_fpkm","C2_fpkm","C3_fpkm","D1_fpkm","D2_fpkm","D3_fpkm")
      data = row_data[,a]
      colnames(data)
    }
    #去掉列表中中表达量异常的基因输入输出都为data
    {
      data[data==0] = NA
      data = na.omit(data)
    }
    #给列表重命名输入输出都为data
    {
      #重命名
      {
        library(reshape)
        #使用该函数时变量名经常出问题。使用dplyr::rename重命名了变量，以便名称不再取决于变量，而取决于字符。后面如果再出现问题，跟换变量前后的位置即可
        data = dplyr::rename(data,c("A1" = A1_fpkm,
                                    "A2" = A2_fpkm,
                                    "A3" = A3_fpkm,
                                    "B1" = B1_fpkm,
                                    "B2" = B2_fpkm,
                                    "B3" = B3_fpkm,
                                    "C1" = C1_fpkm,
                                    "C2" = C2_fpkm,
                                    "C3" = C3_fpkm,
                                    "D1" = D1_fpkm,
                                    "D2" = D2_fpkm,
                                    "D3" = D3_fpkm
        ))
        head(data,1)#查看首行
      }
    }
    #取均值，输入data输出newdata
    #成功#！要比对重复的时勿执行
    {
      i = 1
      newdata = data
      while (i < length(colnames(data))){
        newdata$X <- rowMeans(data[,c(i,i+1,i+2)]) 
        a = c(substr(colnames(data[i]),1,nchar(colnames(data[i]))-1))
        colnames(newdata)[length(newdata)] = a
        i = i+3
      }
      newdata = newdata[,-which(colnames(newdata) %in% colnames(data))]
    }
    #筛选每一行和大于50的行，输入输出都是newdata
    {
      #library(ballgown)
      newdata = subset(newdata,apply(newdata,1,sum)>50,genomesubset=TRUE)
    }
    #加载包
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
      library(enrichplot)
    }
    #处理数据，产出log2FC以及,利用gseGO()以及gseaplot()进行分析以及画图。设置了一个循环，将产生的csv文件以及png图片保存
    for(i in 1:4)
    {
      for (j in 1:4) 
      {
        if (i != j & i < j) 
        {
          k = c()
          k = paste(colnames(newdata[i]),"vs",sep = "")
          k = paste(k,colnames(newdata[j]),sep = "")
          name = c()
          name =paste("C:/Users/Administrator/Desktop/Gene list of DU/GSEA",k,sep = "/") 
          name = paste(name,".csv",sep = "")
          newdata[,c(i,j)] %>% (function(x){write.csv(x,file = name)})
          A = read.csv(name,
                       header = T,
                       sep=",",
                       row.names = 1,
                       fill=TRUE,
                       check.names = T,
                       stringsAsFactors = FALSE#检查变量名（列）是否符合规则，不符合的话自动修改
          )
          A$log2_fold_change = log2(A[1]/A[2])
          B = row_data[row.names(A),]
          A = cbind(A,B[,"gene_name"])
          colnames(A)[length(A)]= "SYMBOL"
          #转ID从"SYMBOL"到"ENTREZID"
          DEG.entrez_id = bitr(A$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mm.eg.db")
          A = merge(A,DEG.entrez_id,by = "SYMBOL",all = F)
          A = as.matrix(A[order(A$log2_fold_change,decreasing = T),])
          #保存并从新读取，不然其中的FC列一直为list格式，无法进行下一步处理
          A %>% (function(x){write.csv(x,file = name)})
          A = read.csv(name,
                       header = T,
                       sep=",",
                       row.names = 1,
                       fill=TRUE,
                       check.names = T,
                       stringsAsFactors = FALSE#检查变量名（列）是否符合规则，不符合的话自动修改
          )
          #设置需要输入的数据，两列包括基因名和FC
          geneList = A$log2_fold_change
          names(geneList) = as.character(A$ENTREZID)
          name = c()
          name =paste("C:/Users/Administrator/Desktop/Gene list of DU/GSEA",k,sep = "/") 
          geneList <- sort(geneList,decreasing = T)
          gseGO.res.ALL = gseGO(geneList, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont="ALL", minGSSize = 10, maxGSSize = 500, pvalueCutoff=1)
          #再次保存
          name = c()
          name =paste("C:/Users/Administrator/Desktop/Gene list of DU/GSEA",k,sep = "/") 
          name = paste(name,".csv",sep = "")
          data.frame(gseGO.res.ALL) %>% (function(x){write.csv(x,file = name)})
          #读取后面Go分析的结果，输出图片并保存ggsave()
          name_csv = c()
          name_csv =paste("C:/Users/Administrator/Desktop/Gene list of DU/go",k,sep = "/") 
          name_csv = paste(name_csv,".csv",sep = "")
          name_png = c()
          name_png =paste("C:/Users/Administrator/Desktop/Gene list of DU/GSEA",k,sep = "/") 
          B = read.csv(file = name_csv,
                       header = T,
                       sep=",",
                       row.names = 1,
                       fill=TRUE,
                       check.names = F,
                       stringsAsFactors = FALSE#
                       )
          B_MF = rownames(B[1:5,])
          B_BP = rownames(B[6:10,])
          B_CC = rownames(B[11:15,])
          name_png_MF = paste(name_png,"_MF",sep = "")
          name_png_MF = paste(name_png_MF,".png",sep = "")
          gseaplot2(gseGO.res.ALL,geneSetID = B_MF, pvalue_table = F,title = paste(k,"MF",sep = "-"),base_size = 40,ES_geom = "line",subplots = 1:3,) %>% (function(x){ggsave(x,filename = name_png_MF,width = 27,height = 20)})
          name_png_BP = paste(name_png,"_BP",sep = "")
          name_png_BP = paste(name_png_BP,".png",sep = "")
          gseaplot2(gseGO.res.ALL,geneSetID = B_BP, pvalue_table = F,title = paste(k,"BP",sep = "-"),base_size = 40,ES_geom = "line",subplots = 1:3,) %>% (function(x){ggsave(x,filename = name_png_BP,width = 27,height = 20)})
          name_png_CC = paste(name_png,"_CC",sep = "")
          name_png_CC = paste(name_png_CC,".png",sep = "")
          gseaplot2(gseGO.res.ALL,geneSetID = B_CC, pvalue_table = F,title = paste(k,"CC",sep = "-"),base_size = 40,ES_geom = "line",subplots = 1:3,) %>% (function(x){ggsave(x,filename = name_png_CC,width = 27,height = 20)})
        }
      }
      i = i+1
    }
  }
}

  


#GO分析———在线，先做总表的GO分析，此处为MAXrow为1以上的,文件名为RNA in nemateode mutants.csv
if (FALSE)
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
    #a = data.frame(number = c("1：————","2：————","Value：————","Screening Criteria：———— ","Big Cluster(four)：————","Small Cluster(twenty four)：————"),'————Specification————' = c("This specification is for #Differentially expressed RNA in nemateode mutants.zip#","The data are not raw data,but were screened for genes that were considered the most likely to be differentially expressed","log2FPKM(FPKM: Fragments Per Kilobase of exon model per Million mapped fragments)","MAX per row > 20,MUT(max)/CTRL>2 or CTRL/MUT(min)>2","Sort by log2FPKM from large to small ,and cluster according to the maximum value of each group","4321 4231 3421 3241 2431 2341 4312 4132 3412 3142 1432 1342 4213 4123 2413 2143 1423 1243 3214 3124 2314 2134 1324 1234,(1 for CTRL,2 for MUT_qx666,3 for MUT_epg.5,4 for MUT_dpy.7,and The four positions represent the size of the value,From left to right means bigger and bigger"),row.names = 1,)
    #write.table(a,file = "Specification.txt",quote = F)
  }
  #-----
}

#取均值，输出newdata
#成功#！要比对重复的时勿执行
{
  i = 1
  newdata = data
  while (i < length(colnames(data))){
    newdata$X <- rowMeans(data[,c(i,i+1,i+2)]) 
    a = c(substr(colnames(data[i]),1,nchar(colnames(data[i]))-1))
    colnames(newdata)[length(newdata)] = a
    i = i+3
  }
  newdata = newdata[,-which(colnames(newdata) %in% colnames(data))]
}

#筛选每一行和大于50的行，输入输出都是newdata
{
  #library(ballgown)
  newdata = subset(newdata,apply(newdata,1,sum)>50,genomesubset=TRUE)
}

#总：差异性表达输入newdata,输出newdata_filt
if (FALSE)
{
  newdata_FC = newdata
  newdata_FC$FC_up = apply(newdata,1,max)/apply(newdata,1,mean)
  newdata_FC$FC_down = apply(newdata,1,mean)/apply(newdata,1,min)
  newdata_filt = subset(newdata_FC,newdata_FC$FC_up>2| newdata_FC$FC_down>2)
  newdata_filt = newdata_filt[,c(-length(newdata_filt),-(length(newdata_filt)-1))]
}
#取Log2
#newdata_filt = log2(newdata_filt)

#分开差异表达，输入newdata,输出newdata_filt
data_x = newdata[,c(1,2)]
data_x = newdata[,c(1,3)]
data_x = newdata[,c(1,4)]
data_x = newdata[,c(2,3)]
data_x = newdata[,c(2,4)]
data_x = newdata[,c(3,4)]
{
  newdata_FC = data_x
  newdata_FC$FC= newdata_FC[1]/newdata_FC[2]
  newdata_filt = subset(newdata_FC,newdata_FC$FC>2| newdata_FC$FC<0.5)
  newdata_filt = newdata_filt[,-length(newdata_filt)]
}
#取Log2
newdata_filt = log2(newdata_filt)

#聚类II按最大值最优聚类 kmeans法输出结果cluster_data
#分开聚类，输入newdata_filt输出cluster_data
data_x = newdata_filt
library("ggplot2")
library ("factoextra")
library ("cluster")
{
  {
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
#集合聚类，输入newdata_filt输出cluster_data
#data_x = newdata_filt
if (FALSE)
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
if (FALSE)
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
#if (FALSE)
{
  ##准备工作
  {
    #读取基因名，读取文件夹"C:/Users/Administrator/Desktop/Gene list of DU/csv"中的每一个文件，取其基因名，产一个包含所有名称的l列表
    {
      i = 1
      j = 1
      l = list()
      name = c()
      for (i in list.files("C:/Users/Administrator/Desktop/Gene list of DU/csv")) 
      {
        a = paste("C:/Users/Administrator/Desktop/Gene list of DU/csv",#第一个向量
                  i,#第二个向量
                  sep = "/" #连接符
        )
        a = read.csv(a,
                     header = T,
                     sep=",",
                     row.names = 1,
                     fill=TRUE,
                     check.names = T,
                     stringsAsFactors = FALSE#检查变量名（列）是否符合规则，不符合的话自动修改
        )
        name = c(name,substr(i,1,nchar(i)-4))
        l[[j]] = rownames(a)
        names(l)=name
        j = j+1
      }
    }
  }
  ##开始作图
  #加载包
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
  }
  i = 1
  for (i in 1:6) {
    genename = unlist(l[i])
  title = paste("GO analyze about Mus of",names(l[i]),sep = " ")
  {
    keytypes(org.Mm.eg.db)
    DEG.entrez_id = bitr(genename,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mm.eg.db")#需要转化为ENTREZID，才可以进行后续分析
    #查看所有的支持和可转化类型
    #keytypes(org.Mm.eg.db)
    #将基因名的SYMBOL格式转化为ENTREZID格式，其中org.Mm.eg.db为小鼠数据库
    #DEG.entrez_id = bitr(genename,fromType = "SYMBOL",toType = "GOALL",OrgDb = "org.Hs.eg.db")
    ##GO分析
    #分子生物学功能（Moleular Function,MF，该基因在分子层面的功能是什么，它催化哪些功能）,生物学过程（Biological Process,BP，该基因参与了哪些生物学过程）,细胞学组分（Cellular Components,CC,基因存在于哪里）
    erich.go.ALL=enrichGO(gene = DEG.entrez_id$ENTREZID,
                          OrgDb = org.Mm.eg.db,
                          ont = "ALL",
                          pvalueCutoff = 1,
                          qvalueCutoff = 1,
                          readable = T)
    erich.go.ALL@result$GeneRatio = apply(data.frame(name=erich.go.ALL@result$GeneRatio),1,function(x) eval(parse(text=x)))
    #分别选子集
    #MF
    erich.go.MF = erich.go.ALL[which(erich.go.ALL@result$ONTOLOGY =="MF",)][1:5,]#which返回下标,因为erich.go默认是根据Pvalue升序排序，故默认选取前五个
    erich.go.BP = erich.go.ALL[which(erich.go.ALL@result$ONTOLOGY =="BP",)][1:5,]
    erich.go.CC = erich.go.ALL[which(erich.go.ALL@result$ONTOLOGY =="CC",)][1:5,]
    erich.go.ALL_result = rbind(erich.go.MF,erich.go.BP,erich.go.CC)
    erich.go.ALL_result$Description = capitalize(erich.go.ALL_result$Description)
    erich.go.ALL_result$Description = factor(erich.go.ALL_result$Description,levels=erich.go.ALL_result$Description[length(erich.go.ALL_result$Description):1])
    erich.go.ALL_result$pvalue = -log10(erich.go.ALL_result$pvalue)
    
    #绘图
    p = ggplot(erich.go.ALL_result,
               aes(GeneRatio,Description),
               fil) +
      geom_point(aes(size=Count,colour=pvalue)) +
      scale_colour_gradient(low="blue",high="red") + 
      labs(colour=expression(("pvalue")),
           size="Gene counts",  
           x="Gene Ratio",y="",title=title) +
      theme_bw() +
      scale_x_continuous(limits = c(0,max(erich.go.ALL_result$GeneRatio) * 1.2))
  }#输出p
  #保存
  library
  file = c()
  file = paste("C:/Users/Administrator/Desktop/Gene list of DU/go",#第一个向量
               names(l[i]),#第二个向量
               sep = "/" ) %>% (function(x){file = paste(x,".png",sep = "")})
  ggsave(filename = file,plot = p)
  file = c()
  file = paste("C:/Users/Administrator/Desktop/Gene list of DU/go",#第一个向量
               names(l[i]),#第二个向量
               sep = "/" ) %>% (function(x){file = paste(x,".csv",sep = "")})
  write.csv(erich.go.ALL_result,file = file)
  
  }
}

#ALL数据读出.csv
{
  write.csv(cluster_data,file = "all high-quality RNA of mut.csv")
}
#Discrepant数据读出.csv
{
  file_name = paste("C:/Users/Administrator/Desktop/Gene list of DU/",colnames(cluster_data[1]),sep = "")
  file_name = paste(file_name,"vs",sep = "")
  file_name = paste(file_name,colnames(cluster_data[2]),sep = "")
  file_name = paste(file_name,".csv",sep = "")
  csv = row_data[rownames(cluster_data),]
  csv = cbind(cluster_data,csv)
  colname = c("gene_name",colnames(csv[,c(1,2)]))
  csv = csv[,colname]
  write.csv(csv,file = file_name,row.names = F)
}

##画venn图
#前期准备，产出一个命名为l的列表，并于C:/Users/Administrator/Desktop/Gene list of DU文件夹中产出基因名的.txt文件
{
  options(stringsAsFactors = F)
  options(warn = -1)
  # Check R packages required here
  i = 1
  j = 1
  l = list()
  name = c()
  for (i in list.files("C:/Users/Administrator/Desktop/Gene list of DU/csv")) 
  {
    a = paste("C:/Users/Administrator/Desktop/Gene list of DU/csv",#第一个向量
              i,#第二个向量
              sep = "/" #连接符
    )
    a = read.csv(a,
                 header = T,
                 sep=",",
                 row.names = NULL,
                 fill=TRUE,
                 check.names = T,
                 stringsAsFactors = FALSE#检查变量名（列）是否符合规则，不符合的话自动修改
    )
    file_name = paste("C:/Users/Administrator/Desktop/Gene list of DU/txt/",colnames(a[2]),sep = "")
    file_name = paste(file_name,"vs",sep = "")
    file_name = paste(file_name,colnames(a[3]),sep = "")
    file_name = paste(file_name,".txt",sep = "")
    txt = a[1]
    write.table(txt,file = file_name, quote = F,sep = "\t",col.names = F,row.names = F)
    name = c(name,substr(i,1,nchar(i)-4))
    l[[j]] = rownames(a)
    names(l)=name
    j = j+1
  }
}
#VennDiagram包画图，！！！仅可以画5个及其以下的样本！！！本例使用BMKCloud在线（http://www.biocloud.net/），中的venn工具（https://console.biocloud.net/static/index.html#/drawtools/intoDrawTools/venn）
{
  if(!requireNamespace("VennDiagram"))
  {
    BiocManager::install()
    BiocManager::install(c("cowplot","showtext"))
  }
  library(VennDiagram)
  v = venn.diagram(l,
                   height = 4000 ,
                   width = 4000 ,
                   #resolution=500,
                   #alpha = 0.1,#设置透明度
                   filename = NULL,#imagetype="tiff",
                   #category = name,
                   lwd=1 ,#设置圆弧宽度
                   #lty=2 ,#设置圆弧线性
                   col="black",
                   fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
                   #cat.col = c("black", "black"),#集合名称的显示颜色
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


#绘小热图
main = c("Differentially expressed RNA in sample of A and B ")
main = c("Differentially expressed RNA in sample of A and C ")
main = c("Differentially expressed RNA in sample of A and D ")
main = c("Differentially expressed RNA in sample of B and C ")
main = c("Differentially expressed RNA in sample of B and D ")
main = c("Differentially expressed RNA in sample of C and D ")
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
                 cellwidth = 100, cellheight = 0.5,
                 fontsize = 20,
                 
        )
      }
    }
  }
  
}

#绘大热图#输入cluster_data
main = c("Discrepant gene list in human(mean value)")
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
                 cellwidth = 50, cellheight = 0.5,
                 fontsize = 20
        )
      }
    }
  }
  
}





