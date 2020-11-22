#��ջ���
rm(list=ls())

#��������Ϊ�б���ʽ���Ե�һ��Ϊ������,����Ϊrow_data���Ϊdata
{
  getwd()
  #���ô洢·��
  setwd("E:/scientific_research/du_human/du_human/du_human/X101SC19051602-Z01-J010-B10-16-result/0.SupFile")
  #���²鿴���Ĺ���ĵ�ǰ·��
  getwd()
  row_data = read.csv("all_compare.csv",
                    header = T,
                    sep=",",
                    row.names = 1,
                    fill=TRUE,
                    check.names = T,
                  stringsAsFactors = FALSE#�����������У��Ƿ���Ϲ��򣬲����ϵĻ��Զ��޸�
  )
  a = c("A1_fpkm","A2_fpkm","A3_fpkm","B1_fpkm","B2_fpkm","B3_fpkm","C1_fpkm","C2_fpkm","C3_fpkm","D1_fpkm","D2_fpkm","D3_fpkm")
  data = row_data[,a]
  colnames(data)
}

#ȥ���б����б������쳣�Ļ������������Ϊdata
{
  data[data==0] = NA
  data = na.omit(data)
}

#���б����������������Ϊdata
{
  #������
  {
    library(reshape)
    #ʹ�øú���ʱ���������������⡣ʹ��dplyr::rename�������˱������Ա����Ʋ���ȡ���ڱ�������ȡ�����ַ�����������ٳ������⣬��������ǰ���λ�ü���
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
    head(data,1)#�鿴����
  }
}

# DEseq2����ɸѡ�������,Ҫ�����Ϊrowcount,�����ܴ�������FPKM/TPM�ȣ�����mycounts ����Ҫ�������񣬷ֱ���mycounts��colData;���A_vs_control,B_vs_control,C_vs_control,D_vs_control
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
  #��������Ϊ�б���ʽ���Ե�һ��Ϊ������,����Ϊmycount_row���Ϊmycounts
  {
    getwd()
    #���ô洢·��
    setwd("E:/scientific_research/du_human/du_human/du_human/X101SC19051602-Z01-J010-B10-16-result/0.SupFile")
    #���²鿴���Ĺ���ĵ�ǰ·��
    mycount_row = read.csv("all_compare.csv",
                           header = T,
                           sep=",",
                           row.names = 1,
                           fill=TRUE,
                           check.names = T,
                           stringsAsFactors = FALSE#�����������У��Ƿ���Ϲ��򣬲����ϵĻ��Զ��޸�
    )
    a = c("A1_count","A2_count","A3_count","B1_count","B2_count","B3_count","C1_count","C2_count","C3_count","D1_count","D2_count","D3_count")
    mycount_row = mycount_row[,a]
    colnames(mycount_row)
  }
  #ȥ���б����б������쳣�Ļ������������Ϊmycount_row
  {
    mycount_row[mycount_row==0] = NA
    mycount_row = na.omit(mycount_row)
  }
  #���б����������������Ϊmycounts
  {
    #������
    {
      library(reshape)
      #ʹ�øú���ʱ���������������⡣ʹ��dplyr::rename�������˱������Ա����Ʋ���ȡ���ڱ�������ȡ�����ַ�����������ٳ������⣬��������ǰ���λ�ü���
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
      head(mycount_row,1)#�鿴����
    }
  }
  #����mycount_row�����mycounts����������,С�������ж�����ʹ����飬�������ظ�,��������ǵľ�ֵ��Ϊ������
  mycounts = mycount_row
  mycounts$mean_1 = round(apply(mycounts[,1:12],1,mean)) 
  mycounts$mean_2 = round(apply(mycounts[,1:12],1,mean))
  mycounts$mean_3 = round(apply(mycounts[,1:12],1,mean))
  #����������,Ĭ������£�R�������ĸ��˳�����������ͱ�����������ǰ���������Ϊ����,�����������mycounts
  {
    library(dplyr,warn.conflicts = F)
    mycounts = mycounts %>% dplyr::select(colnames(mycounts[13:15]),colnames(mycounts[1:12]))
  }
  # ����colData����һ���ܹؼ���Ҫ����condition���������ӣ������������ƣ�������ϢcolData��ÿһ�ж�Ӧһ��������������countData������˳��һһ��Ӧ����Ϊ���ַ�����Ϣ��
  condition <- factor(c(rep("control",3),rep("treat_A",3),rep("treat_B",3),rep("treat_C",3),rep("treat_D",3)), levels = c("control","treat_A","treat_B","treat_C","treat_D"))
  colData <- data.frame(row.names=colnames(mycounts), condition)
  #����dds����,��ʼDESeq���̣�dds=DESeqDataSet Object; design�������ù�ʽ��Ĭ������£��˰��еĺ�����ʹ�ù�ʽ�е����һ������������������ͻ�ͼ��Ĭ������£�R�������ĸ��˳�����������ͱ�����������ǰ���������Ϊ���ա�
  dds = DESeqDataSetFromMatrix(countData = mycounts, colData = colData, design= ~ condition)
  #ɸѡ�еĺʹ���5�Ļ���
  dds <- dds[rowSums(counts(dds)) > 50, ] 
  #У�����ݣ�������ɢ�̶ȵȣ�
  dds = DESeq(dds)
  #��ȡ���
  #res <- results(dds)Ĭ��ʹ��������Ϣ�����һ���������һ�����ӽ��б�
  resultsNames = resultsNames(dds)[-1]
  
  #####�����ȶ�
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
      #�߱���Ļ����У�ѡ��padj��С��10��
      up_genes = head(A$`B[, "gene_name"]`[which(A$Group == "up-regulated")], 10)
      #�ͱ���Ļ����У�ѡ��padj��С��10��
      down_genes = head(A$`B[, "gene_name"]`[which(A$Group == "down-regulated")], 10)
      #��up_genes��down. genes�ϲ��������뵽Label��
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
#��ɽͼչʾ����Ч��
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
            font.label = 10,#�����С
            repel = T,
            xlab = "log2FoldChange",
            ylab = "-log10(Adjust P-value)",)+theme_base()+
    geom_hline(yintercept = 1.30,linetype="dashed")+
    geom_vline(xintercept = c(-1.5,1.5),linetype="dashed")
  
}
#ggsave()

###GSEA�����򼯸�����������GSEA������������������ᱣ�������ȫ��Ĺؼ���Ϣ�����԰��������ҵ���Щ�������첻�Ǻ����Ե�����������ƺ�һ�¡����Ĺ��ܻ���
{
  #��ȡ����������ȡ�ļ���"C:/Users/Administrator/Desktop/Gene list of DU/csv"�е�ÿһ���ļ���ȡ�����������һ�������������Ƶ�l�б�
  {
    #��������Ϊ�б���ʽ���Ե�һ��Ϊ������,����Ϊrow_data���Ϊdata
    {
      getwd()
      #���ô洢·��
      setwd("E:/scientific_research/du_human/du_human/du_human/X101SC19051602-Z01-J010-B10-16-result/0.SupFile")
      #���²鿴���Ĺ���ĵ�ǰ·��
      getwd()
      row_data = read.csv("all_compare.csv",
                          header = T,
                          sep=",",
                          row.names = 1,
                          fill=TRUE,
                          check.names = T,
                          stringsAsFactors = FALSE#�����������У��Ƿ���Ϲ��򣬲����ϵĻ��Զ��޸�
      )
      a = c("A1_fpkm","A2_fpkm","A3_fpkm","B1_fpkm","B2_fpkm","B3_fpkm","C1_fpkm","C2_fpkm","C3_fpkm","D1_fpkm","D2_fpkm","D3_fpkm")
      data = row_data[,a]
      colnames(data)
    }
    #ȥ���б����б������쳣�Ļ������������Ϊdata
    {
      data[data==0] = NA
      data = na.omit(data)
    }
    #���б����������������Ϊdata
    {
      #������
      {
        library(reshape)
        #ʹ�øú���ʱ���������������⡣ʹ��dplyr::rename�������˱������Ա����Ʋ���ȡ���ڱ�������ȡ�����ַ�����������ٳ������⣬��������ǰ���λ�ü���
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
        head(data,1)#�鿴����
      }
    }
    #ȡ��ֵ������data���newdata
    #�ɹ�#��Ҫ�ȶ��ظ���ʱ��ִ��
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
    #ɸѡÿһ�кʹ���50���У������������newdata
    {
      #library(ballgown)
      newdata = subset(newdata,apply(newdata,1,sum)>50,genomesubset=TRUE)
    }
    #���ذ�
    {
      #BiocManager::install("clusterProfiler")
      #BiocManager::install("topGO")
      #BiocManager::install("Rgraphviz")
      #BiocManager::install("pathview")
      #BiocManager::install("org.Mm.eg.db")
      ##���ػ�ͼ��
      #BiocManager::install("Hmisc")
      #BiocManager::install("ggplot2")
      
      library(clusterProfiler)
      library(topGO)
      library(Rgraphviz)
      library(pathview)
      library(org.Mm.eg.db)#С���#�������EntreZ���л����ע�ͷ�������ͬ������ѡȡ��ͬ��DroBD��,�ο�https://www.omicsclass.com/article/262����http://www.bioconductor.org/packages/release/BiocViews.html#___AnnotationData
      library(ggplot2)
      library(Hmisc)
      library(enrichplot)
    }
    #�������ݣ�����log2FC�Լ�,����gseGO()�Լ�gseaplot()���з����Լ���ͼ��������һ��ѭ������������csv�ļ��Լ�pngͼƬ����
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
                       stringsAsFactors = FALSE#�����������У��Ƿ���Ϲ��򣬲����ϵĻ��Զ��޸�
          )
          A$log2_fold_change = log2(A[1]/A[2])
          B = row_data[row.names(A),]
          A = cbind(A,B[,"gene_name"])
          colnames(A)[length(A)]= "SYMBOL"
          #תID��"SYMBOL"��"ENTREZID"
          DEG.entrez_id = bitr(A$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mm.eg.db")
          A = merge(A,DEG.entrez_id,by = "SYMBOL",all = F)
          A = as.matrix(A[order(A$log2_fold_change,decreasing = T),])
          #���沢���¶�ȡ����Ȼ���е�FC��һֱΪlist��ʽ���޷�������һ������
          A %>% (function(x){write.csv(x,file = name)})
          A = read.csv(name,
                       header = T,
                       sep=",",
                       row.names = 1,
                       fill=TRUE,
                       check.names = T,
                       stringsAsFactors = FALSE#�����������У��Ƿ���Ϲ��򣬲����ϵĻ��Զ��޸�
          )
          #������Ҫ��������ݣ����а�����������FC
          geneList = A$log2_fold_change
          names(geneList) = as.character(A$ENTREZID)
          name = c()
          name =paste("C:/Users/Administrator/Desktop/Gene list of DU/GSEA",k,sep = "/") 
          geneList <- sort(geneList,decreasing = T)
          gseGO.res.ALL = gseGO(geneList, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont="ALL", minGSSize = 10, maxGSSize = 500, pvalueCutoff=1)
          #�ٴα���
          name = c()
          name =paste("C:/Users/Administrator/Desktop/Gene list of DU/GSEA",k,sep = "/") 
          name = paste(name,".csv",sep = "")
          data.frame(gseGO.res.ALL) %>% (function(x){write.csv(x,file = name)})
          #��ȡ����Go�����Ľ�������ͼƬ������ggsave()
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

  


#GO�������������ߣ������ܱ���GO�������˴�ΪMAXrowΪ1���ϵ�,�ļ���ΪRNA in nemateode mutants.csv
if (FALSE)
{
  write.table(data.frame(geneid = row.names(data)),file = "geneid1.txt",row.names = F,col.names = F,quote = F)#����������Ϊ���������������Ϊ.txt��ʽ����Ҫ��������Ҫ����,��Ҫ����
  #-----
  #�����ݻ���id��Դ��Wormbase���ݿ�,��ʽΪ������WBGene00009237��
  #��������DAVID  https://david.ncifcrf.gov/,����ճ���ύ��Identifierѡ��WORMBASE_GENE_ID.������ѡ�����е� Gene ID Conversion Tool��ת��OFFICAL_GENE_SYMBOL������ͨ�����ƣ����ؼ��ɣ�����Ϊall.txt
  #���벢������
  {
    GO_gene_name = read.table("all.txt",header = T,sep ="\t",#Tap�ָ�
                              fill=TRUE,#�ò����ɺ������ݲ�����Ĵ���ֱ�Ӷ�ȡ
                              stringsAsFactors = FALSE,quote = ""#��Ҫ����
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
  #ɾ��ȱʧֵ
  {
    GO_gene_name[GO_gene_name==""] = NA
    GO_gene_name = na.omit(GO_gene_name)
  }
  #�ϲ�
  {
    cluster_go = merge(cluster,GO_gene_name,sort = F,all = T)
    rownames(cluster_go) = cluster_go[,"Wormbase_id"]
    #����
    #������
    library(dplyr,warn.conflicts = F)
    colnames(cluster_go)
    cluster_go = cluster_go %>% select("Wormbase_id","Gene Name","Species","Discription","CTRL","MUT_qx666","MUT_epg.5","MUT_dpy.7","MAX_MUT","MIN_MUT")
    
  }
  #�����б���ע��
  {
    write.csv(cluster_go,file = "RNA in nemateode mutants.csv")
    #a = data.frame(number = c("1����������","2����������","Value����������","Screening Criteria���������� ","Big Cluster(four)����������","Small Cluster(twenty four)����������"),'��������Specification��������' = c("This specification is for #Differentially expressed RNA in nemateode mutants.zip#","The data are not raw data,but were screened for genes that were considered the most likely to be differentially expressed","log2FPKM(FPKM: Fragments Per Kilobase of exon model per Million mapped fragments)","MAX per row > 20,MUT(max)/CTRL>2 or CTRL/MUT(min)>2","Sort by log2FPKM from large to small ,and cluster according to the maximum value of each group","4321 4231 3421 3241 2431 2341 4312 4132 3412 3142 1432 1342 4213 4123 2413 2143 1423 1243 3214 3124 2314 2134 1324 1234,(1 for CTRL,2 for MUT_qx666,3 for MUT_epg.5,4 for MUT_dpy.7,and The four positions represent the size of the value,From left to right means bigger and bigger"),row.names = 1,)
    #write.table(a,file = "Specification.txt",quote = F)
  }
  #-----
}

#ȡ��ֵ�����newdata
#�ɹ�#��Ҫ�ȶ��ظ���ʱ��ִ��
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

#ɸѡÿһ�кʹ���50���У������������newdata
{
  #library(ballgown)
  newdata = subset(newdata,apply(newdata,1,sum)>50,genomesubset=TRUE)
}

#�ܣ������Ա�������newdata,���newdata_filt
if (FALSE)
{
  newdata_FC = newdata
  newdata_FC$FC_up = apply(newdata,1,max)/apply(newdata,1,mean)
  newdata_FC$FC_down = apply(newdata,1,mean)/apply(newdata,1,min)
  newdata_filt = subset(newdata_FC,newdata_FC$FC_up>2| newdata_FC$FC_down>2)
  newdata_filt = newdata_filt[,c(-length(newdata_filt),-(length(newdata_filt)-1))]
}
#ȡLog2
#newdata_filt = log2(newdata_filt)

#�ֿ�����������newdata,���newdata_filt
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
#ȡLog2
newdata_filt = log2(newdata_filt)

#����II�����ֵ���ž��� kmeans��������cluster_data
#�ֿ����࣬����newdata_filt���cluster_data
data_x = newdata_filt
library("ggplot2")
library ("factoextra")
library ("cluster")
{
  {
    #�����ֵ���ž���
    {
      i = 1
      cluster_data = data.frame()
      while(i <= length(data_x))
      {
        a =  subset(data_x,data_x[,i] > data_x[,colnames(data_x)[-i]])
        #ѡ����෽��
        #�ھ�����������ľ���ķ�������(���ڶ��塰���롱��ͳ����)�����Ծ��룺manhattan;ŷ�Ͼ��룺euclideanĬ��;�ɿɷ�˹�����룺minkowski;�б�ѩ����룺chebyshev;���Ͼ��룺mahalanobis;���Ͼ��룺canberra
        dist_method = "euclidean"
        #����������ķ�����
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
#���Ͼ��࣬����newdata_filt���cluster_data
#data_x = newdata_filt
if (FALSE)
{
  #BiocManager::install()
  #BiocManager::install(c("factoextra","cluster"))
  library ("factoextra")
  library ("cluster")
  #�����ֵ���ž���
  {
    {
      #BiocManager::install()
      #BiocManager::install(c("factoextra","cluster"))
      library ("factoextra")
      library ("cluster")
      #�����ֵ���ž���
      {
        i = 1
        cluster_data = data.frame()
        while(i <= length(data_x))
        {
          a =  subset(data_x,data_x[,i] > apply(data_x[,colnames(data_x)[-i]],1,max))
          #ѡ����෽��
          #�ھ�����������ľ���ķ�������(���ڶ��塰���롱��ͳ����)�����Ծ��룺manhattan;ŷ�Ͼ��룺euclideanĬ��;�ɿɷ�˹�����룺minkowski;�б�ѩ����룺chebyshev;���Ͼ��룺mahalanobis;���Ͼ��룺canberra
          dist_method = "euclidean"
          #����������ķ�������ƽ������average;���ķ���centroid;�м���뷨��median;����뷨��complete(Ĭ��);��̾��뷨��single;���ƽ���ͷ���ward;�ܶȹ��Ʒ���density
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

#�����б���GO��ʽ
if (FALSE)
{
  a = read.csv("RNA in nemateode mutants.csv",
               header = T,
               row.names = 1,
               check.names = F#�������������У��Ƿ���Ϲ���
  )
  a = a[row.names(cluster_data),]
  write.csv(a,file = "Differentially expressed RNA in nemateode mutants(kmeans cluster).csv",row.names = F)
}

####GOע��
#if (FALSE)
{
  ##׼������
  {
    #��ȡ����������ȡ�ļ���"C:/Users/Administrator/Desktop/Gene list of DU/csv"�е�ÿһ���ļ���ȡ�����������һ�������������Ƶ�l�б�
    {
      i = 1
      j = 1
      l = list()
      name = c()
      for (i in list.files("C:/Users/Administrator/Desktop/Gene list of DU/csv")) 
      {
        a = paste("C:/Users/Administrator/Desktop/Gene list of DU/csv",#��һ������
                  i,#�ڶ�������
                  sep = "/" #���ӷ�
        )
        a = read.csv(a,
                     header = T,
                     sep=",",
                     row.names = 1,
                     fill=TRUE,
                     check.names = T,
                     stringsAsFactors = FALSE#�����������У��Ƿ���Ϲ��򣬲����ϵĻ��Զ��޸�
        )
        name = c(name,substr(i,1,nchar(i)-4))
        l[[j]] = rownames(a)
        names(l)=name
        j = j+1
      }
    }
  }
  ##��ʼ��ͼ
  #���ذ�
  {
    #BiocManager::install("clusterProfiler")
    #BiocManager::install("topGO")
    #BiocManager::install("Rgraphviz")
    #BiocManager::install("pathview")
    #BiocManager::install("org.Mm.eg.db")
    ##���ػ�ͼ��
    #BiocManager::install("Hmisc")
    #BiocManager::install("ggplot2")
    
    library(clusterProfiler)
    library(topGO)
    library(Rgraphviz)
    library(pathview)
    library(org.Mm.eg.db)#С���#�������EntreZ���л����ע�ͷ�������ͬ������ѡȡ��ͬ��DroBD��,�ο�https://www.omicsclass.com/article/262����http://www.bioconductor.org/packages/release/BiocViews.html#___AnnotationData
    library(ggplot2)
    library(Hmisc)
  }
  i = 1
  for (i in 1:6) {
    genename = unlist(l[i])
  title = paste("GO analyze about Mus of",names(l[i]),sep = " ")
  {
    keytypes(org.Mm.eg.db)
    DEG.entrez_id = bitr(genename,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mm.eg.db")#��Ҫת��ΪENTREZID���ſ��Խ��к�������
    #�鿴���е�֧�ֺͿ�ת������
    #keytypes(org.Mm.eg.db)
    #����������SYMBOL��ʽת��ΪENTREZID��ʽ������org.Mm.eg.dbΪС�����ݿ�
    #DEG.entrez_id = bitr(genename,fromType = "SYMBOL",toType = "GOALL",OrgDb = "org.Hs.eg.db")
    ##GO����
    #��������ѧ���ܣ�Moleular Function,MF���û����ڷ��Ӳ���Ĺ�����ʲô�����߻���Щ���ܣ�,����ѧ���̣�Biological Process,BP���û����������Щ����ѧ���̣�,ϸ��ѧ��֣�Cellular Components,CC,������������
    erich.go.ALL=enrichGO(gene = DEG.entrez_id$ENTREZID,
                          OrgDb = org.Mm.eg.db,
                          ont = "ALL",
                          pvalueCutoff = 1,
                          qvalueCutoff = 1,
                          readable = T)
    erich.go.ALL@result$GeneRatio = apply(data.frame(name=erich.go.ALL@result$GeneRatio),1,function(x) eval(parse(text=x)))
    #�ֱ�ѡ�Ӽ�
    #MF
    erich.go.MF = erich.go.ALL[which(erich.go.ALL@result$ONTOLOGY =="MF",)][1:5,]#which�����±�,��Ϊerich.goĬ���Ǹ���Pvalue�������򣬹�Ĭ��ѡȡǰ���
    erich.go.BP = erich.go.ALL[which(erich.go.ALL@result$ONTOLOGY =="BP",)][1:5,]
    erich.go.CC = erich.go.ALL[which(erich.go.ALL@result$ONTOLOGY =="CC",)][1:5,]
    erich.go.ALL_result = rbind(erich.go.MF,erich.go.BP,erich.go.CC)
    erich.go.ALL_result$Description = capitalize(erich.go.ALL_result$Description)
    erich.go.ALL_result$Description = factor(erich.go.ALL_result$Description,levels=erich.go.ALL_result$Description[length(erich.go.ALL_result$Description):1])
    erich.go.ALL_result$pvalue = -log10(erich.go.ALL_result$pvalue)
    
    #��ͼ
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
  }#���p
  #����
  library
  file = c()
  file = paste("C:/Users/Administrator/Desktop/Gene list of DU/go",#��һ������
               names(l[i]),#�ڶ�������
               sep = "/" ) %>% (function(x){file = paste(x,".png",sep = "")})
  ggsave(filename = file,plot = p)
  file = c()
  file = paste("C:/Users/Administrator/Desktop/Gene list of DU/go",#��һ������
               names(l[i]),#�ڶ�������
               sep = "/" ) %>% (function(x){file = paste(x,".csv",sep = "")})
  write.csv(erich.go.ALL_result,file = file)
  
  }
}

#ALL���ݶ���.csv
{
  write.csv(cluster_data,file = "all high-quality RNA of mut.csv")
}
#Discrepant���ݶ���.csv
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

##��vennͼ
#ǰ��׼��������һ������Ϊl���б�������C:/Users/Administrator/Desktop/Gene list of DU�ļ����в�����������.txt�ļ�
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
    a = paste("C:/Users/Administrator/Desktop/Gene list of DU/csv",#��һ������
              i,#�ڶ�������
              sep = "/" #���ӷ�
    )
    a = read.csv(a,
                 header = T,
                 sep=",",
                 row.names = NULL,
                 fill=TRUE,
                 check.names = T,
                 stringsAsFactors = FALSE#�����������У��Ƿ���Ϲ��򣬲����ϵĻ��Զ��޸�
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
#VennDiagram����ͼ�������������Ի�5���������µ���������������ʹ��BMKCloud���ߣ�http://www.biocloud.net/�����е�venn���ߣ�https://console.biocloud.net/static/index.html#/drawtools/intoDrawTools/venn��
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
                   #alpha = 0.1,#����͸����
                   filename = NULL,#imagetype="tiff",
                   #category = name,
                   lwd=1 ,#����Բ������
                   #lty=2 ,#����Բ������
                   col="black",
                   fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
                   #cat.col = c("black", "black"),#�������Ƶ���ʾ��ɫ
                   label.col="black",
                   cex =1.3,#�����ڲ����ֵ������С 
                   cat.cex=1.3,#�������Ƶ������С
                   #cat.dist = 0.09,   #�������ƾ���ߵľ��� ʵ�ʵ��� 
                   #cat.just = list(c(-1, -1), c(1, 1)),  #�������Ƶ�λ��  ��Ȧ�ڻ���Ȧ��
                   #ext.pos = 30,  #�ߵĽǶ� Ĭ�������Ϸ�12��λ�� 
                   #ext.dist = -0.05,   #�ⲿ�ߵľ���  ������ԲȦ�Ĵ�С�ʵ�����
                   #ext.length = 0.85,  #�ⲿ�߳��� 
                   ext.line.lwd = 0.1,  #�ⲿ�ߵĿ��� 
                   main = "The commonality of Treat A and Treat B",
                   main.cex = 2
  )
  dev.off()
  grid.draw(v)
  
}


#��С��ͼ
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
  while (i <= length(colnames(cluster_data))) #���Ʋ�ͬ��ͼ�ǵø�����
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
  #�������ͼ
  {
    #����ע��
    annotation_col= data.frame(Group = colnames(cluster_data))
    annotation_row = data.frame(Cluster = rep(Cluster,length_Cluster))
    color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")#��ɫ
    #����û��Ӧʱ������dev.off()
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
                 #�ָ�
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

#�����ͼ#����cluster_data
main = c("Discrepant gene list in human(mean value)")
{
  library(pheatmap)
  gaps_row = c()
  sum = 0
  Cluster = c()
  length_Cluster = c()
  i = 1
  while (i <= length(colnames(cluster_data))) #���Ʋ�ͬ��ͼ�ǵø�����
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
  #�������ͼ
  {
    #����ע��
    #annotation_col= data.frame(Group = colnames(cluster_data))
    annotation_row = data.frame(Cluster = rep(Cluster,length_Cluster))
    color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")#��ɫ
    #����û��Ӧʱ������dev.off()
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
                 #�ָ�
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




