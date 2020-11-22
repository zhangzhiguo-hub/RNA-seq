#��ջ���
rm(list=ls())


###��ѡһ
###ballgown���������ݣ�maxrow > 20;���bg_filt_trans
{
  #����ballgown��
  library(ballgown)
  library(stringr)
  #�鿴��ǰ·��
  getwd()
  #���ô洢·��
  setwd("E:/scientific_research/Transcriptome_analysis/chendan_celegans_RNAseq/mut")
  #���²鿴���Ĺ���ĵ�ǰ·��
  getwd()
  #��ȡ���������ı������ݣ����������ͣ�
  list.files("./ballgown")#�鿴�ļ���
  #�ַ�������
  a = paste(rep("ballgown",length(list.files("./ballgown"))),#��һ������
            list.files("./ballgown"),#�ڶ�������
            sep = "/" #���ӷ�
  )
  bg = ballgown(dataDir = "./ballgown",samples = a )
  #���˵��������͵Ļ���,��ÿһ���������Ļ������ĸ�ֵ��Ϊ0��
  bg_filt=as.data.frame(gexpr(bg))
  bg_filt_trans = subset(bg_filt,apply(bg_filt,1,max)>20,genomesubset=TRUE)
  #gexpr()ѡȡ�Ľ��Ϊ����,texpr()ѡȡ�Ľ��Ϊת¼��
  #feature:���Ǹ�ˮƽ���㣬��ѡ��gene��,"transcript","exon","intron"
  #Э����covariate�Ƕ��������Ӱ��ı������������о����о����Ա�������Ϊʵ���߲���������������Ӱ�졣
  #������adjustvars
  #getFC����ָ����������ʾ����������foldchange(���챶����FC=Exp x1/x2)
  #meas:measurement�������ּ��㷽ʽ
  #mod0=model.matrix(~?, pData(bg_filt)$treat)
  #mod=model.matrix(~?, pData(bg_filt)$treat)
  #results_tanscripts=stattest(bg_filt,feature = "transcript",covariat = ?,adjustvars = ?,mod0 = mod0,mod = mod,meas = "FPKM")���������ظ����˷��������������컯������
}

###ballgown���������ݣ�maxrow > 1;���bg_filt_trans
{
  #����ballgown��
  library(ballgown)
  library(stringr)
  #�鿴��ǰ·��
  getwd()
  #���ô洢·��
  setwd("E:/scientific_research/Transcriptome_analysis/chendan_celegans_RNAseq/mut")
  #���²鿴���Ĺ���ĵ�ǰ·��
  getwd()
  #��ȡ���������ı������ݣ����������ͣ�
  list.files("./ballgown")#�鿴�ļ���
  #�ַ�������
  a = paste(rep("ballgown",length(list.files("./ballgown"))),#��һ������
            list.files("./ballgown"),#�ڶ�������
            sep = "/" #���ӷ�
  )
  bg = ballgown(dataDir = "./ballgown",samples = a )
  #���˵��������͵Ļ���,��ÿһ���������Ļ������ĸ�ֵ��Ϊ0��
  bg_filt=as.data.frame(gexpr(bg))
  bg_filt_trans = subset(bg_filt,apply(bg_filt,1,max)>1,genomesubset=TRUE)
  #gexpr()ѡȡ�Ľ��Ϊ����,texpr()ѡȡ�Ľ��Ϊת¼��
  #feature:���Ǹ�ˮƽ���㣬��ѡ��gene��,"transcript","exon","intron"
  #Э����covariate�Ƕ��������Ӱ��ı������������о����о����Ա�������Ϊʵ���߲���������������Ӱ�졣
  #������adjustvars
  #getFC����ָ����������ʾ����������foldchange(���챶����FC=Exp x1/x2)
  #meas:measurement�������ּ��㷽ʽ
  #mod0=model.matrix(~?, pData(bg_filt)$treat)
  #mod=model.matrix(~?, pData(bg_filt)$treat)
  #results_tanscripts=stattest(bg_filt,feature = "transcript",covariat = ?,adjustvars = ?,mod0 = mod0,mod = mod,meas = "FPKM")���������ظ����˷��������������컯������
}

####��������������mut.txt
{
  write.table(bg_filt_trans,file = "mut.csv",sep = ",",row.names = T) #�б�����Ϊ.csv��ʽ
  write.table(bg_filt_trans,file = "mut.txt",sep = " ",row.names = T) #�б�����Ϊ.txt��ʽ
}

#��������Ϊ�б���ʽ���Ե�һ��Ϊ������
{
  getwd()
  #���ô洢·��
  setwd("E:/scientific_research/Transcriptome_analysis/chendan_celegans_RNAseq/mut")
  #���²鿴���Ĺ���ĵ�ǰ·��
  getwd()
  data = read.table("mut.txt",
                    header = T,
                    row.names = 1,
                    check.names = TRUE#�����������У��Ƿ���Ϲ��򣬲����ϵĻ��Զ��޸�
  )
}

#ȥ���б����б������쳣�Ļ���
{
  data[data==0] = NA
  data = na.omit(data)
}

#���б���������ȡ
{
  #������
  {
    library(reshape)
    #ʹ�øú���ʱ���������������⡣ʹ��dplyr::rename�������˱������Ա����Ʋ���ȡ���ڱ�������ȡ�����ַ�����������ٳ������⣬��������ǰ���λ�ü���
    data = dplyr::rename(data,c("CTRL" = FPKM.4_N2_12_h_Ctr,
                                "MUT_qx666" = FPKM.6_qx666_Day_1,
                                "MUT_epg.5" = FPKM.7_epg.5_Day_1,
                                "MUT_dpy.7" = FPKM.8_dpy.7_Day_1
    ))
    head(data,1)#�鿴����
  }
}

#GO�������������ߣ������ܱ���GO�������˴�ΪMAXrowΪ1���ϵ�,�ļ���ΪRNA in nemateode mutants.csv
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
    cluster_go[cluster_go == NA] = c("Not")
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

#����������
#��Ҫԭʼ����ʱ������˲�
{
  #���������ÿ������ı������Сֵ�����ֵ
  {
    #���������ÿ������ı������Сֵ
    treat = c(colnames(data[-1]))
    ctrl = c(colnames(data[1]))
    data$treat_min = as.matrix(apply(data[,treat],1,min))
    #���������ÿ������ı�������ֵ
    data$treat_max = as.matrix(apply(data[,treat],1,max))
  }
  #ɸѡ����������data_filt
  {
    data$log2FC1 = log2((data$treat_max)/(data$CTRL))
    data$log2FC2 = log2((data$CTRL)/(data$treat_min))
    data_filt = subset(data,data$log2FC1>1|data$log2FC2>1)
    data = data_filt[,c(ctrl,treat)]
  }
}
 
data = log2(data)
 
#����I ��ʮ�ĸ�С����
{
  #����1:4�����������ϲ���С��������,num����
  #order()����������ԭ�򣬣���С����С���������
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
  #ȡlog2
  cluster_data = log2(cluster_data[,1:base::ncol(cluster_data)])
}

#����II�����ֵ���ž��� kmeansa��������cluster_data
#�ֿ�����
data_x = data[,c(1,2)]
data_x = data[,c(1,3)]
data_x = data[,c(1,4)]
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

#���Ͼ���
data_x = data
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
if (FALSE)
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
  genename = c(rownames(cluster_data))
  #�鿴���е�֧�ֺͿ�ת������
  keytypes(org.Hs.eg.db)
  #����������SYMBOL��ʽת��ΪENTREZID��ʽ������org.Hs.eg.dbΪ�������ݿ�
  DEG.entrez_id = bitr(genename,fromType = "SYMBOL",toType = "GOALL",OrgDb = "org.Hs.eg.db")
  ##GO����
  #��������ѧ���ܣ�Moleular Function,MF���û����ڷ��Ӳ���Ĺ�����ʲô�����߻���Щ���ܣ�,����ѧ���̣�Biological Process,BP���û����������Щ����ѧ���̣�,ϸ��ѧ��֣�Cellular Components,CC,������������
  erich.go.ALL=enrichGO(gene = DEG.entrez_id$ENTREZID,
                        OrgDb = org.Mm.eg.db,
                        ont = "ALL",
                        pvalueCutoff = 1,
                        qvalueCutoff = 1,
                        readable = T)
  #�ֱ�ѡ�Ӽ�
  erich.go.MF = erich.go.ALL[which(erich.go.ALL@result$ONTOLOGY =="MF",)]#which�����±�
  erich.go.BP = erich.go.ALL[which(erich.go.ALL@result$ONTOLOGY =="BP",)]
  erich.go.CC = erich.go.ALL[which(erich.go.ALL@result$ONTOLOGY =="CC",)]
  #��GeneRatioΪС��
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
      return("����! ������ȷ�Ĺ��� [category].")
    }
    dat1 = dat[c(1:top.number),c("ID","Description","GeneRatio","pvalue","geneID","Count")]#�������ݵĵ�һ�е���top.number�У�"ID","Description","GeneRadio","pvalue","geneID"��
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

#ALL���ݶ���.csv
{
  write.csv(cluster_data,file = "all high-quality RNA of mut.csv")
}
#Discrepant���ݶ���.csv
{
  write.csv(cluster_data,file = "Discrepant of mut.csv")
}


#��С��ͼ
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
                 cellwidth = 80, cellheight = 0.5,
                 fontsize = 20,
                 
        )
      }
    }
  }
  
}

#�����ͼ
main = c("Differentially expressed RNA in all nemateode mutants")
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
                 cellwidth = 80, cellheight = 0.8,
                 fontsize = 20
        )
      }
    }
  }
  
}




