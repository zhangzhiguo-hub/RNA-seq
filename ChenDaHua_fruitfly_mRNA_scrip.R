###ballgown����������,���bg_filt_trans
{
  #����ballgown��
  library(ballgown)
  #�鿴��ǰ·��
  getwd()
  #���ô洢·��
  setwd("E:/scientific_research/Transcriptome_analysis/fruitfly_maternal")
  #���²鿴���Ĺ���ĵ�ǰ·��
  getwd()
  #��ȡ���������ı������ݣ����������ͣ�
  bg = ballgown(dataDir = "./ballgown", samplePattern = "C", meas = "all")
  #���˵��������͵Ļ���,��ÿһ���������Ļ������ĸ�ֵ��Ϊ0��
  bg_filt=as.data.frame(gexpr(bg))
  bg_filt_row=as.data.frame(bg@expr$trans)
  bg_filt_trans = subset(bg_filt,apply(bg_filt,1,max)>20,genomesubset=TRUE)
  colnames(bg_filt_trans)
  
  #feature:���Ǹ�ˮƽ���㣬��ѡ��gene��,"transcript","exon","intron"
  #Э����covariate�Ƕ��������Ӱ��ı������������о����о����Ա�������Ϊʵ���߲���������������Ӱ�졣
  #������adjustvars
  #getFC����ָ����������ʾ����������foldchange(���챶����FC=Exp x1/x2)
  #meas:measurement�������ּ��㷽ʽ
  #mod0=model.matrix(~?, pData(bg_filt)$treat)
  #mod=model.matrix(~?, pData(bg_filt)$treat)
  #results_tanscripts=stattest(bg_filt,feature = "transcript",covariat = ?,adjustvars = ?,mod0 = mod0,mod = mod,meas = "FPKM")���������ظ����˷��������������컯������
}
#�����б�bg_filt_row���������������򣬶���ԭʼ�б�bg_filt_row_1
{
  library(dplyr)
  bg_filt_row_1= select(bg_filt_row,c(t_name,gene_name,FPKM.C29_L4_340340,FPKM.C29_L4_340340,FPKM.C30_L4_341341,FPKM.C31_L4_342342,FPKM.C32_L4_343343,FPKM.C33_L4_345345,FPKM.C34_L4_346346,FPKM.C35_L4_347347,FPKM.C36_L1_349349,FPKM.C37_L4_350350,FPKM.C38_L4_351351,FPKM.C39_L4_352352,FPKM.C40_L4_353353))
  #������
  {
    library(reshape)
    #ʹ�øú���ʱ���������������⡣ʹ��dplyr::rename�������˱������Ա����Ʋ���ȡ���ڱ�������ȡ�����ַ�����������ٳ������⣬��������ǰ���λ�ü���
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
    head(bg_filt_row_1,1)#�鿴����
  }
  #����������
  {
    #dplyr�����ú������ǿ�������ݿ�ѧtidyverse���ϵĺ��Ĳ���֮һ��ѧϰ����https://bookdown.org/wangminjie/R4DS/colwise.html
    library(dplyr,warn.conflicts = F)
    bg_filt_row_1 = bg_filt_row_1 %>% select("t_name","gene_name","WT0.5-1h_1","WT0.5-1h_2","WT2-3h_1","WT2-3h_2","WT5.5-6h_1","WT5.5-6h_2","MT0.5-1h_1","MT0.5-1h_2","MT2-3h_1","MT2-3h_2","MT5.5-6h_1","MT5.5-6h_2")
  }
}
####��������������
{
  write.table(bg_filt_trans,file = "fruitfly_maternal_bg_filt_trans.csv",sep = ",",row.names = T) #�б�����Ϊ.csv��ʽ
  write.table(bg_filt_trans,file = "fruitfly_maternal_bg_filt_trans.txt",sep = " ",row.names = T) #�б�����Ϊ.txt��ʽ
}
getwd()
#���ô洢·��
setwd("E:/scientific_research/Transcriptome_analysis/fruitfly_maternal")
#���²鿴���Ĺ���ĵ�ǰ·��
getwd()
#��������Ϊ�б���ʽ���Ե�һ��Ϊ������
{
  data = read.table("fruitfly_maternal_bg_filt_trans.txt",
                    header = T,
                    row.names = 1)
}
#ȥ���б����б������쳣�Ļ���
{
  data[data==0] = NA
  data = na.omit(data)
}
#���б���������������
{
  #������
  {
    library(reshape)
    #ʹ�øú���ʱ���������������⡣ʹ��dplyr::rename�������˱������Ա����Ʋ���ȡ���ڱ�������ȡ�����ַ�����������ٳ������⣬��������ǰ���λ�ü���
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
    head(data,1)#�鿴����
  }
  #����
  {
    WT = c("WT0.5-1h_1","WT0.5-1h_2","WT2-3h_1","WT2-3h_2","WT5.5-6h_1","WT5.5-6h_2")
    MT = c("MT0.5-1h_1","MT0.5-1h_2","MT2-3h_1","MT2-3h_2","MT5.5-6h_1","MT5.5-6h_2")
  }
  #����������
  {
    #dplyr�����ú������ǿ�������ݿ�ѧtidyverse���ϵĺ��Ĳ���֮һ��ѧϰ����https://bookdown.org/wangminjie/R4DS/colwise.html
    library(dplyr,warn.conflicts = F)
    data = data %>% select("WT0.5-1h_1","WT0.5-1h_2","WT2-3h_1","WT2-3h_2","WT5.5-6h_1","WT5.5-6h_2","MT0.5-1h_1","MT0.5-1h_2","MT2-3h_1","MT2-3h_2","MT5.5-6h_1","MT5.5-6h_2")
  }
}
#ȡ��ֵ�����newdata
#�ɹ�#��Ҫ�ȶ��ظ���ʱ��ִ��
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
#WTԭʼ����,���data_wt_cluster����ȡlog2������Ϊdata(newdata)
data = newdata#���ݾ�ֵ����
WT = colnames(newdata)[c(1,2,3)]
{
  data_wt = data[,which(colnames(data) %in% WT)]
  #data_wt_cluster
  i = 1
  data_wt_cluster = data.frame()
  while(i <= length(data_wt))
  {
    a =  subset(data_wt,data_wt[,i] > apply(data_wt[,colnames(data_wt)[-i]],1,max))
    #���а���colnames(data_wt)[i]��ֵ��������
    a = a[order(-a[i]),]
    data_wt_cluster = rbind(data_wt_cluster,a)
    i <- i + 1
  }
  #��ȡlog2ֵ
  data_wt_cluster = log2(data_wt_cluster[,1:base::ncol(data_wt_cluster)])
}
#WT��MT�����࣬ɸѡ���data_cluster����ȡlog2
data = newdata#���ݾ�ֵ����
{
  data_cluster = data[rownames(data_wt_cluster),]
  #ȡlog2
  data_cluster = log2(data_cluster)
}

#ɸѡ������ ͻ���͡����ֵ������һ�еĲС�� Ұ���͡����ֵ������һ�еĲ��list
#Maternal_mRNA :WT(0.5~1h)-WT(5.5~6h) > MUT(0.5~1h)-MUT(5.5~6h)##���data_cluster1,���ȡlog2
{
  data_cluster1 = subset(data_wt,data_wt[,1]>20)
  data_cluster1 = subset(data_wt,data_wt[,1] > data_wt[,2] & data_wt[,2] >= data_wt[,3])
  data_cluster1 = newdata[rownames(data_cluster1),]
  data_cluster1 = subset(data_cluster1,(data_cluster1[,1] - data_cluster1[,3]) > (data_cluster1[,4] - data_cluster1[,6]),genomesubset=TRUE)
  data_cluster1 = log2(data_cluster1)
  a = tibble::as_tibble(bg_filt_row_1[which(bg_filt_row_1$gene_name %in% rownames(data_cluster1)),])
  #ȥ���б����б������쳣�Ļ���
  {
    a[a==0] = NA
    a = na.omit(a)
  }
  a$gene_name = factor(a$gene_name)
  a = a %>% group_by(gene_name) %>% top_n(1,`WT0.5-1h_1`)
  write.csv(a,file = paste(getwd(),c("Maternal_RNA.csv"),sep = "/"))
}
#Zygote_mRNA ##δ֪##���data_cluster3
{
  data_cluster3 = subset(data_wt,apply(data_wt[,colnames(data_wt)[-1]],1,max)>20)
  data_cluster3 = subset(data_cluster3,apply(data_cluster3[,colnames(data_cluster3)[-1]],1,min)/data_cluster3[,1]>1.5)
  data_cluster3 = newdata[rownames(data_cluster3),]
  a = tibble::as_tibble(bg_filt_row_1[which(bg_filt_row_1$gene_name %in% rownames(data_cluster3)),])
  #ȥ���б����б������쳣�Ļ���
  {
    a[a==0] = NA
    a = na.omit(a)
  }
  a$gene_name = factor(a$gene_name)
  a = a %>% group_by(gene_name) %>% top_n(1,`WT5.5-6h_1`)
  write.csv(a,file = paste(getwd(),c("Zygote.csv"),sep = "/"))
}

#����ɸѡ�Ļ��򣬷������ε�bedgraph�ļ����鿴��ͬת¼����3'5'�˵ķֲ�
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
    #ȥ���б����б�����Ϊ0��ת¼��
    {
      A[A==0] = NA
      A = na.omit(A)
    }
    A = data.frame(tapply(A$V3,A$V1, median)) %>% (function(x){cbind(x,V1 = rownames(x))}) %>% (function(x){merge(A,x)})
    colnames(A)[length(A)]="Median"
    A$V3 = (A$V3/A$Median)
    A = left_join(A,reference,by = "V1",copy = T)
    #��bin
    A$site = (A$V2.x/A$V2.y)*100
    A$bin = factor(ceiling(A$site))
    #counts������ת¼������,Ϊ����bin�н��б�׼��
    names(A)[3]=c("The number of reads")
    A$`The number of reads`=A$`The number of reads`/A$V2.y
    A = data.frame(aggregate(A$`The number of reads`~A$bin, A,sum)) %>% (function(x){cbind(x,site=rownames(x))}) 
    A$group = substr(i,1,nchar(i)-3)
    A$Type1=ifelse(unlist(strsplit(i,"_"))[2]=="Zygote","Zygote","Maternal")
    A$Type2=ifelse(substr(i,1,2)=="WT","WT","MT")
    a = rbind(a,A)
  }
  #��ͼ
  library(ggplot2)
    A$A.bin=factor(A$A.bin)
    title="mRNA"
    ggplot(data=a,mapping=aes(x=A.bin,y=A..The.number.of.reads.,group=group,colour = group,linetype=Type1,shape=Type2))+
    geom_line(size = 0.5)+
    geom_point(size = 1.24)+
    labs(x = "Bins",y = "The number of reads",title = "mRNA")+
    theme(panel.grid = element_blank())+#ȥ��������
    scale_x_discrete(breaks=seq(0,100,10))+#�޸�x��Ŀ̶ȼ��,ʹ��discreteʱ���б���������
    scale_y_discrete(breaks=seq(0,15,3))
}   

#�´��;���,���cluster_data,��Ҫ����720*6678�Ρ�������ʧ��
data = newdata
{
  a = c(1:6)
  #����1:6�����������ϲ���С��������,num����
  #order()����������ԭ�򣬣���С����С���������
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
  #ȡlog2
  cluster_data = log2(cluster_data[,1:base::ncol(cluster_data)])
}

#δ���
#GO�������������ߣ������ܱ���GO�������˴�ΪMAXrowΪ1���ϵ�,�ļ���ΪRNA in nemateode mutants.csv#DAVID�޷�ʶ�𳬹�3000��ID��
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
  }
  #-----
}

#�ֿ�GO
{
  #data_cluster1,GO�������������ߣ�hclust���࣬����б�Ϊ��cluster1_go,����ļ���ΪWT(0.5~1h)_MAX of fruitfly.csv
  {
    write.table(data.frame(geneid = row.names(data_cluster1)),file = "data_cluster1.txt",row.names = F,col.names = F,quote = F)#����������Ϊ���������������Ϊ.txt��ʽ����Ҫ��������Ҫ����,��Ҫ����
    #-----
    #�����ݻ���idΪsympol
    #��������DAVID  https://david.ncifcrf.gov/,����ճ���ύ��Identifierѡ��WORMBASE_GENE_ID.������ѡ�����е� Gene ID Conversion Tool��ת��OFFICAL_GENE_SYMBOL������ͨ�����ƣ����ؼ��ɣ�����Ϊall.txt
    #���벢������
    {
      GO_gene_name = read.table("data_cluster1.txt",header = T,sep ="\t",#Tap�ָ�
                                fill=TRUE,#�ò����ɺ������ݲ�����Ĵ���ֱ�Ӷ�ȡ
                                stringsAsFactors = FALSE,quote = ""#��Ҫ����
      )
      GO_gene_name = GO_gene_name[which(GO_gene_name[,"Species"]=="Drosophila melanogaster"),]
      cluster1 = data_cluster1#log2�Ѿ�ȡ��
      #����
      dist_method = "euclidean"
      #����������ķ�������ƽ������average;���ķ���centroid;�м���뷨��median;����뷨��complete(Ĭ��);��̾��뷨��single;���ƽ���ͷ���ward;�ܶȹ��Ʒ���density
      cluster_method = "complete"
      b = dist(cluster1,method = dist_method)
      b = hclust(b,method = cluster_method)
      order = b$order
      cluster1 = cluster1[order,]
      cluster1$SYMBOL_ID = c(SYMBOL_ID  = row.names(cluster1))
      #������
      library(reshape)
      GO_gene_name = dplyr::rename(GO_gene_name,c("SYMBOL_ID" = From,
                                                  "NCBI_ID" = To,
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
      cluster1_go = merge(cluster1,GO_gene_name,sort = F,all = T)
      #rownames(cluster1_go) = cluster1_go[,"SYMBOL_ID"]
      #����
      #������
      library(dplyr,warn.conflicts = F)
      colnames(cluster1_go)
      cluster1_go = cluster1_go %>% select("SYMBOL_ID","NCBI_ID","Species","Discription","WT0.5-1h","WT2-3h","WT5.5-6h","MT0.5-1h","MT2-3h","MT5.5-6h")
    }
    #�����б���ע��
    {
      write.csv(cluster1_go,file = "WT(0.5~1h)_MAX of fruitfly.csv")
    }
    #-----
  }
  #data_cluster2,GO�������������ߣ�hclust���࣬����б�Ϊ��cluster2_go,����ļ���ΪWT(2~3h)_MAX of fruitfly.csv
  {
    write.table(data.frame(geneid = row.names(data_cluster2)),file = "data_cluster2.txt",row.names = F,col.names = F,quote = F)#����������Ϊ���������������Ϊ.txt��ʽ����Ҫ��������Ҫ����,��Ҫ����
    #-----
    #�����ݻ���idΪsympol
    #��������DAVID  https://david.ncifcrf.gov/,����ճ���ύ��Identifierѡ��WORMBASE_GENE_ID.������ѡ�����е� Gene ID Conversion Tool��ת��OFFICAL_GENE_SYMBOL������ͨ�����ƣ����ؼ��ɣ�����Ϊall.txt
    #���벢������
    {
      GO_gene_name = read.table("data_cluster2.txt",header = T,sep ="\t",#Tap�ָ�
                                fill=TRUE,#�ò����ɺ������ݲ�����Ĵ���ֱ�Ӷ�ȡ
                                stringsAsFactors = FALSE,quote = ""#��Ҫ����
      )
      GO_gene_name = GO_gene_name[which(GO_gene_name[,"Species"]=="Drosophila melanogaster"),]
      cluster2 = data_cluster2#log2�Ѿ�ȡ��
      #����
      dist_method = "euclidean"
      #����������ķ�������ƽ������average;���ķ���centroid;�м���뷨��median;����뷨��complete(Ĭ��);��̾��뷨��single;���ƽ���ͷ���ward;�ܶȹ��Ʒ���density
      cluster_method = "complete"
      b = dist(cluster2,method = dist_method)
      b = hclust(b,method = cluster_method)
      order = b$order
      cluster2 = cluster2[order,]
      cluster2$SYMBOL_ID = c(SYMBOL_ID  = row.names(cluster2))
      #������
      library(reshape)
      GO_gene_name = dplyr::rename(GO_gene_name,c("SYMBOL_ID" = From,
                                                  "NCBI_ID" = To,
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
      cluster2_go = merge(cluster2,GO_gene_name,sort = F,all = T)
      #rownames(cluster1_go) = cluster1_go[,"SYMBOL_ID"]
      #����
      #������
      library(dplyr,warn.conflicts = F)
      colnames(cluster2_go)
      cluster2_go = cluster2_go %>% select("SYMBOL_ID","NCBI_ID","Species","Discription","WT0.5-1h","WT2-3h","WT5.5-6h","MT0.5-1h","MT2-3h","MT5.5-6h")
    }
    #�����б���ע��
    {
      write.csv(cluster2_go,file = "WT(2~3h)_MAX of fruitfly.csv")
    }
    #-----
  }
}


#�����ͼ
{
  library(pheatmap)
  gaps_row = c()
  sum = 0
  Cluster = c()
  length_Cluster = c()
  i = 1
  while (i <= length(colnames(data_wt_cluster))) #���Ʋ�ͬ��ͼ�ǵø�����
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
  #��data_cluster�����ͼ
  {
    #����ע��
    annotation_col= data.frame(Time = c(0.75,2.5,5.75,0.75,2.5,5.75))
    annotation_row = data.frame(Cluster = rep(Cluster,length_Cluster))
    color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")#��ɫ
    #����û��Ӧʱ������dev.off()
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
                 #�ָ�
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
#����С��ͼ
data_cluster=data_cluster1
main = "WT(0.5~1h)_MAX of fruitfly:
WT(0.5~1h)-WT(2~3h) > MUT(0.5~1h)-MUT(2~3h)"
data_cluster=data_cluster2
main = "WT(2~3h)_MAX of fruitfly:
WT(2~3h)-WT(5.5~6h) > MUT(2~3h)-MUT(5.5~6h)"
{
  color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")#��ɫ
  #����û��Ӧʱ������dev.off()
  {
    #data_cluster
    {
      library(pheatmap)
      pheatmap(data_cluster,
               color = colorRampPalette(color.key)(50),
               border_color = NA,
               cluster_cols = F,
               cluster_rows = T,
               #�ָ�
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



##����ͼ
{
  library(pheatmap)
  #��dataת��Ϊmatrix��ʽ
  #�ֱ�AB�����ͼ
  {
    #����ע��
    annotation_col = data.frame(Grop= rep(c("Control","Treat"),c(1,5)),Time = c(0,1:5))
    color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")#��ɫ
    #����û��Ӧʱ������dev.off()
    
    #��A��Ļ�����ͼ
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
                     #�ָ�
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
    
    #��B��Ļ�����ͼ
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
                     #�ָ�
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
  #��AB���鹲�е�ͼ
  {
    #����ע��
    annotation_col_ab = data.frame(Grop= rep(c("Control","TreatA","TreatB"),c(1,5,5)),Dose = c(0,1:5,1:5))
    annotation_row_ab = data.frame(Cluster = rep(c("Max_CTRL_33(7%)","Max_A15_259(53%)","Max_A30_23(5%)","Max_A60_17(3%)","Max_A90_17(3%)","Max_A120_52(10%)","Max_B15_18(4%)","Max_B30_19(4%)","Max_B60_22(5%)","Max_B90_17(3%)","Max_B120_16(3%)"),c(33,259,23,17,17,52,18,19,22,17,16)))
    color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")#��ɫ
    #����û��Ӧʱ������dev.off()
    #��ab�鹲�еĻ�����ͼ
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
                      #�ָ�
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
##��vennͼ
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
                   #alpha = 0.1,#����͸����
                   filename = NULL,#imagetype="tiff",
                   category = c("A Treat", "B Treat"),
                   lwd=1 ,#����Բ������
                   #lty=2 ,#����Բ������
                   col=c('cornflowerblue','yellow'),
                   fill = c("cornflowerblue", "yellow"),
                   cat.col = c("black", "black"),#�������Ƶ���ʾ��ɫ
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
##UpSetR����ͼ
{
  #1��������ʽ����R������������ݿ��ˡ��б�ʾԪ�أ��б�ʾ���ݼ�����Ͷ�����Ϣ��
  #2��Ԫ�����ļ���(û��������֪������)fromList
  #3��venneuler������������������Ͻ���������fromExpression
  library(UpSetR)
  #���庯��
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
  #��������
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
  #ת��
  data_2 <- fromExpression(input)
  upset(data_2,nsets = 10,nintersects = 26,
        text.scale = c(1.3, 1.3, 1, 1, 1.5, 1)) #�������֣��ֱ����c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars))
  
  
}

#ƴ��
if(FALSE)#ȫ��ע��
{
  #��װ�����ذ�
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
  #ת����ͼΪggplot��ʽ
  {
    library(ggplotify)
    a = as.ggplot(a)
    b = as.ggplot(b)
    ab = as.ggplot(ab)
    v = as.ggplot(v)
  }
  #������
  {
    grid.newpage() # �½�����
    layout_1 <- grid.layout(nrow = 2, ncol = 2,widths = c(10.7,80), 
                            heights = c(base::nrow(data_filt_a_cluster1) * 0.5,
                                        base::nrow(data_filt_a_cluster2) * 0.5,
                                        base::nrow(data_filt_a_cluster3) * 0.5,
                                        base::nrow(data_filt_a_cluster4) * 0.5,
                                        base::nrow(data_filt_a_cluster5) * 0.5,
                                        base::nrow(data_filt_a_cluster6) * 0.5
                            )
    ) # �ֳ�4��1�й�4�����
    pushViewport(viewport(layout = layout_1)) # �Ƴ���Ϊ4�������Ӵ�
  }
  #��ӡͼ
  {
    print(a1, vp = viewport(layout.pos.row = 1,
                            layout.pos.col = 2
    )) # ��a1�������һ��
    print(a2, vp = viewport(layout.pos.row = 2,
                            layout.pos.col = c(1,2)
    )) # ��a2������ڶ���
    print(a3, vp = viewport(layout.pos.row = 3,
                            layout.pos.col = c(1,2),
                            just = bottom
    )) # ��a3�����������
    print(a4, vp = viewport(layout.pos.row = 4,
                            layout.pos.col = c(1,2)
    )) # ��a4�����������
    print(a5, vp = viewport(layout.pos.row = 5,
                            layout.pos.col = c(1,2)
    )) # ��a5�����������
    print(a6, vp = viewport(layout.pos.row = 6,
                            layout.pos.col = c(1,2)
    )) # ��a6�����������
    # ����������ͷֱ���
  }
}
#�ӱ�ǩ
if(FALSE)#ȫ��ע��
{
  #ת����ͼΪggblot��ʽ
  ��
  grid.text(label="        Clucster1 78(17%)
          Low to high", x = 0.39,y = 0.85,gp=gpar(col="black",fontsize = 12,draw=TRUE))
  grid.text(label="        Clucster2 211(46%)
          Low to high to lower", x = 0.39,y =0.6,gp=gpar(col="black",fontsize = 12,draw=TRUE,just = c("left", "top")))
  grid.text(label="        Clucster3 44(10%)
          high to low", x = 0.39,y =0.34,gp=gpar(col="black",fontsize = 12,draw=TRUE))
  grid.text(label="        Clucster4 122(27%)
          high to low to higher", x = 0.39,y =0.16,gp=gpar(col="black",fontsize = 12,draw=TRUE)) 
}

#dev.off()#�رյ�ǰ�豸


####GOע��
#��������
data_cluster = read.csv("WT(0.5~1h)_MAX of fruitfly.csv",header = T,row.names = 3)
{
  #BiocManager::install("clusterProfiler")
  #BiocManager::install("topGO")
  #BiocManager::install("Rgraphviz")
  #BiocManager::install("pathview")
  #BiocManager::install("org.Dm.eg.db"),��Ӭ�İ�
  ##���ػ�ͼ��
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
  #�鿴���е�֧�ֺͿ�ת������
  keytypes(org.Dm.eg.db)
  #����������SYMBOL��ʽת��ΪENTREZID��ʽ������org.Dm.eg.dbΪfly���ݿ�
  DEG.entrez_id = data.frame()
  DEG.entrez_id = bitr(genename,fromType = "ENSEMBLTRANS",toType = "FLYBASE",OrgDb = "org.Dm.eg.db")
  ##GO����
  #��������ѧ���ܣ�Moleular Function,MF���û����ڷ��Ӳ���Ĺ�����ʲô�����߻���Щ���ܣ�,����ѧ���̣�Biological Process,BP���û����������Щ����ѧ���̣�,ϸ��ѧ��֣�Cellular Components,CC,������������
  erich.go.ALL=enrichGO(gene = DEG.entrez_id$ENTREZID,
                        OrgDb = org.Hs.eg.db,
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
  
}
# Read in data and generate the plots
# Biological Process
DrawGOBubblePlot(erich.go.BP,"BP",10,"blue")
DrawGOBubblePlot(erich.go.CC,"CC",10,"blue")
DrawGOBubblePlot(erich.go.MF,"MF",10,"blue")

