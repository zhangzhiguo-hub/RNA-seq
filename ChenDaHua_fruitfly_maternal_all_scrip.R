###ballgown����������,������Ϊbg_filt_trans��maxrowΪ>20
{
  #����ballgown��
  library(ballgown)
  library(stringr)
  #�鿴��ǰ·��
  getwd()
  #���ô洢·��
  setwd("E:/scientific_research/Transcriptome_analysis/fruitfly_maternal_allrna")
  #���²鿴���Ĺ���ĵ�ǰ·��
  getwd()
  #��ȡ���������ı������ݣ����������ͣ�
  list.files("./ballgown2")#�鿴�ļ���
  #�ַ�������
  a = paste(rep("ballgown2",length(list.files("./ballgown2"))),#��һ������
            list.files("./ballgown2"),#�ڶ�������
            sep = "/" #���ӷ�
            )
  bg = ballgown(dataDir = "./ballgown2",samples = a )
  #���˵��������͵Ļ���,��ÿһ���������Ļ������ĸ�ֵ��Ϊ0��
  bg_filt=as.data.frame(gexpr(bg))
  bg_filt_trans = subset(bg_filt,apply(bg_filt,1,max)>20,genomesubset=TRUE)
  colnames(bg_filt_trans)
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

####��������������
{
  write.table(bg_filt_trans,file = "fruitfly_maternal_allrna_bg_filt_trans.csv",sep = ",",row.names = T) #�б�����Ϊ.csv��ʽ
  write.table(bg_filt_trans,file = "fruitfly_maternal_allrna_bg_filt_trans.txt",sep = " ",row.names = T) #�б�����Ϊ.txt��ʽ
}
getwd()
#���ô洢·��
setwd("E:/scientific_research/Transcriptome_analysis/fruitfly_maternal_allrna")
#���²鿴���Ĺ���ĵ�ǰ·��
getwd()

#��������Ϊ�б���ʽ���Ե�һ��Ϊ������
{
  data = read.table("fruitfly_maternal_allrna_bg_filt_trans.txt",
                    header = T,
                    row.names = 1,
                    check.names = T#�����������У��Ƿ���Ϲ��򣬲����ϵĻ��Զ��޸�
                    )
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
    data = dplyr::rename(data,c("WT0_1h_1" = FPKM.1.2N,
                                "WT2_3h_1" = FPKM.2.2N,
                                "WT5_6h_1" = FPKM.3.2N,
                                "WT0_1h_2" = FPKM.A.0.1,
                                "WT0_1h_3" = FPKM.A.0.2,
                                "WT2_3h_2" = FPKM.A.2.1,
                                "WT2_3h_3" = FPKM.A.2.2,
                                "WT5_6h_2" = FPKM.A.5.1,
                                "WT5_6h_3" = FPKM.A.5.2
                                ))
    head(data,1)#�鿴����
  }
  #����������
  {
    #dplyr�����ú������ǿ�������ݿ�ѧtidyverse���ϵĺ��Ĳ���֮һ��ѧϰ����https://bookdown.org/wangminjie/R4DS/colwise.html
    library(dplyr,warn.conflicts = F)
    data = data %>% select("WT0_1h_1","WT0_1h_2","WT0_1h_3","WT2_3h_1","WT2_3h_2","WT2_3h_3","WT5_6h_1","WT5_6h_2","WT5_6h_3")
  }
  }
  
#ȡ��ֵ�����newdata
#�ɹ�#��Ҫ�ȶ��ظ���ʱ��ִ��
{
  i = 1
  newdata = data
  while (i < length(colnames(data))){
    newdata$A <- rowMeans(data[,c(i,i+1,i+2)]) 
    a = c(substr(colnames(data[i]),1,nchar(colnames(data[i]))-2))
    colnames(newdata)[length(newdata)] = a
    i = i+3
  }
  newdata = newdata[,-which(colnames(newdata) %in% colnames(data))]
}

#�����ͻ���ͱȽϣ�ɸ��ͻ���͵�genelist,�鿴ʣ�µķǱ���RNA
{
  a = read.table("E:/scientific_research/Transcriptome_analysis/fruitfly_maternal/fruitfly_maternal_bg_filt_trans.txt",row.names = 1,header = T)
  a[a==0] = NA
  a = na.omit(a)
  a = intersect(rownames(a),rownames(data))
  data_cluster = newdata[!rownames(newdata) %in% a,]
  #����
  write.csv(data_cluster,file = "non-coding RNAs.csv",row.names = T)
  data_cluster = log2(data_cluster)
  data = data_cluster
}

#ԭʼ����
{
    #attach(data)
    #cluster
  i <- 1
  data_cluster = data.frame()
  while(i <= length(data))
    {
      a =  subset(data,data[,i] > apply(data[,colnames(data)[-i]],1,max))
      #���а���colnames(data)[1]��ֵ��������
      a = a[order(-a[i]),]
      data_cluster = rbind(data_cluster,a)
      i <- i + 1
    }
    #��ȡlog2ֵ
    data_cluster = log2(data_cluster[,1:base::ncol(data_cluster)])
}


#�¾���?δ��
{
  #����1:4�����������ϲ���С��������,num����
  #order()����������ԭ�򣬣���С����С���������
  {
    num = c() 
    b=c()
    a = (1:4)
    for (i in a) {
      b[1]=a[i]
      for (j in a[-i]) {
        b[2]=a[j]
        for (k in a[c(-i,-j)]) {
          b[3]=a[k]
          for (l in a[c(-i,-j,-k)]) {
            b[4]=a[l]
            num = c(num,paste(b,collapse = ""))
          }
        }
      }
    }
    num = num[c(grep("...1",num),grep("...2",num),grep("...3",num),grep("...4",num))]
  }
  
            n=n+1
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


#����ͼ������data_cluster
main = "non-coding RNAs"
{
  library(pheatmap)
  gaps_row = c()
  sum = 0
  Cluster = c()
  length_Cluster = c()
  i = 1
  while (i <= length(colnames(data_cluster))) #���Ʋ�ͬ��ͼ�ǵø�����
  {
    a = nrow(subset(data,data[,i] > apply(data[,colnames(data)[-i]],1,max)))
    b = paste (colnames(data)[i],
               a,
               sep = " ")
    b = paste (b,
               round(a/nrow(data)*100),
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
  #��WTallRNA�����ͼ
  {
    #����ע��
    annotation_col= data.frame(Time = c(0.5,2.5,5.5))
    annotation_row = data.frame(Cluster = rep(Cluster,length_Cluster))
    color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")#��ɫ
    #����û��Ӧʱ������dev.off()
    {
      rownames(annotation_col)=colnames(data_cluster)
      rownames(annotation_row)=rownames(data_cluster)
      #data_cluster
      {
        ab = pheatmap(data_cluster,
                      color = colorRampPalette(color.key)(50),
                      border_color = NA,
                      cluster_cols = F,
                      cluster_rows = F,
                      #�ָ�
                      gaps_row = gaps_row,
                      annotation_col = annotation_col,
                      annotation_row = annotation_row,
                      main = main,
                      #scale = "row",
                      show_rownames = F,show_colnames = T,
                      cellwidth = 20, cellheight = 10,
                      fontsize = 20
        )
      }
    }
  }
  
}







#ɸѡMT,MT�����������в������Ļ���
{
  #���a,b����������ÿ������ı������Сֵ�����ֵ
  {
    #���a,b����������ÿ������ı������Сֵ
    data$min_a = as.matrix(apply(data[,treat_a],1,min))
    data$min_b = as.matrix(apply(data[,treat_b],1,min))
    #���a,b����������ÿ������ı�������ֵ
    data$max_a = as.matrix(apply(data[,treat_a],1,max))
    data$max_b = as.matrix(apply(data[,treat_b],1,max))
  }
  #ɸѡab�����������в������Ļ���
  {
    #ɸѡA��Ĳ���������
    {
      data$log2FC_a_con_1 = log2((data$max_a)/(data$CTRL))
      data$log2FC_a_con_2 = log2((data$CTRL)/(data$min_a))
      data_filt = subset(data,data$log2FC_a_con_1>2|data$log2FC_a_con_2>2)
      data_filt_a = data_filt[,c(control,treat_a)]
      ##A�����
      {
        attach(data_filt_a)
        #cluster1-max_CTRL
        {
          data_filt_a_cluster1 =  subset(data_filt_a,CTRL>apply(data_filt_a[,treat_a],1,max))
          #��control��0���ο�
          #data_filt_a_cluster1 = data_filt_a_cluster1 - data_filt_a_cluster1 $ CTRL
          #���а���CTRL��ֵ��������
          {
            data_filt_a_cluster1 = data_filt_a_cluster1[order(-data_filt_a_cluster1$CTRL),]
          }
        }
        #cluster2-max_A15
        {
          data_filt_a_cluster2 =  subset(data_filt_a,A15>apply(data_filt_a[,c("CTRL","A30","A60","A90","A120")],1,max))
          #��control��0���ο�
          #data_filt_a_cluster2 = data_filt_a_cluster2 - data_filt_a_cluster2 $ CTRL
          #���а���CTRL��ֵ��������
          {
            data_filt_a_cluster2 = data_filt_a_cluster2[order(-data_filt_a_cluster2 $ A15),]
          }
        }
        #cluster3-max_A30
        {
          data_filt_a_cluster3 = subset(data_filt_a,A30>apply(data_filt_a[,c("CTRL","A15","A60","A90","A120")],1,max))
          #��control��0���ο�
          #data_filt_a_cluster2 = data_filt_a_cluster2 - data_filt_a_cluster2 $ CTRL
          #���а���CTRL��ֵ��������
          {
            data_filt_a_cluster3 = data_filt_a_cluster3[order(-data_filt_a_cluster3 $ A30),]
          }
        }
        #cluster4-max_A60
        {
          data_filt_a_cluster4 =  subset(data_filt_a,A60>apply(data_filt_a[,c("CTRL","A30","A15","A90","A120")],1,max))
          #��control��0���ο�
          #data_filt_a_cluster4 = data_filt_a_cluster4 - data_filt_a_cluster4 $ CTRL
          #���а���CTRL��ֵ��������
          {
            data_filt_a_cluster4 = data_filt_a_cluster4[order(-data_filt_a_cluster4 $ A60),]
          }
        }
        #cluster5-max_A90
        {
          data_filt_a_cluster5 =  subset(data_filt_a,A90>apply(data_filt_a[,c("CTRL","A30","A15","A60","A120")],1,max))
          #��control��0���ο�
          #data_filt_a_cluster4 = data_filt_a_cluster4 - data_filt_a_cluster4 $ CTRL
          #���а���A90��ֵ��������
          {
            data_filt_a_cluster5 = data_filt_a_cluster5[order(-data_filt_a_cluster5 $ A90),]
          }
        }
        #cluster6-max_A120
        {
          data_filt_a_cluster6 =  subset(data_filt_a,A120>apply(data_filt_a[,c("CTRL","A30","A15","A90","A60")],1,max))
          #��control��0���ο�
          #data_filt_a_cluster6 = data_filt_a_cluster6 - data_filt_a_cluster6 $ CTRL
          #���а���A120��ֵ��������
          {
            data_filt_a_cluster6 = data_filt_a_cluster6[order(-data_filt_a_cluster6 $ A120),]
          }
        }
        detach(data_filt_a)
        #�ϲ�
        data_filt_a_cluster = rbind(data_filt_a_cluster1,data_filt_a_cluster2,data_filt_a_cluster3,data_filt_a_cluster4,data_filt_a_cluster5,data_filt_a_cluster6)
        #��ȡlog2ֵ
        data_filt_a_cluster = log2(data_filt_a_cluster[,1:base::ncol(data_filt_a_cluster)])
      }
    }
    #ɸѡB��Ĳ���������
    {
      data$log2FC_b_con_1 = log2((data$max_b)/(data$CTRL))
      data$log2FC_b_con_2 = log2((data$CTRL)/(data$min_b))
      data_filt = subset(data,data$log2FC_b_con_1>2|data$log2FC_b_con_2>2)
      data_filt_b = data_filt[,c(control,treat_b)]
      ##B�����
      {
        attach(data_filt_b)
        #cluster1-max_CTRL
        {
          data_filt_b_cluster1 =  subset(data_filt_b,CTRL>apply(data_filt_b[,treat_b],1,max))
          #��control��0���ο�
          #data_filt_b_cluster1 = data_filt_b_cluster1 - data_filt_b_cluster1 $ CTRL
          #���а���CTRL��ֵ��������
          {
            data_filt_b_cluster1 = data_filt_b_cluster1[order(-data_filt_b_cluster1$CTRL),]
          }
        }
        #cluster2-max_B15
        {
          data_filt_b_cluster2 =  subset(data_filt_b,B15>apply(data_filt_b[,c("CTRL","B30","B60","B90","B120")],1,max))
          #��control��0���ο�
          #data_filt_b_cluster2 = data_filt_b_cluster2 - data_filt_b_cluster2 $ CTRL
          #���а���B15��ֵ��������
          {
            data_filt_b_cluster2 = data_filt_b_cluster2[order(-data_filt_b_cluster2 $ B15),]
          }
        }
        #cluster3-max_B30
        {
          data_filt_b_cluster3 = subset(data_filt_b,B30>apply(data_filt_b[,c("CTRL","B15","B60","B90","B120")],1,max))
          #��control��0���ο�
          #data_filt_b_cluster2 = data_filt_b_cluster2 - data_filt_b_cluster2 $ CTRL
          #���а���B30��ֵ��������
          {
            data_filt_b_cluster3 = data_filt_b_cluster3[order(-data_filt_b_cluster3 $ B30),]
          }
        }
        #cluster4-max_A60
        {
          data_filt_b_cluster4 =  subset(data_filt_b,B60>apply(data_filt_b[,c("CTRL","B30","B15","B90","B120")],1,max))
          #��control��0���ο�
          #data_filt_b_cluster4 = data_filt_b_cluster4 - data_filt_b_cluster4 $ CTRL
          #���а���B60��ֵ��������
          {
            data_filt_b_cluster4 = data_filt_b_cluster4[order(-data_filt_b_cluster4 $ B60),]
          }
        }
        #cluster5-max_A90
        {
          data_filt_b_cluster5 =  subset(data_filt_b,B90>apply(data_filt_b[,c("CTRL","B30","B15","B60","B120")],1,max))
          #��control��0���ο�
          #data_filt_b_cluster4 = data_filt_b_cluster4 - data_filt_b_cluster4 $ CTRL
          #���а���A90��ֵ��������
          {
            data_filt_b_cluster5 = data_filt_b_cluster5[order(-data_filt_b_cluster5 $ B90),]
          }
        }
        #cluster6-max_A120
        {
          data_filt_b_cluster6 =  subset(data_filt_b,B120>apply(data_filt_b[,c("CTRL","B30","B15","B90","B60")],1,max))
          #��control��0���ο�
          #data_filt_b_cluster6 = data_filt_b_cluster6 - data_filt_b_cluster6 $ CTRL
          #���а���B120��ֵ��������
          {
            data_filt_b_cluster6 = data_filt_b_cluster6[order(-data_filt_b_cluster6 $ B120),]
          }
        }
        detach(data_filt_b)
        #�ϲ�
        data_filt_b_cluster = rbind(data_filt_b_cluster1,data_filt_b_cluster2,data_filt_b_cluster3,data_filt_b_cluster4,data_filt_b_cluster5,data_filt_b_cluster6)
        #��ȡlog2ֵ
        data_filt_b_cluster = log2(data_filt_b_cluster[,1:base::ncol(data_filt_b_cluster)])
      }
    }
    #ɸѡAB������Ĳ���������
    {
      data_filt = subset(data,(data$log2FC_b_con_1>2|data$log2FC_b_con_2>2)
                         &
                           (data$log2FC_a_con_1>2|data$log2FC_a_con_2>2)
      )
      data_filt_ab = data_filt[,c(control,treat_a,treat_b)]
      ##��A���������
      {
        attach(data_filt_ab)
        #cluster1-max_CTRL
        {
          data_filt_ab_cluster1 =  subset(data_filt_ab,CTRL>apply(data_filt_ab[,c(treat_a,treat_b)],1,max))
          #��control��0���ο�
          #data_filt_ab_cluster1 = data_filt_ab_cluster1 - data_filt_ab_cluster1 $ CTRL
          #���а���CTRL��ֵ��������
          {
            data_filt_ab_cluster1 = data_filt_ab_cluster1[order(-data_filt_ab_cluster1$CTRL),]
          }
        }
        #cluster2-max_A15
        {
          data_filt_ab_cluster2 =  subset(data_filt_ab,A15>apply(data_filt_ab[,c("CTRL","A30","A60","A90","A120",treat_b)],1,max))
          #��control��0���ο�
          #data_filt_ab_cluster2 = data_filt_ab_cluster2 - data_filt_ab_cluster2 $ CTRL
          #���а���A15��ֵ��������
          {
            data_filt_ab_cluster2 = data_filt_ab_cluster2[order(-data_filt_ab_cluster2 $ A15),]
          }
        }
        #cluster3-max_A30
        {
          data_filt_ab_cluster3 = subset(data_filt_ab,A30>apply(data_filt_ab[,c("CTRL","A15","A60","A90","A120",treat_b)],1,max))
          #��control��0���ο�
          #data_filt_a_cluster2 = data_filt_a_cluster2 - data_filt_a_cluster2 $ CTRL
          #���а���A30��ֵ��������
          {
            data_filt_ab_cluster3 = data_filt_ab_cluster3[order(-data_filt_ab_cluster3 $ A30),]
          }
        }
        #cluster4-max_A60
        {
          data_filt_ab_cluster4 =  subset(data_filt_ab,A60>apply(data_filt_ab[,c("CTRL","A30","A15","A90","A120",treat_b)],1,max))
          #��control��0���ο�
          #data_filt_a_cluster4 = data_filt_a_cluster4 - data_filt_a_cluster4 $ CTRL
          #���а�A60��ֵ��������
          {
            data_filt_ab_cluster4 = data_filt_ab_cluster4[order(-data_filt_ab_cluster4 $ A60),]
          }
        }
        #cluster5-max_A90
        {
          data_filt_ab_cluster5 =  subset(data_filt_ab,A90>apply(data_filt_ab[,c("CTRL","A30","A15","A60","A120",treat_b)],1,max))
          #��control��0���ο�
          #data_filt_ab_cluster4 = data_filt_ab_cluster4 - data_filt_ab_cluster4 $ CTRL
          #���а���A90��ֵ��������
          {
            data_filt_ab_cluster5 = data_filt_ab_cluster5[order(-data_filt_ab_cluster5 $ A90),]
          }
        }
        #cluster6-max_A120
        {
          data_filt_ab_cluster6 =  subset(data_filt_ab,A120>apply(data_filt_ab[,c("CTRL","A30","A15","A90","A60",treat_b)],1,max))
          #��control��0���ο�
          #data_filt_ab_cluster6 = data_filt_ab_cluster6 - data_filt_ab_cluster6 $ CTRL
          #���а���A120��ֵ��������
          {
            data_filt_ab_cluster6 = data_filt_ab_cluster6[order(-data_filt_ab_cluster6 $ A120),]
          }
        }
        detach(data_filt_ab)
        #�ϲ�
        data_filt_ab_cluster_a = rbind(data_filt_ab_cluster1,data_filt_ab_cluster2,data_filt_ab_cluster3,data_filt_ab_cluster4,data_filt_ab_cluster5,data_filt_ab_cluster6)
      }
      ##��B���������
      {
        attach(data_filt_ab)
        #cluster7-max_B15
        {
          data_filt_ab_cluster7 =  subset(data_filt_ab,B15>apply(data_filt_ab[,c("CTRL","B30","B60","B90","B120",treat_a)],1,max))
          #���а���B15��ֵ��������
          {
            data_filt_ab_cluster7 = data_filt_ab_cluster7[order(-data_filt_ab_cluster7 $ B15),]
          }
        }
        #cluster8-max_B30
        {
          data_filt_ab_cluster8 = subset(data_filt_ab,B30>apply(data_filt_ab[,c("CTRL","B15","B60","B90","B120",treat_a)],1,max))
          #���а���B30��ֵ��������
          {
            data_filt_ab_cluster8 = data_filt_ab_cluster8[order(-data_filt_ab_cluster8 $ B30),]
          }
        }
        #cluster9-max_B60
        {
          data_filt_ab_cluster9 =  subset(data_filt_ab,B60>apply(data_filt_ab[,c("CTRL","B30","B15","B90","B120",treat_a)],1,max))
          #���а���B60��ֵ��������
          {
            data_filt_ab_cluster9 = data_filt_ab_cluster9[order(-data_filt_ab_cluster9 $ B60),]
          }
        }
        #cluster10-max_B90
        {
          data_filt_ab_cluster10 =  subset(data_filt_ab,B90>apply(data_filt_ab[,c("CTRL","B30","B15","B60","B120",treat_a)],1,max))
          #���а���B90��ֵ��������
          {
            data_filt_ab_cluster10 = data_filt_ab_cluster10[order(-data_filt_ab_cluster10 $ B90),]
          }
        }
        #cluster11-max_A120
        {
          data_filt_ab_cluster11 =  subset(data_filt_ab,B120>apply(data_filt_ab[,c("CTRL","B30","B15","B90","B60",treat_a)],1,max))
          #���а���B120��ֵ��������
          {
            data_filt_ab_cluster11 = data_filt_ab_cluster11[order(-data_filt_ab_cluster11 $ B120),]
          }
        }
        detach(data_filt_ab)
        #�ϲ�
        data_filt_ab_cluster_b = rbind(data_filt_ab_cluster7,data_filt_ab_cluster8,data_filt_ab_cluster9,data_filt_ab_cluster10,data_filt_ab_cluster11)
      }
      #�ϲ�����ȡlog2ֵ
      {
        #�ϲ�
        data_filt_ab_cluster = rbind(data_filt_ab_cluster_a,data_filt_ab_cluster_b)
        #��ȡlog2ֵ
        data_filt_ab_cluster = log2(data_filt_ab_cluster[,1:base::ncol(data_filt_ab_cluster)])
      }
      
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
{
  #BiocManager::install("clusterProfiler")
  #BiocManager::install("topGO")
  #BiocManager::install("Rgraphviz")
  #BiocManager::install("pathview")
  #BiocManager::install("org.Hs.eg.db")
  ##���ػ�ͼ��
  #BiocManager::install("Hmisc")
  #BiocManager::install("ggplot2")
  
  library(clusterProfiler)
  library(topGO)
  library(Rgraphviz)
  library(pathview)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(Hmisc)
  genename = c(rownames(data_filt_ab_up))
  #�鿴���е�֧�ֺͿ�ת������
  keytypes(org.Hs.eg.db)
  #����������SYMBOL��ʽת��ΪENTREZID��ʽ������org.Hs.eg.dbΪ�������ݿ�
  DEG.entrez_id = bitr(genename,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
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
  # Read in data and generate the plots
  # Biological Process
  DrawGOBubblePlot(erich.go.BP,"BP",10,"blue")
  DrawGOBubblePlot(erich.go.CC,"CC",10,"blue")
  DrawGOBubblePlot(erich.go.MF,"MF",10,"blue")
  
  
  DrawGOBubblePlot = function(erich.go.ALL,category = "Mf",)
    
    dotplot(erich.go.ALL)
  plotGOgraph(erich.go.ALL)
}



