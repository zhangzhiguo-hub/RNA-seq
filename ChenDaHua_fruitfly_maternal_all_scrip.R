###ballgown包处理数据,输出结果为bg_filt_trans，maxrow为>20
{
  #加载ballgown包
  library(ballgown)
  library(stringr)
  #查看当前路径
  getwd()
  #设置存储路径
  setwd("E:/scientific_research/Transcriptome_analysis/fruitfly_maternal_allrna")
  #重新查看更改过后的当前路径
  getwd()
  #读取定量分析的表型数据（批量操作型）
  list.files("./ballgown2")#查看文件名
  #字符串连接
  a = paste(rep("ballgown2",length(list.files("./ballgown2"))),#第一个向量
            list.files("./ballgown2"),#第二个向量
            sep = "/" #连接符
            )
  bg = ballgown(dataDir = "./ballgown2",samples = a )
  #过滤掉表达量低的基因,即每一行所代表的基因，其后的各值都为0的
  bg_filt=as.data.frame(gexpr(bg))
  bg_filt_trans = subset(bg_filt,apply(bg_filt,1,max)>20,genomesubset=TRUE)
  colnames(bg_filt_trans)
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

####读出并保存数据
{
  write.table(bg_filt_trans,file = "fruitfly_maternal_allrna_bg_filt_trans.csv",sep = ",",row.names = T) #列表保存为.csv格式
  write.table(bg_filt_trans,file = "fruitfly_maternal_allrna_bg_filt_trans.txt",sep = " ",row.names = T) #列表保存为.txt格式
}
getwd()
#设置存储路径
setwd("E:/scientific_research/Transcriptome_analysis/fruitfly_maternal_allrna")
#重新查看更改过后的当前路径
getwd()

#读入数据为列表格式，以第一列为行名称
{
  data = read.table("fruitfly_maternal_allrna_bg_filt_trans.txt",
                    header = T,
                    row.names = 1,
                    check.names = T#检查变量名（列）是否符合规则，不符合的话自动修改
                    )
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
    head(data,1)#查看首行
  }
  #给列重排序
  {
    #dplyr包，该宏包功能强大。是数据科学tidyverse集合的核心部件之一。学习链接https://bookdown.org/wangminjie/R4DS/colwise.html
    library(dplyr,warn.conflicts = F)
    data = data %>% select("WT0_1h_1","WT0_1h_2","WT0_1h_3","WT2_3h_1","WT2_3h_2","WT2_3h_3","WT5_6h_1","WT5_6h_2","WT5_6h_3")
  }
  }
  
#取均值，输出newdata
#成功#！要比对重复的时勿执行
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

#与隔壁突变型比较，筛除突变型的genelist,查看剩下的非编码RNA
{
  a = read.table("E:/scientific_research/Transcriptome_analysis/fruitfly_maternal/fruitfly_maternal_bg_filt_trans.txt",row.names = 1,header = T)
  a[a==0] = NA
  a = na.omit(a)
  a = intersect(rownames(a),rownames(data))
  data_cluster = newdata[!rownames(newdata) %in% a,]
  #读出
  write.csv(data_cluster,file = "non-coding RNAs.csv",row.names = T)
  data_cluster = log2(data_cluster)
  data = data_cluster
}

#原始聚类
{
    #attach(data)
    #cluster
  i <- 1
  data_cluster = data.frame()
  while(i <= length(data))
    {
      a =  subset(data,data[,i] > apply(data[,colnames(data)[-i]],1,max))
      #给行按照colnames(data)[1]的值升序排序
      a = a[order(-a[i]),]
      data_cluster = rbind(data_cluster,a)
      i <- i + 1
    }
    #求取log2值
    data_cluster = log2(data_cluster[,1:base::ncol(data_cluster)])
}


#新聚类?未改
{
  #生成1:4的所有随机组合并由小到大排列,num向量
  #order()函数的排序原则，（最小，次小，。。最大）
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
  #取log2
  cluster_data = log2(cluster_data[,1:base::ncol(cluster_data)])
}


#绘热图，输入data_cluster
main = "non-coding RNAs"
{
  library(pheatmap)
  gaps_row = c()
  sum = 0
  Cluster = c()
  length_Cluster = c()
  i = 1
  while (i <= length(colnames(data_cluster))) #绘制不同的图记得改名儿
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
  #画WTallRNA聚类的图
  {
    #添加注释
    annotation_col= data.frame(Time = c(0.5,2.5,5.5))
    annotation_row = data.frame(Cluster = rep(Cluster,length_Cluster))
    color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")#调色
    #运行没反应时，运行dev.off()
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
                      #分隔
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







#筛选MT,MT两个处理组中差异表达的基因
{
  #求出a,b两个处理组每个基因的表达的最小值和最大值
  {
    #求出a,b两个处理组每个基因的表达的最小值
    data$min_a = as.matrix(apply(data[,treat_a],1,min))
    data$min_b = as.matrix(apply(data[,treat_b],1,min))
    #求出a,b两个处理组每个基因的表达的最大值
    data$max_a = as.matrix(apply(data[,treat_a],1,max))
    data$max_b = as.matrix(apply(data[,treat_b],1,max))
  }
  #筛选ab两个处理组中差异表达的基因
  {
    #筛选A组的差异表达基因
    {
      data$log2FC_a_con_1 = log2((data$max_a)/(data$CTRL))
      data$log2FC_a_con_2 = log2((data$CTRL)/(data$min_a))
      data_filt = subset(data,data$log2FC_a_con_1>2|data$log2FC_a_con_2>2)
      data_filt_a = data_filt[,c(control,treat_a)]
      ##A组聚类
      {
        attach(data_filt_a)
        #cluster1-max_CTRL
        {
          data_filt_a_cluster1 =  subset(data_filt_a,CTRL>apply(data_filt_a[,treat_a],1,max))
          #将control设0做参考
          #data_filt_a_cluster1 = data_filt_a_cluster1 - data_filt_a_cluster1 $ CTRL
          #给行按照CTRL的值升序排序
          {
            data_filt_a_cluster1 = data_filt_a_cluster1[order(-data_filt_a_cluster1$CTRL),]
          }
        }
        #cluster2-max_A15
        {
          data_filt_a_cluster2 =  subset(data_filt_a,A15>apply(data_filt_a[,c("CTRL","A30","A60","A90","A120")],1,max))
          #将control设0做参考
          #data_filt_a_cluster2 = data_filt_a_cluster2 - data_filt_a_cluster2 $ CTRL
          #给行按照CTRL的值升序排序
          {
            data_filt_a_cluster2 = data_filt_a_cluster2[order(-data_filt_a_cluster2 $ A15),]
          }
        }
        #cluster3-max_A30
        {
          data_filt_a_cluster3 = subset(data_filt_a,A30>apply(data_filt_a[,c("CTRL","A15","A60","A90","A120")],1,max))
          #将control设0做参考
          #data_filt_a_cluster2 = data_filt_a_cluster2 - data_filt_a_cluster2 $ CTRL
          #给行按照CTRL的值升序排序
          {
            data_filt_a_cluster3 = data_filt_a_cluster3[order(-data_filt_a_cluster3 $ A30),]
          }
        }
        #cluster4-max_A60
        {
          data_filt_a_cluster4 =  subset(data_filt_a,A60>apply(data_filt_a[,c("CTRL","A30","A15","A90","A120")],1,max))
          #将control设0做参考
          #data_filt_a_cluster4 = data_filt_a_cluster4 - data_filt_a_cluster4 $ CTRL
          #给行按照CTRL的值升序排序
          {
            data_filt_a_cluster4 = data_filt_a_cluster4[order(-data_filt_a_cluster4 $ A60),]
          }
        }
        #cluster5-max_A90
        {
          data_filt_a_cluster5 =  subset(data_filt_a,A90>apply(data_filt_a[,c("CTRL","A30","A15","A60","A120")],1,max))
          #将control设0做参考
          #data_filt_a_cluster4 = data_filt_a_cluster4 - data_filt_a_cluster4 $ CTRL
          #给行按照A90的值升序排序
          {
            data_filt_a_cluster5 = data_filt_a_cluster5[order(-data_filt_a_cluster5 $ A90),]
          }
        }
        #cluster6-max_A120
        {
          data_filt_a_cluster6 =  subset(data_filt_a,A120>apply(data_filt_a[,c("CTRL","A30","A15","A90","A60")],1,max))
          #将control设0做参考
          #data_filt_a_cluster6 = data_filt_a_cluster6 - data_filt_a_cluster6 $ CTRL
          #给行按照A120的值升序排序
          {
            data_filt_a_cluster6 = data_filt_a_cluster6[order(-data_filt_a_cluster6 $ A120),]
          }
        }
        detach(data_filt_a)
        #合并
        data_filt_a_cluster = rbind(data_filt_a_cluster1,data_filt_a_cluster2,data_filt_a_cluster3,data_filt_a_cluster4,data_filt_a_cluster5,data_filt_a_cluster6)
        #求取log2值
        data_filt_a_cluster = log2(data_filt_a_cluster[,1:base::ncol(data_filt_a_cluster)])
      }
    }
    #筛选B组的差异表达基因
    {
      data$log2FC_b_con_1 = log2((data$max_b)/(data$CTRL))
      data$log2FC_b_con_2 = log2((data$CTRL)/(data$min_b))
      data_filt = subset(data,data$log2FC_b_con_1>2|data$log2FC_b_con_2>2)
      data_filt_b = data_filt[,c(control,treat_b)]
      ##B组聚类
      {
        attach(data_filt_b)
        #cluster1-max_CTRL
        {
          data_filt_b_cluster1 =  subset(data_filt_b,CTRL>apply(data_filt_b[,treat_b],1,max))
          #将control设0做参考
          #data_filt_b_cluster1 = data_filt_b_cluster1 - data_filt_b_cluster1 $ CTRL
          #给行按照CTRL的值升序排序
          {
            data_filt_b_cluster1 = data_filt_b_cluster1[order(-data_filt_b_cluster1$CTRL),]
          }
        }
        #cluster2-max_B15
        {
          data_filt_b_cluster2 =  subset(data_filt_b,B15>apply(data_filt_b[,c("CTRL","B30","B60","B90","B120")],1,max))
          #将control设0做参考
          #data_filt_b_cluster2 = data_filt_b_cluster2 - data_filt_b_cluster2 $ CTRL
          #给行按照B15的值升序排序
          {
            data_filt_b_cluster2 = data_filt_b_cluster2[order(-data_filt_b_cluster2 $ B15),]
          }
        }
        #cluster3-max_B30
        {
          data_filt_b_cluster3 = subset(data_filt_b,B30>apply(data_filt_b[,c("CTRL","B15","B60","B90","B120")],1,max))
          #将control设0做参考
          #data_filt_b_cluster2 = data_filt_b_cluster2 - data_filt_b_cluster2 $ CTRL
          #给行按照B30的值升序排序
          {
            data_filt_b_cluster3 = data_filt_b_cluster3[order(-data_filt_b_cluster3 $ B30),]
          }
        }
        #cluster4-max_A60
        {
          data_filt_b_cluster4 =  subset(data_filt_b,B60>apply(data_filt_b[,c("CTRL","B30","B15","B90","B120")],1,max))
          #将control设0做参考
          #data_filt_b_cluster4 = data_filt_b_cluster4 - data_filt_b_cluster4 $ CTRL
          #给行按照B60的值升序排序
          {
            data_filt_b_cluster4 = data_filt_b_cluster4[order(-data_filt_b_cluster4 $ B60),]
          }
        }
        #cluster5-max_A90
        {
          data_filt_b_cluster5 =  subset(data_filt_b,B90>apply(data_filt_b[,c("CTRL","B30","B15","B60","B120")],1,max))
          #将control设0做参考
          #data_filt_b_cluster4 = data_filt_b_cluster4 - data_filt_b_cluster4 $ CTRL
          #给行按照A90的值升序排序
          {
            data_filt_b_cluster5 = data_filt_b_cluster5[order(-data_filt_b_cluster5 $ B90),]
          }
        }
        #cluster6-max_A120
        {
          data_filt_b_cluster6 =  subset(data_filt_b,B120>apply(data_filt_b[,c("CTRL","B30","B15","B90","B60")],1,max))
          #将control设0做参考
          #data_filt_b_cluster6 = data_filt_b_cluster6 - data_filt_b_cluster6 $ CTRL
          #给行按照B120的值升序排序
          {
            data_filt_b_cluster6 = data_filt_b_cluster6[order(-data_filt_b_cluster6 $ B120),]
          }
        }
        detach(data_filt_b)
        #合并
        data_filt_b_cluster = rbind(data_filt_b_cluster1,data_filt_b_cluster2,data_filt_b_cluster3,data_filt_b_cluster4,data_filt_b_cluster5,data_filt_b_cluster6)
        #求取log2值
        data_filt_b_cluster = log2(data_filt_b_cluster[,1:base::ncol(data_filt_b_cluster)])
      }
    }
    #筛选AB两个组的差异表达基因
    {
      data_filt = subset(data,(data$log2FC_b_con_1>2|data$log2FC_b_con_2>2)
                         &
                           (data$log2FC_a_con_1>2|data$log2FC_a_con_2>2)
      )
      data_filt_ab = data_filt[,c(control,treat_a,treat_b)]
      ##按A的情况聚类
      {
        attach(data_filt_ab)
        #cluster1-max_CTRL
        {
          data_filt_ab_cluster1 =  subset(data_filt_ab,CTRL>apply(data_filt_ab[,c(treat_a,treat_b)],1,max))
          #将control设0做参考
          #data_filt_ab_cluster1 = data_filt_ab_cluster1 - data_filt_ab_cluster1 $ CTRL
          #给行按照CTRL的值升序排序
          {
            data_filt_ab_cluster1 = data_filt_ab_cluster1[order(-data_filt_ab_cluster1$CTRL),]
          }
        }
        #cluster2-max_A15
        {
          data_filt_ab_cluster2 =  subset(data_filt_ab,A15>apply(data_filt_ab[,c("CTRL","A30","A60","A90","A120",treat_b)],1,max))
          #将control设0做参考
          #data_filt_ab_cluster2 = data_filt_ab_cluster2 - data_filt_ab_cluster2 $ CTRL
          #给行按照A15的值升序排序
          {
            data_filt_ab_cluster2 = data_filt_ab_cluster2[order(-data_filt_ab_cluster2 $ A15),]
          }
        }
        #cluster3-max_A30
        {
          data_filt_ab_cluster3 = subset(data_filt_ab,A30>apply(data_filt_ab[,c("CTRL","A15","A60","A90","A120",treat_b)],1,max))
          #将control设0做参考
          #data_filt_a_cluster2 = data_filt_a_cluster2 - data_filt_a_cluster2 $ CTRL
          #给行按照A30的值升序排序
          {
            data_filt_ab_cluster3 = data_filt_ab_cluster3[order(-data_filt_ab_cluster3 $ A30),]
          }
        }
        #cluster4-max_A60
        {
          data_filt_ab_cluster4 =  subset(data_filt_ab,A60>apply(data_filt_ab[,c("CTRL","A30","A15","A90","A120",treat_b)],1,max))
          #将control设0做参考
          #data_filt_a_cluster4 = data_filt_a_cluster4 - data_filt_a_cluster4 $ CTRL
          #给行按A60的值升序排序
          {
            data_filt_ab_cluster4 = data_filt_ab_cluster4[order(-data_filt_ab_cluster4 $ A60),]
          }
        }
        #cluster5-max_A90
        {
          data_filt_ab_cluster5 =  subset(data_filt_ab,A90>apply(data_filt_ab[,c("CTRL","A30","A15","A60","A120",treat_b)],1,max))
          #将control设0做参考
          #data_filt_ab_cluster4 = data_filt_ab_cluster4 - data_filt_ab_cluster4 $ CTRL
          #给行按照A90的值升序排序
          {
            data_filt_ab_cluster5 = data_filt_ab_cluster5[order(-data_filt_ab_cluster5 $ A90),]
          }
        }
        #cluster6-max_A120
        {
          data_filt_ab_cluster6 =  subset(data_filt_ab,A120>apply(data_filt_ab[,c("CTRL","A30","A15","A90","A60",treat_b)],1,max))
          #将control设0做参考
          #data_filt_ab_cluster6 = data_filt_ab_cluster6 - data_filt_ab_cluster6 $ CTRL
          #给行按照A120的值升序排序
          {
            data_filt_ab_cluster6 = data_filt_ab_cluster6[order(-data_filt_ab_cluster6 $ A120),]
          }
        }
        detach(data_filt_ab)
        #合并
        data_filt_ab_cluster_a = rbind(data_filt_ab_cluster1,data_filt_ab_cluster2,data_filt_ab_cluster3,data_filt_ab_cluster4,data_filt_ab_cluster5,data_filt_ab_cluster6)
      }
      ##按B的情况聚类
      {
        attach(data_filt_ab)
        #cluster7-max_B15
        {
          data_filt_ab_cluster7 =  subset(data_filt_ab,B15>apply(data_filt_ab[,c("CTRL","B30","B60","B90","B120",treat_a)],1,max))
          #给行按照B15的值升序排序
          {
            data_filt_ab_cluster7 = data_filt_ab_cluster7[order(-data_filt_ab_cluster7 $ B15),]
          }
        }
        #cluster8-max_B30
        {
          data_filt_ab_cluster8 = subset(data_filt_ab,B30>apply(data_filt_ab[,c("CTRL","B15","B60","B90","B120",treat_a)],1,max))
          #给行按照B30的值升序排序
          {
            data_filt_ab_cluster8 = data_filt_ab_cluster8[order(-data_filt_ab_cluster8 $ B30),]
          }
        }
        #cluster9-max_B60
        {
          data_filt_ab_cluster9 =  subset(data_filt_ab,B60>apply(data_filt_ab[,c("CTRL","B30","B15","B90","B120",treat_a)],1,max))
          #给行按照B60的值升序排序
          {
            data_filt_ab_cluster9 = data_filt_ab_cluster9[order(-data_filt_ab_cluster9 $ B60),]
          }
        }
        #cluster10-max_B90
        {
          data_filt_ab_cluster10 =  subset(data_filt_ab,B90>apply(data_filt_ab[,c("CTRL","B30","B15","B60","B120",treat_a)],1,max))
          #给行按照B90的值升序排序
          {
            data_filt_ab_cluster10 = data_filt_ab_cluster10[order(-data_filt_ab_cluster10 $ B90),]
          }
        }
        #cluster11-max_A120
        {
          data_filt_ab_cluster11 =  subset(data_filt_ab,B120>apply(data_filt_ab[,c("CTRL","B30","B15","B90","B60",treat_a)],1,max))
          #给行按照B120的值升序排序
          {
            data_filt_ab_cluster11 = data_filt_ab_cluster11[order(-data_filt_ab_cluster11 $ B120),]
          }
        }
        detach(data_filt_ab)
        #合并
        data_filt_ab_cluster_b = rbind(data_filt_ab_cluster7,data_filt_ab_cluster8,data_filt_ab_cluster9,data_filt_ab_cluster10,data_filt_ab_cluster11)
      }
      #合并并且取log2值
      {
        #合并
        data_filt_ab_cluster = rbind(data_filt_ab_cluster_a,data_filt_ab_cluster_b)
        #求取log2值
        data_filt_ab_cluster = log2(data_filt_ab_cluster[,1:base::ncol(data_filt_ab_cluster)])
      }
      
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
{
  #BiocManager::install("clusterProfiler")
  #BiocManager::install("topGO")
  #BiocManager::install("Rgraphviz")
  #BiocManager::install("pathview")
  #BiocManager::install("org.Hs.eg.db")
  ##下载绘图包
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
  #查看所有的支持和可转化类型
  keytypes(org.Hs.eg.db)
  #将基因名的SYMBOL格式转化为ENTREZID格式，其中org.Hs.eg.db为人类数据库
  DEG.entrez_id = bitr(genename,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
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
  # Read in data and generate the plots
  # Biological Process
  DrawGOBubblePlot(erich.go.BP,"BP",10,"blue")
  DrawGOBubblePlot(erich.go.CC,"CC",10,"blue")
  DrawGOBubblePlot(erich.go.MF,"MF",10,"blue")
  
  
  DrawGOBubblePlot = function(erich.go.ALL,category = "Mf",)
    
    dotplot(erich.go.ALL)
  plotGOgraph(erich.go.ALL)
}




