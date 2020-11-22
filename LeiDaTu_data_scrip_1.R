##2020/9/27�ﹶ���������״�ͼ�����������������
#1�����ذ���������bioconductor�������غ�ͨ��devtools��GitHub���ء�
#2.�����ļ��������������ļ�ʱ�����Բ���ֱ�Ӹ���
#3.��������ʱ����������UTF-8��������ĵĻ����޷�ʶ�������ڱ�����������ʱҪ���⡣
#4.���ַ�ת��Ϊ��ֵ����ʱ�������б���Ҫ��ת��Ϊ���ݿ�ͳһ��ʽ��Ȼ��ͨ�������mode()=,apply()�ȷ�������ת����
#5.ggradar��������ͼʱֻ����Բ�εģ�������εĻ���ǩ�������⡣


#ÿ�ΰ�װ��ǰ����Ҫʹ��BiocManager::install()�����°����ǲ�ס��ô��װ������鿴��https://bioconductor.org/install/index.html#find-bioconductor-packages����
BiocManager::install()
BiocManager::install(c("rdevtools","usethis","ggradar"))
#"ggradar"��һ�����ƻ�ɽͼ������������Ҫ��github�����ػ�ȡ��BiocManger��û��
#��rdevtools����R������GitHub��ȡ�������İ�����Ҫ�ֶ���װ�����ص�ַΪ��https://CRAN.R-project.org/package=devtools������װ��һ���ļ����Ժ󣬴�Rstudio�Ҳ��·���Packages--install,�ڵ����Ľ����У���install fromΪ��package Archive.File(.tar.gz)������Browseѡ��ոհ�װ���ĵ�ַ��install��
library(rJava)
library(devtools)
library(usethis)
devtools::install_github("ricardo-bion/ggradar", dependencies=TRUE)
library(ggradar)
#ɾ�����߰���Ļ�������
rm(list = ls())
gc()
#�鿴��ǰ·��
getwd()
#���ô洢·��
setwd("E:/scientific_research/RData")
mydata = read.csv( "HUOSHANTU_data.csv",header = F) #��������Ϊ�б���ʽ���ޱ�ͷ���Ե�һ��Ϊ�������˴������в��������ġ�
mydata = mydata[-1,-1]#ȥ����һ��Ӣ�ı�ͷ
mydata = as.matrix(mydata)#��mydataΪ������˿���ͳһ�������ͣ�������һ�����ַ�����Ϊ��ֵ����ʱ�ᱨ����list' object cannot be coerced to type 'double'��
#��mydataΪ��ֵ���ͣ��Ա������м���
#����һ��
mode(mydata) = "integer"#��ȷ��������ȫΪ��ֵʱ���ô�

#���������˷���ֱ�Ӳ����༭�õľ���
#mydata2= matrix��as.numeric(mydata��,nrow = nrow(mydata))#������mydata2ȱ��������
#rownames(mydata2) = rownames(mydata2)#��������
#colnames(mydata2) = colnames(mydata2)#��������

#���������д�̽��
#apply(madata,c(1,2),as.numeric())

colnames(mydata)=c("��","��","ͪ","������","����","ȩ","֬","��","����","��ϩ") #ʹ����Ҫ����������Ϊ����������
#��Ϊggradar��ʶ������ݿ��ʽ
mynewdata<-data.frame(mydata)
#����������
Name = c("7��","14��","21��")
mynewdata<-data.frame(Name,mynewdata)

#���ӡ�������groud�������ݿ�
nx = dim(mynewdata)[2]-1#�鿴������[1]Ϊ��ά��
da1 = matrix(c(rep(3,nx),rep(6,nx),rep(9,nx),rep(12,nx),rep(15,nx)),5,byrow = T)
da2 = data.frame(Name = c("U1","U2","U3","U4","U5"),da1)
colnames(da2) = colnames(mynewdata)
dat2 = rbind(da2,mynewdata)
dat2$Name = factor(dat2$Name,levels = dat2$Name,ordered = T)

ggradar(dat2,grid.min =0.1,grid.mid =7 ,grid.max = 16,
        #baseΪ�������֣�axisΪȦ�����֣�gridΪȦ������
        base.size = 1,
        values.radar = c("0", "10", "20"),
        #�������ִ�С
        grid.label.size = 4,axis.label.size = 5,
        axis.labels = colnames(dat2)[-1],
        #���ñ��������͵�Ĵ�С��ϸ
        group.point.size = c(rep(0.5,dim(dat2)[2]*5),rep(2,dim(dat2)[2]*3)),
        group.line.width = c(rep(0.5,dim(dat2)[2]*5),rep(0.8,dim(dat2)[2]*3)),
        background.circle.transparency = 0.1,
        #���ñ�ǩ��С��λ��
        legend.text.size = 12,legend.position = "top",
        background.circle.colour = "white",
        gridline.min.colour = "grey80",
        gridline.mid.colour = "white",
        gridline.max.colour = "white",
        #���ñ�����ɫ
        group.colours = c(rep("grey80",5),"red","blue","yellow"),
        plot.legend = F,
        )