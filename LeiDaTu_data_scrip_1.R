##2020/9/27帮苟晓刚做“雷达图”。遇到的问题包括
#1，下载包，包括在bioconductor官网下载和通过devtools从GitHub下载。
#2.表格文件保存其他类型文件时，绝对不能直接复制
#3.读入数据时，数据中有UTF-8编码的中文的话会无法识别，所以在保存数据文献时要避免。
#4.将字符转化为数值类型时，若是列表需要先转化为数据框统一格式。然后通过下面的mode()=,apply()等方法进行转换。
#5.ggradar包做卫星图时只能做圆形的，做多边形的话标签会有问题。


#每次安装包前，都要使用BiocManager::install()来跟新包。记不住怎么安装包，请查看“https://bioconductor.org/install/index.html#find-bioconductor-packages”。
BiocManager::install()
BiocManager::install(c("rdevtools","usethis","ggradar"))
#"ggradar"是一个绘制火山图的软件包，需要在github上下载获取，BiocManger中没有
#“rdevtools”是R用来从GitHub获取软件包的包，需要手动安装，下载地址为“https://CRAN.R-project.org/package=devtools”。安装到一个文件夹以后，打开Rstudio右侧下方的Packages--install,在弹出的界面中，改install from为“package Archive.File(.tar.gz)”，在Browse选择刚刚安装包的地址，install。
library(rJava)
library(devtools)
library(usethis)
devtools::install_github("ricardo-bion/ggradar", dependencies=TRUE)
library(ggradar)
#删除乱七八糟的环境变量
rm(list = ls())
gc()
#查看当前路径
getwd()
#设置存储路径
setwd("E:/scientific_research/RData")
mydata = read.csv( "HUOSHANTU_data.csv",header = F) #读入数据为列表格式，无标头，以第一列为行名，此处表格中不能有中文。
mydata = mydata[-1,-1]#去掉第一行英文表头
mydata = as.matrix(mydata)#改mydata为矩阵，如此可以统一数据类型，否则下一步改字符类型为数值类型时会报错“list' object cannot be coerced to type 'double'”
#改mydata为数值类型，以便最后进行计算
#方法一：
mode(mydata) = "integer"#当确定矩阵中全为数值时可用此

#方法二：此法可直接操作编辑好的矩阵
#mydata2= matrix（as.numeric(mydata）,nrow = nrow(mydata))#出来的mydata2缺少行列名
#rownames(mydata2) = rownames(mydata2)#加行名称
#colnames(mydata2) = colnames(mydata2)#加列名称

#方法三：有待探索
#apply(madata,c(1,2),as.numeric())

colnames(mydata)=c("酸","醇","酮","芳香烃","烷烃","醛","脂","胺","含硫","萜烯") #使用需要的中文向量为矩阵列命名
#改为ggradar可识别的数据框格式
mynewdata<-data.frame(mydata)
#添加行名称
Name = c("7天","14天","21天")
mynewdata<-data.frame(Name,mynewdata)

#增加“背景（groud）”数据框
nx = dim(mynewdata)[2]-1#查看列数，[1]为行维度
da1 = matrix(c(rep(3,nx),rep(6,nx),rep(9,nx),rep(12,nx),rep(15,nx)),5,byrow = T)
da2 = data.frame(Name = c("U1","U2","U3","U4","U5"),da1)
colnames(da2) = colnames(mynewdata)
dat2 = rbind(da2,mynewdata)
dat2$Name = factor(dat2$Name,levels = dat2$Name,ordered = T)

ggradar(dat2,grid.min =0.1,grid.mid =7 ,grid.max = 16,
        #base为所有文字，axis为圈外文字，grid为圈内文字
        base.size = 1,
        values.radar = c("0", "10", "20"),
        #调整文字大小
        grid.label.size = 4,axis.label.size = 5,
        axis.labels = colnames(dat2)[-1],
        #设置背景线条和点的大小粗细
        group.point.size = c(rep(0.5,dim(dat2)[2]*5),rep(2,dim(dat2)[2]*3)),
        group.line.width = c(rep(0.5,dim(dat2)[2]*5),rep(0.8,dim(dat2)[2]*3)),
        background.circle.transparency = 0.1,
        #设置标签大小和位置
        legend.text.size = 12,legend.position = "top",
        background.circle.colour = "white",
        gridline.min.colour = "grey80",
        gridline.mid.colour = "white",
        gridline.max.colour = "white",
        #设置背景颜色
        group.colours = c(rep("grey80",5),"red","blue","yellow"),
        plot.legend = F,
        )
