library(ggplot2)
library(dplyr)
library(tidyverse)
library(NetIndices)
# set wd to data folder and load data

w.webs.1000<-readRDS("w.webs.1000.RDS")# weighted webs 1000/year
uw.webs.1000<-readRDS("index.arr.RDS") # unweighted webs 1000/year
uwR50.list<-readRDS("unweighted_R50_list.RDS")# list of matrixes with unweighted R50 values. [[1]] is the binary scenario, [[2]]= 10% scenario etc.
R50.list<-readRDS("Weighted_R50_list.RDS")# list of matrixes with with weighted values. [[1]] is the binary scenario, [[2]]= 10% scenario etc.


#### Calculate AMI and A/C###
uw.AMIind.34000<-matrix(NA,34,1000)
uw.rel.asc<-matrix(NA,34,1000)
for (i in 1:1000) {
  a<-uw.webs.1000[[i]]
  for (j in 1:34) {
    b<-a[[j]]
    ami.list<-UncInd(b)
    uw.AMIind.34000[j,i]<-ami.list$AMI/ami.list$HR # normalize the AMI to be comparable
    asc<-AscInd(b)
    uw.rel.asc[j,i]<-asc[1,4]
    #w.AC.34000[j,i]<-asc[1,4]
  }
}
uw.AMIind.34000
uw.rel.asc

####weighted version####
w.AMIind.34000<-matrix(NA,34,1000)
rel.asc<-matrix(NA,34,1000)
for (i in 1:1000) {
  a<-w.webs.1000[[i]]
  ami.list<-list()
  for (j in 1:34) {
    b<-a[[j]]
    ami.list<-UncInd(b)
    w.AMIind.34000[j,i]<-ami.list$AMI/ami.list$HR # normalize the AMI to be comparable
    asc<-AscInd(b)
    rel.asc[j,i]<-asc[1,4]
    #w.AC.34000[j,i]<-asc[1,4]
  }
}
w.AMIind.34000
rel.asc


####Make df for plot####
w.AMIind.34000
rel.asc
m<-matrix(NA,1000,34)
colnames(m)<-c(1981:2014)
row.names(m)<-c(1:1000)
mlist<-list()
for (j in 1:10){
  a<-R50.list[[j]]#all R50 values by th
  for (i in 1:34) {
    m[,i]<-a[i,]
  }
  mlist[[j]]<-m
}

THlist<-list("binary","10%", "20%","30%","40%","50%", "60%", "70%", "80%", "90%")
df.R50list<-list()
for (i in 1:10) {
  m2<-as.data.frame(as.table(mlist[[i]]))
  durp<-as_tibble(m2)
  durp$Var2<-as.numeric(as.character(durp$Var2))
  colnames(durp)<-c("iter", "year", THlist[[i]][1])
  df.R50list[[i]]<-durp
}

w.con<-readRDS("W.con_34000.RDS")
m3<-matrix(NA,1000,34)
colnames(m3)<-c(1981:2014)
row.names(m3)<-c(1:1000)
for (i in 1:34) {
  m3[,i]<-w.con[i,]
}
m3
m4<-as.data.frame(as.table(m3))
durp2<-as_tibble(m4)
durp2$Var2<-as.numeric(as.character(durp2$Var2))
colnames(durp2)<-c("iter", "year", "w.con")

urp<-full_join(durp2, df.R50list[[1]])
for (i in 1:10) {
  urp<-left_join(urp, df.R50list[[i]])
}
urp

w.AMIind.34000
m3<-matrix(NA,1000,34)
colnames(m3)<-c(1981:2014)
row.names(m3)<-c(1:1000)
for (i in 1:34) {
  m3[,i]<-w.AMIind.34000[i,]
}
m3
m4<-as.data.frame(as.table(m3))
AMI<-as_tibble(m4)
AMI$Var2<-as.numeric(as.character(AMI$Var2))
colnames(AMI)<-c("iter", "year", "AMI")


DF<-right_join(urp, AMI)

rel.asc
m3<-matrix(NA,1000,34)
colnames(m3)<-c(1981:2014)
row.names(m3)<-c(1:1000)
for (i in 1:34) {
  m3[,i]<-rel.asc[i,]
}
m3
m4<-as.data.frame(as.table(m3))
AC<-as_tibble(m4)
AC$Var2<-as.numeric(as.character(AC$Var2))
colnames(AC)<-c("iter", "year", "A/C")

DF<-right_join(DF, AC)


DF.final<-pivot_longer(DF, 4:13)
colnames(DF.final)[6]<-c("threshold")
colnames(DF.final)[7]<-c("R50")

DF.final
ggplot(DF.final, aes(x=AMI, y=R50, colour=year))+
  geom_point()+
  stat_smooth(method = "glm",
              formula = "y~x",
              geom = "smooth",
              col="Red")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 18))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  scale_color_gradient(low="yellow", high="Blue")+
  theme(plot.title = element_text(hjust = 0))+xlab("AMI")+ylab("R50")+
  theme(strip.text.x = element_text(size = 15))+
  theme(axis.title=element_text(size=28))+theme(legend.text=element_text(size=25),
                                                legend.title=element_text(size=25))+
  facet_wrap(~factor(threshold, levels=c("binary","10%", "20%","30%","40%",
                                    "50%", "60%", "70%", "80%", "90%")),ncol=5)+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  theme(strip.text.x = element_text(size = 18))
#  ggtitle(element_text("weighted webs normalized AMI & R50  at diffrent THs"))


ggplot(DF.final, aes(x=`A/C`, y=R50, colour=year))+
  geom_point()+
  stat_smooth(method = "glm",
              formula = "y~x",
              geom = "smooth",
              col="Red")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 18))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  scale_color_gradient(low="yellow", high="Blue")+
  theme(plot.title = element_text(hjust = 0))+xlab("Relative accendancy")+ylab("R50")+
  theme(strip.text.x = element_text(size = 15))+
  theme(axis.title=element_text(size=28))+theme(legend.text=element_text(size=25),
                                                legend.title=element_text(size=25))+
  facet_wrap(~factor(threshold, levels=c("binary","10%", "20%","30%","40%",
                                         "50%", "60%", "70%", "80%", "90%")),ncol=5)+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  theme(strip.text.x = element_text(size = 18))
  #ggtitle(element_text("weighted webs A/C & R50  at diffrent THs"))



####Make df for unweighted version####
#overrides the weighted version dataframe, so save weighted plot before you run
uw.AMIind.34000
uw.rel.asc


m<-matrix(NA,1000,34)
colnames(m)<-c(1981:2014)
row.names(m)<-c(1:1000)
mlist<-list()
for (j in 1:10){
  a<-uwR50.list[[j]]#all R50 values by th
  for (i in 1:34) {
    m[,i]<-a[i,]
  }
  mlist[[j]]<-m
}

THlist<-list("binary","10%", "20%","30%","40%","50%", "60%", "70%", "80%", "90%")
df.R50list<-list()
for (i in 1:10) {
  m2<-as.data.frame(as.table(mlist[[i]]))
  durp<-as_tibble(m2)
  durp$Var2<-as.numeric(as.character(durp$Var2))
  colnames(durp)<-c("iter", "year", THlist[[i]][1])
  df.R50list[[i]]<-durp
}

Con<-readRDS("Con_34000.RDS")
m3<-matrix(NA,1000,34)
colnames(m3)<-c(1981:2014)
row.names(m3)<-c(1:1000)
for (i in 1:34) {
  m3[,i]<-Con[i,]
}
m3
m4<-as.data.frame(as.table(m3))
durp2<-as_tibble(m4)
durp2$Var2<-as.numeric(as.character(durp2$Var2))
colnames(durp2)<-c("iter", "year", "con")

urp<-full_join(durp2, df.R50list[[1]])
for (i in 1:10) {
  urp<-left_join(urp, df.R50list[[i]])
}
urp

uw.AMIind.34000
m3<-matrix(NA,1000,34)
colnames(m3)<-c(1981:2014)
row.names(m3)<-c(1:1000)
for (i in 1:34) {
  m3[,i]<-uw.AMIind.34000[i,]
}
m3
m4<-as.data.frame(as.table(m3))
AMI<-as_tibble(m4)
AMI$Var2<-as.numeric(as.character(AMI$Var2))
colnames(AMI)<-c("iter", "year", "AMI")


DF<-right_join(urp, AMI)

uw.rel.asc
m3<-matrix(NA,1000,34)
colnames(m3)<-c(1981:2014)
row.names(m3)<-c(1:1000)
for (i in 1:34) {
  m3[,i]<-uw.rel.asc[i,]
}
m3
m4<-as.data.frame(as.table(m3))
AC<-as_tibble(m4)
AC$Var2<-as.numeric(as.character(AC$Var2))
colnames(AC)<-c("iter", "year", "A/C")

DF<-right_join(DF, AC)


DF.final<-pivot_longer(DF, 4:13)
colnames(DF.final)[6]<-c("threshold")
colnames(DF.final)[7]<-c("R50")

DF.final
ggplot(DF.final, aes(x=AMI, y=R50, colour=year))+
  geom_point()+
  stat_smooth(method = "glm",
              formula = "y~x",
              geom = "smooth",
              col="Red")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 18))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  scale_color_gradient(low="yellow", high="Blue")+
  theme(plot.title = element_text(hjust = 0))+xlab("AMI")+ylab("R50")+
  theme(strip.text.x = element_text(size = 15))+
  theme(axis.title=element_text(size=28))+theme(legend.text=element_text(size=25),
                                                legend.title=element_text(size=25))+
  facet_wrap(~factor(threshold, levels=c("binary","10%", "20%","30%","40%",
                                         "50%", "60%", "70%", "80%", "90%")),ncol=5)+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  theme(strip.text.x = element_text(size = 18))
  #ggtitle(element_text("unweighted webs normalized AMI & R50  at diffrent THs"))


ggplot(DF.final, aes(x=`A/C`, y=R50, colour=year))+
  geom_point()+
  stat_smooth(method = "glm",
              formula = "y~x",
              geom = "smooth",
              col="Red")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 18))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  scale_color_gradient(low="yellow", high="Blue")+
  theme(plot.title = element_text(hjust = 0))+xlab("Relative accendancy")+ylab("R50")+
  theme(strip.text.x = element_text(size = 15))+
  theme(axis.title=element_text(size=28))+theme(legend.text=element_text(size=25),
                                                legend.title=element_text(size=25))+
  facet_wrap(~factor(threshold, levels=c("binary","10%", "20%","30%","40%",
                                         "50%", "60%", "70%", "80%", "90%")),ncol=5)+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  theme(strip.text.x = element_text(size = 18))
