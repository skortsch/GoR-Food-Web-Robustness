####uw new plots####
uwR50.list<-readRDS()

m<-matrix(NA,1000,34)
colnames(m)<-c(1981:2014)
row.names(m)<-c(1:1000)
mlist<-list()
for (j in 1:10){
  a<-uwR50.list[[j]]
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
colnames(durp2)<-c("iter", "year", "Con")

urp<-full_join(durp2, df.R50list[[1]])
for (i in 1:10) {
  urp<-left_join(urp, df.R50list[[i]])
}
urp
urp3<-pivot_longer(urp, 4:13)
colnames(urp3)[4]<-c("threshold")
colnames(urp3)[5]<-c("R50")
urp3

#hist plot
ggplot(urp3, aes(x=R50))+
  geom_histogram()+
  facet_wrap(~factor(threshold, levels=c("binary","10%", "20%","30%","40%",
                                                           "50%", "60%", "70%", "80%", "90%")),ncol=5)+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  ggtitle(element_text("unweighted webs R50 histplot at diffrent THs"))

ggplot(urp3, aes(sample=R50))+
  stat_qq()+
  stat_qq_line()+
  facet_wrap(~factor(threshold,levels=c("binary","10%", "20%","30%","40%",
                                        "50%", "60%", "70%", "80%", "90%")),ncol=5)+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  ggtitle(element_text("unweighted webs R50 qqplot at diffrent THs"))


ggplot(urp3, aes(x=Con))+
  geom_histogram()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())+
 # facet_wrap(~factor(threshold, levels=c("binary","10%", "20%","30%","40%",
#                                     "50%", "60%", "70%", "80%", "90%")),ncol=5)+
 # theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  ggtitle(element_text("unweighted con distribution"))



ggplot(urp3, aes(x=Con, y=R50, colour=year))+
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
  theme(plot.title = element_text(hjust = 0))+xlab("Connectance")+ylab("R50")+
  theme(strip.text.x = element_text(size = 15))+
  theme(axis.title=element_text(size=28))+theme(legend.text=element_text(size=25),
                                                legend.title=element_text(size=25))+
  facet_wrap(~factor(threshold, levels=c("binary","10%", "20%","30%","40%",
                                         "50%", "60%", "70%", "80%", "90%")),ncol=5)+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  theme(strip.text.x = element_text(size = 18))+
  theme(panel.spacing.x = unit(1, "lines"))


####new uw R50 over time####
df.r50<-matrix(NA,34,7)
colnames(df.r50)<-c("Median","upperCI95","lowerCI95","upperCI50","lowerCI50","Year","threshold")
df.r50[,6]<-c(1981:2014)
dflist<-list()
THlist<-list("binary","10%", "20%","30%","40%","50%", "60%", "70%", "80%", "90%")
for (j in 1:10){
  a<-uwR50.list[[j]]
  df.r50[,7]<-THlist[[j]][1]
  #a[1,]
  for (i in 1:34){
    df.r50[i,1]<-median(a[i,])
    df.r50[i,2]<-quantile(a[i,], probs = 0.975)
    df.r50[i,3]<-quantile(a[i,], probs = 0.05/2)
    df.r50[i,4]<-quantile(a[i,], probs = 0.75)
    df.r50[i,5]<-quantile(a[i,], probs = 0.25/2)
    dflist[[j]]<-df.r50
  }
}
dflist[[1]]
df.r50.time<-as.data.frame(dflist[[1]])
for (i in 1:10) {
  df.r50.time<-full_join(df.r50.time, as.data.frame(dflist[[i]]))
}
df.r50.time
#df.r50.time<-pivot_longer(df.r50.time, 1:5)

df.r50.time$Median<-as.numeric(as.character(df.r50.time$Median))
df.r50.time$Year<-as.numeric(as.character(df.r50.time$Year))
df.r50.time$upperCI95<-as.numeric(as.character(df.r50.time$upperCI95))
df.r50.time$lowerCI95<-as.numeric(as.character(df.r50.time$lowerCI95))
df.r50.time$lowerCI50<-as.numeric(as.character(df.r50.time$lowerCI50))
df.r50.time$upperCI50<-as.numeric(as.character(df.r50.time$upperCI50))

df.r50.time<-pivot_longer(df.r50.time, 6)
df.r50.time$value<-as.numeric(as.character(df.r50.time$value))
#filter(df.r50.time, name == "Median")

ggplot(df.r50.time, aes(x=value, y=Median, colour=threshold, fill=threshold))+
  geom_line(linewidth=1.1)+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5, size = 16))+
  geom_ribbon(aes(ymin=upperCI95, ymax=lowerCI95),alpha=0.3,colour = NA)+
  geom_ribbon(aes(ymin=upperCI50, ymax=lowerCI50),alpha=0.5,colour = NA)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Year", breaks = scales::pretty_breaks(n = 5))+
  xlab("Year")+ylab("R50")+theme(axis.title=element_text(size=25))+
  theme(legend.text=element_text(size=25),legend.title=element_text(size=25))+
  theme(legend.title = element_text("Thresholds"))+theme(legend.position = "none")+
  facet_wrap(~factor(threshold, levels=c("binary","10%", "20%","30%","40%",
                                         "50%", "60%", "70%", "80%", "90%")),ncol=5)+
  theme(strip.text.x = element_text(size = 16))+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  scale_color_brewer(palette = "Spectral", direction =1)+scale_fill_brewer(palette = "Spectral", direction =1)
# scale_fill_manual(values=c("green","orange","blue","black","red","yellow","grey","orange","dark green","dark blue"))

