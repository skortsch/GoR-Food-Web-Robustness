####ploting R50####
uwR50<-readRDS("uw.R50.output.1000")
w.R50th90<-readRDS("w.R50th90.RDS")
w.R50th80<-readRDS("w.R50th80.RDS")
w.R50th70<-readRDS("w.R50th70.RDS")
w.R50th60<-readRDS("w.R50th60.RDS")
w.R50th50<-readRDS("w.R50th50.RDS")
w.R50th40<-readRDS("w.R50th40.RDS")
w.R50th30<-readRDS("w.R50.output.1000.RDS")
w.R50th20<-readRDS("w.R50th20.RDS")
w.R50th10<-readRDS("w.R50th10.RDS")
w.R50Binary<-readRDS("w.R50Binary.RDS")
R50.list<-list(w.R50Binary, w.R50th10, w.R50th20, w.R50th30, w.R50th40,
               w.R50th50, w.R50th60, w.R50th70, w.R50th80, w.R50th90)
saveRDS(R50.list,file = "Weighted_R50_list.RDS")
R50.list

df<-matrix(NA,34,6)
colnames(df)<-c("year","bin_median","upCI95","lowCI95", "upCI50", "lowCI50")
df[,1]<-c(seq(1981,2014))
df[,2]<-apply(w.R50Binary,1,median)
#df[,3]<-apply(w.R50th30,1,median)

# conf intervals 

for (i in 1:34){
  df[i,3]<-quantile(w.R50Binary[i,], probs = 0.975)
  df[i,4]<-quantile(w.R50Binary[i,], probs = 0.05/2)
  df[i,5]<-quantile(w.R50Binary[i,], probs = 0.75)
  df[i,6]<-quantile(w.R50Binary[i,], probs = 0.25/2)
}


#plot
df<-as.data.frame(df)
g<-ggplot(df, aes(x=year))+
  geom_ribbon(aes(ymin=upCI95, ymax=lowCI95), fill="blue",alpha=0.2)+
  geom_ribbon(aes(ymin=upCI50, ymax=lowCI50), fill="blue",alpha=0.2)+
  geom_line(aes(x=year, y=bin_median ), col="blue", linewidth=1.1)+
  xlim(c(1981,2014))+ylim(c(0, 0.55))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5,size = 5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Year", breaks = 1981:2014)+
  xlab("Year")+ylab("R50 index")+
  ggtitle("Binary")


g+ggtitle("R50 (loss of 70%= extinction) over time")+theme(plot.title = element_text(hjust = 0.5))

g+geom_line(aes(y=bin_median), col="red",linewidth=1.1)+
  geom_ribbon(aes(ymin=bin_min , ymax=bin_max), fill="red", alpha=0.2)+
  ggtitle("R50 (loss of 70%= extinction) over time")+theme(plot.title = element_text(hjust = 0.5))




range<-matrix(NA,34,2)
for (i in 1:34) {
  range[i,1]<-min(con[i,])
  range[i,2]<-max(con[i,])
}

g+geom_line(aes(x=year, y=apply(w.con,1,median) ), col="orange", linewidth=1.1)+
  geom_ribbon(aes(ymin=range[,1], ymax=range[,2]), fill="orange",alpha=0.4)+
  ggtitle("R50 (loss of 70%= extinction) and wConnectance")+theme(plot.title = element_text(hjust = 0.5))



g+geom_line(aes(x=year, y=apply(con,1,median) ), col="orange", linewidth=1.1)+
  geom_ribbon(aes(ymin=range[,1], ymax=range[,2]), fill="orange",alpha=0.4)+
  ggtitle("R50 (loss of 70%= extinction) and Connectance")+theme(plot.title = element_text(hjust = 0.5))

  
#all in one
#ggplot(df, aes(x=year))+
#  geom_line(aes(x=year, y=w_median ), col="blue", linewidth=1.1)+
#  xlim(c(1981,2014))+ylim(c(0, 0.5))+
#  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
#  scale_x_continuous(name = "Year", breaks = 1981:2014)+
#  xlab("Year")+ylab("R50 index")+
#  geom_line(aes(x=year, y=apply(w.R50Binary,1,median)), col=1, linewidth=1.1)+
#  geom_line(aes(x=year, y=apply(w.R50th20,1,median)), col=2, linewidth=1.1)+
#  geom_line(aes(x=year, y=apply(w.R50th40,1,median)), col=3, linewidth=1.1)+
#  geom_line(aes(x=year, y=apply(w.R50th50,1,median)), col=4, linewidth=1.1)+
#  geom_line(aes(x=year, y=apply(w.R50th60,1,median)), col=5, linewidth=1.1)+
#  geom_line(aes(x=year, y=apply(w.R50th70,1,median)), col=6, linewidth=1.1)+
#  geom_line(aes(x=year, y=apply(w.R50th80,1,median)), col=7, linewidth=1.1)+
#  geom_line(aes(x=year, y=apply(w.R50th90,1,median)), col=8, linewidth=1.1)


# Combined R50 plot 1
df.allTH<-matrix(NA,34,11)
df.allTH[,11]<-c(seq(1981,2014))
colnames(df.allTH)<-c("Binary","TH 10%","TH 20%", "TH 30%", "TH 40%", "TH 50%",
                      "TH 60%", "TH 70%", "TH 80%", "TH 90%", "Year")
for (i in 1:10){
df.allTH[,i]<-apply(R50.list[[i]],1, median)
}
as.data.frame(df.allTH)
d<-pivot_longer(as.data.frame(df.allTH), cols = 1:10)

ggplot(d, aes(x=Year, y=value, color = variable))+
  geom_line(aes(col=name),linewidth=1.1)+
  xlim(c(1981,2014))+ylim(c(0, 0.5))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 25))+
  theme(axis.text.y = element_text(vjust = 0.5, size = 25))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Year", breaks = scales::pretty_breaks(n = 5))+
  xlab("Year")+ylab("R50")+theme(axis.title=element_text(size=28))+
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))+
  scale_color_brewer(palette="Spectral")+
  ggtitle("Weighted food webs R50 with Thresholds")+theme(plot.title = element_text(hjust = 0.5, size = 25))+
  labs(color = "Thresholds")+
  theme(legend.text=element_text(size=25),
        legend.title=element_text(size=25))

# Combined R50 plot 2
df.allTH<-matrix(NA,34,12)
df.allTH[,12]<-c(seq(1981,2014))
colnames(df.allTH)<-c("Binary","TH 10%","TH 20%", "TH 30%", "TH 40%", "TH 50%",
                      "TH 60%", "TH 70%", "TH 80%", "TH 90%","uw", "Year")
for (i in 1:11){
  df.allTH[,i]<-apply(R50.list[[i]],1, median)
}
as.data.frame(df.allTH)
d<-pivot_longer(as.data.frame(df.allTH), cols = 1:11)

ggplot(d, aes(x=Year, y=value, color = variable))+
  geom_line(aes(col=name),linewidth=1.1)+
  xlim(c(1981,2014))+ylim(c(0, 0.5))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Year", breaks = 1981:2014)+
  xlab("Year")+ylab("R50 index")+
  scale_color_brewer(palette="Spectral")+
  ggtitle("Median R50 at diffrent thresholds")+theme(plot.title = element_text(hjust = 0.5))+
  labs(color = "Thresholds")

#without bin
ggplot(d, aes(x=Year, y=value, color = variable))+
  geom_line(aes(col=name),linewidth=1.1)+
  xlim(c(1981,2014))+ylim(c(0.06, 0.3))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Year", breaks = 1981:2014)+
  xlab("Year")+ylab("R50 index")+
  scale_color_brewer(palette="Spectral")+
  ggtitle("Median R50 at diffrent thresholds")+theme(plot.title = element_text(hjust = 0.5))+
  labs(color = "Thresholds")


####metrics plots####
df.wcon<-matrix(NA,34,6)
colnames(df.wcon)<-c("Median","upperCI95","lowerCI95","upperCI50","lowerCI50","Year")
df.wcon[,6]<-c(1981:2014)
df.wcon[,1]<-apply(w.con,1,median)
#for (i in 1:34) {
#  df.wcon[i,1]<-min(w.con[i,])
#  df.wcon[i,2]<-max(w.con[i,])
#} repalce with Conf intervals 95% and 50%

for (i in 1:34){
  df.wcon[i,2]<-quantile(w.con[i,], probs = 0.975)
  df.wcon[i,3]<-quantile(w.con[i,], probs = 0.05/2)
  df.wcon[i,4]<-quantile(w.con[i,], probs = 0.75)
  df.wcon[i,5]<-quantile(w.con[i,], probs = 0.25/2)
}
df.wcon<-as.data.frame(df.wcon)

g.wcon<-ggplot(df.wcon, aes(x=Year))+
  geom_line(aes(y=Median), col="orange", linewidth=1.1)+
  geom_ribbon(aes(ymin=upperCI95, ymax=lowerCI95), fill="orange",alpha=0.2)+
  geom_ribbon(aes(ymin=upperCI50, ymax=lowerCI50), fill="orange",alpha=0.3)+
  xlim(c(1981,2014))+ylim(c(0.05, 0.15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 15))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Year", breaks = scales::pretty_breaks(n = 5))+
  xlab("Year")+ylab("Weighted Connectance")+theme(axis.title=element_text(size=15))




df.qG.w<-matrix(NA,34,6)
colnames(df.qG.w)<-c("Median","upperCI95","lowerCI95","upperCI50","lowerCI50","Year")
df.qG.w[,6]<-c(1981:2014)
df.qG.w[,1]<-apply(qG.w,1,median)

for (i in 1:34){
  df.qG.w[i,2]<-quantile(qG.w[i,], probs = 0.975)
  df.qG.w[i,3]<-quantile(qG.w[i,], probs = 0.05/2)
  df.qG.w[i,4]<-quantile(qG.w[i,], probs = 0.75)
  df.qG.w[i,5]<-quantile(qG.w[i,], probs = 0.25/2)
}
df.qG.w<-as.data.frame(df.qG.w)

g.qG.w<-ggplot(df.qG.w, aes(x=Year))+
  geom_line(aes(y=Median), col="green", linewidth=1.1)+
  geom_ribbon(aes(ymin=upperCI95, ymax=lowerCI95), fill="green",alpha=0.2)+
  geom_ribbon(aes(ymin=upperCI50, ymax=lowerCI50), fill="green",alpha=0.3)+
  xlim(c(1981,2014))+ylim(c(1, 1.6))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Year", breaks = 1981:2014)+
  xlab("Year")+ylab("weighted Generality")


df.qV.w
df.qV.w<-matrix(NA,34,6)
colnames(df.qV.w)<-c("Median","upperCI95","lowerCI95","upperCI50","lowerCI50","Year")
df.qV.w[,6]<-c(1981:2014)
df.qV.w[,1]<-apply(qV.w,1,median)

for (i in 1:34){
  df.qV.w[i,2]<-quantile(qV.w[i,], probs = 0.975)
  df.qV.w[i,3]<-quantile(qV.w[i,], probs = 0.05/2)
  df.qV.w[i,4]<-quantile(qV.w[i,], probs = 0.75)
  df.qV.w[i,5]<-quantile(qV.w[i,], probs = 0.25/2)
}
df.qV.w<-as.data.frame(df.qV.w)

g.qV.w<-ggplot(df.qV.w, aes(x=Year))+
  geom_line(aes(y=Median), col="red", linewidth=1.1)+
  geom_ribbon(aes(ymin=upperCI95, ymax=lowerCI95), fill="red",alpha=0.2)+
  geom_ribbon(aes(ymin=upperCI50, ymax=lowerCI50), fill="red",alpha=0.3)+
  xlim(c(1981,2014))+ylim(c(1.5, 5.5))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Year", breaks = 1981:2014)+
  xlab("Year")+ylab("Weighted Vulnerability")




df.no.spec<-matrix(NA,34,6)
colnames(df.no.spec)<-c("Median","upperCI95","lowerCI95","upperCI50","lowerCI50","Year")
df.no.spec[,6]<-c(1981:2014)
df.no.spec[,1]<-apply(no.spec,1,median)

for (i in 1:34){
  df.no.spec[i,2]<-quantile(no.spec[i,], probs = 0.975)
  df.no.spec[i,3]<-quantile(no.spec[i,], probs = 0.05/2)
  df.no.spec[i,4]<-quantile(no.spec[i,], probs = 0.75)
  df.no.spec[i,5]<-quantile(no.spec[i,], probs = 0.25/2)
}
df.no.spec<-as.data.frame(df.no.spec)

g.no.spec<-ggplot(df.no.spec, aes(x=Year))+
  geom_line(aes(y=Median), col="black", linewidth=1.1)+
  geom_ribbon(aes(ymin=upperCI95, ymax=lowerCI95), fill="black",alpha=0.2)+
  geom_ribbon(aes(ymin=upperCI50, ymax=lowerCI50), fill="black",alpha=0.3)+
  xlim(c(1981,2014))+ylim(c(19, 33))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 15))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Year", breaks = scales::pretty_breaks(n = 5))+
  ylab("Species")+theme(axis.title=element_text(size=15))



df.con<-matrix(NA,34,6)
colnames(df.con)<-c("Median","upperCI95","lowerCI95","upperCI50","lowerCI50","Year")
df.con[,6]<-c(1981:2014)
df.con[,1]<-apply(Con,1,median)

for (i in 1:34){
  df.con[i,2]<-quantile(Con[i,], probs = 0.975)
  df.con[i,3]<-quantile(Con[i,], probs = 0.05/2)
  df.con[i,4]<-quantile(Con[i,], probs = 0.75)
  df.con[i,5]<-quantile(Con[i,], probs = 0.25/2)
}
df.con<-as.data.frame(df.con)

g.con<-ggplot(df.con, aes(x=Year))+
  geom_line(aes(y=Median), col="red", linewidth=1.1)+
  geom_ribbon(aes(ymin=upperCI95, ymax=lowerCI95), fill="red",alpha=0.2)+
  geom_ribbon(aes(ymin=upperCI50, ymax=lowerCI50), fill="red",alpha=0.3)+
  xlim(c(1981,2014))+ylim(c(0.155, 0.255))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 15))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Year", breaks = scales::pretty_breaks(n = 5))+
  ylab("Connectance")+theme(axis.title=element_text(size=15))



df.links<-matrix(NA,34,6)
colnames(df.links)<-c("Median","upperCI95","lowerCI95","upperCI50","lowerCI50","Year")
df.links[,6]<-c(1981:2014)
df.links[,1]<-apply(links,1,median)

for (i in 1:34){
  df.links[i,2]<-quantile(links[i,], probs = 0.975)
  df.links[i,3]<-quantile(links[i,], probs = 0.05/2)
  df.links[i,4]<-quantile(links[i,], probs = 0.75)
  df.links[i,5]<-quantile(links[i,], probs = 0.25/2)
}
df.links<-as.data.frame(df.links)

g.L<-ggplot(df.links, aes(x=Year))+
  geom_line(aes(y=Median), col="black", linewidth=1.1)+
  geom_ribbon(aes(ymin=upperCI95, ymax=lowerCI95), fill="black",alpha=0.2)+
  geom_ribbon(aes(ymin=upperCI50, ymax=lowerCI50), fill="black",alpha=0.3)+
  xlim(c(1981,2014))+ylim(c(50, 200))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Year", breaks = 1981:2014)+
  ylab("Links")




df.deg<-matrix(NA,34,6)
colnames(df.deg)<-c("Median","upperCI95","lowerCI95","upperCI50","lowerCI50","Year")
df.deg[,6]<-c(1981:2014)
df.deg[,1]<-apply(deg,1,median)

for (i in 1:34){
  df.deg[i,2]<-quantile(deg[i,], probs = 0.975)
  df.deg[i,3]<-quantile(deg[i,], probs = 0.05/2)
  df.deg[i,4]<-quantile(deg[i,], probs = 0.75)
  df.deg[i,5]<-quantile(deg[i,], probs = 0.25/2)
}
df.deg<-as.data.frame(df.deg)

g.deg<-ggplot(df.deg, aes(x=Year))+
  geom_line(aes(y=Median), col="black", linewidth=1.1)+
  geom_ribbon(aes(ymin=upperCI95, ymax=lowerCI95), fill="black",alpha=0.2)+
  geom_ribbon(aes(ymin=upperCI50, ymax=lowerCI50), fill="black",alpha=0.3)+
  xlim(c(1981,2014))+ylim(c(7, 12.5))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Year", breaks = 1981:2014)+
  xlab(" ")+ylab("Degree")




link.dens
df.link.dens<-matrix(NA,34,6)
colnames(df.link.dens)<-c("Median","upperCI95","lowerCI95","upperCI50","lowerCI50","Year")
df.link.dens[,6]<-c(1981:2014)
df.link.dens[,1]<-apply(link.dens,1,median)

for (i in 1:34){
  df.link.dens[i,2]<-quantile(link.dens[i,], probs = 0.975)
  df.link.dens[i,3]<-quantile(link.dens[i,], probs = 0.05/2)
  df.link.dens[i,4]<-quantile(link.dens[i,], probs = 0.75)
  df.link.dens[i,5]<-quantile(link.dens[i,], probs = 0.25/2)
}
df.link.dens<-as.data.frame(df.link.dens)

g.link.dens<-ggplot(df.link.dens, aes(x=Year))+
  geom_line(aes(y=Median), col="black", linewidth=1.1)+
  geom_ribbon(aes(ymin=upperCI95, ymax=lowerCI95), fill="black",alpha=0.2)+
  geom_ribbon(aes(ymin=upperCI50, ymax=lowerCI50), fill="black",alpha=0.3)+
  xlim(c(1981,2014))+ylim(c(3.5, 6.2))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Year", breaks = 1981:2014)+
  ylab("Link density")



df.qLD.w<-matrix(NA,34,6)
colnames(df.qLD.w)<-c("Median","upperCI95","lowerCI95","upperCI50","lowerCI50","Year")
df.qLD.w[,6]<-c(1981:2014)
df.qLD.w[,1]<-apply(qLD.w,1,median)

for (i in 1:34){
  df.qLD.w[i,2]<-quantile(qLD.w[i,], probs = 0.975)
  df.qLD.w[i,3]<-quantile(qLD.w[i,], probs = 0.05/2)
  df.qLD.w[i,4]<-quantile(qLD.w[i,], probs = 0.75)
  df.qLD.w[i,5]<-quantile(qLD.w[i,], probs = 0.25/2)
}
df.qLD.w<-as.data.frame(df.qLD.w)

g.LDw<-ggplot(df.qLD.w, aes(x=Year))+
  geom_line(aes(y=Median), col="blue", linewidth=1.1)+
  geom_ribbon(aes(ymin=upperCI95, ymax=lowerCI95), fill="blue",alpha=0.2)+
  geom_ribbon(aes(ymin=upperCI50, ymax=lowerCI50), fill="blue",alpha=0.3)+
  xlim(c(1981,2014))+ylim(c(1.2, 3.55))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Year", breaks = 1981:2014)+
  ylab("Weighted link density")


V
df.V<-matrix(NA,34,6)
colnames(df.V)<-c("Median","upperCI95","lowerCI95","upperCI50","lowerCI50","Year")
df.V[,6]<-c(1981:2014)
df.V[,1]<-apply(V,1,median)

for (i in 1:34){
  df.V[i,2]<-quantile(V[i,], probs = 0.975)
  df.V[i,3]<-quantile(V[i,], probs = 0.05/2)
  df.V[i,4]<-quantile(V[i,], probs = 0.75)
  df.V[i,5]<-quantile(V[i,], probs = 0.25/2)
}
df.V<-as.data.frame(df.V)

g.V<-ggplot(df.V, aes(x=Year))+
  geom_line(aes(y=Median), col="black", linewidth=1.1)+
  geom_ribbon(aes(ymin=upperCI95, ymax=lowerCI95), fill="black",alpha=0.2)+
  geom_ribbon(aes(ymin=upperCI50, ymax=lowerCI50), fill="black",alpha=0.3)+
  xlim(c(1981,2014))+ylim(c(3.5, 6.5))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Year", breaks = 1981:2014)+
  ylab("Vulnerability")

df$g.V<-"g.V"


#gen
for (i in 1:34){
  df.G[i,2]<-quantile(G[i,], probs = 0.975)
  df.G[i,3]<-quantile(G[i,], probs = 0.05/2)
  df.G[i,4]<-quantile(G[i,], probs = 0.75)
  df.G[i,5]<-quantile(G[i,], probs = 0.25/2)
}
df.G<-as.data.frame(df.G)

g.G<-ggplot(df.G, aes(x=Year))+
  geom_line(aes(y=Median), col="black", linewidth=1.1)+
  geom_ribbon(aes(ymin=upperCI95, ymax=lowerCI95), fill="black",alpha=0.2)+
  geom_ribbon(aes(ymin=upperCI50, ymax=lowerCI50), fill="black",alpha=0.3)+
  xlim(c(1981,2014))+ylim(c(4, 6.8))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Year", breaks = 1981:2014)+
  ylab("Generality")


####combined metricsplot####
g.no.spec+g.link.dens+g.con+g.V+g.G+g.wcon+g.qV.w+g.qG.w+g.LDw+
  plot_layout(ncol = 3)

g.con+g.V+g.G+g.link.dens+g.wcon+g.qV.w+g.qG.w+g.LDw+g.no.spec+
  plot_layout(ncol = 4)

####species, con &weighted con plot####
df.no.spec<-matrix(NA,34,6)
colnames(df.no.spec)<-c("Median","upperCI95","lowerCI95","upperCI50","lowerCI50","Year")
df.no.spec[,6]<-c(1981:2014)
df.no.spec[,1]<-apply(no.spec,1,median)

for (i in 1:34){
  df.no.spec[i,2]<-quantile(no.spec[i,], probs = 0.975)
  df.no.spec[i,3]<-quantile(no.spec[i,], probs = 0.05/2)
  df.no.spec[i,4]<-quantile(no.spec[i,], probs = 0.75)
  df.no.spec[i,5]<-quantile(no.spec[i,], probs = 0.25/2)
}
df.no.spec<-as.data.frame(df.no.spec)

g.no.spec<-ggplot(df.no.spec, aes(x=Year))+
  geom_line(aes(y=Median), col="black", linewidth=1.1)+
  geom_ribbon(aes(ymin=upperCI95, ymax=lowerCI95), fill="black",alpha=0.2)+
  geom_ribbon(aes(ymin=upperCI50, ymax=lowerCI50), fill="black",alpha=0.3)+
  xlim(c(1981,2014))+ylim(c(19, 33))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 20))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 20))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Year", breaks = scales::pretty_breaks(n = 5))+
  ylab("Species")+theme(axis.title=element_text(size=20))

df.con<-matrix(NA,34,6)
colnames(df.con)<-c("Median","upperCI95","lowerCI95","upperCI50","lowerCI50","Year")
df.con[,6]<-c(1981:2014)
df.con[,1]<-apply(Con,1,median)

for (i in 1:34){
  df.con[i,2]<-quantile(Con[i,], probs = 0.975)
  df.con[i,3]<-quantile(Con[i,], probs = 0.05/2)
  df.con[i,4]<-quantile(Con[i,], probs = 0.75)
  df.con[i,5]<-quantile(Con[i,], probs = 0.25/2)
}
df.con<-as.data.frame(df.con)
g.con<-ggplot(df.con, aes(x=Year))+
  geom_line(aes(y=Median), col="red", linewidth=1.1)+
  geom_ribbon(aes(ymin=upperCI95, ymax=lowerCI95), fill="red",alpha=0.2)+
  geom_ribbon(aes(ymin=upperCI50, ymax=lowerCI50), fill="red",alpha=0.3)+
  xlim(c(1981,2014))+ylim(c(0.155, 0.255))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 20))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 20))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Year", breaks = scales::pretty_breaks(n = 5))+
  ylab("Connectance")+theme(axis.title=element_text(size=20))

df.wcon<-matrix(NA,34,6)
colnames(df.wcon)<-c("Median","upperCI95","lowerCI95","upperCI50","lowerCI50","Year")
df.wcon[,6]<-c(1981:2014)
df.wcon[,1]<-apply(w.con,1,median)
#for (i in 1:34) {
#  df.wcon[i,1]<-min(w.con[i,])
#  df.wcon[i,2]<-max(w.con[i,])
#} repalce with Conf intervals 95% and 50%

for (i in 1:34){
  df.wcon[i,2]<-quantile(w.con[i,], probs = 0.975)
  df.wcon[i,3]<-quantile(w.con[i,], probs = 0.05/2)
  df.wcon[i,4]<-quantile(w.con[i,], probs = 0.75)
  df.wcon[i,5]<-quantile(w.con[i,], probs = 0.25/2)
}
df.wcon<-as.data.frame(df.wcon)

g.wcon<-ggplot(df.wcon, aes(x=Year))+
  geom_line(aes(y=Median), col="orange", linewidth=1.1)+
  geom_ribbon(aes(ymin=upperCI95, ymax=lowerCI95), fill="orange",alpha=0.2)+
  geom_ribbon(aes(ymin=upperCI50, ymax=lowerCI50), fill="orange",alpha=0.3)+
  xlim(c(1981,2014))+ylim(c(0.05, 0.15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 20))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 20))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Year", breaks = scales::pretty_breaks(n = 5))+
  xlab("Year")+ylab("Weighted Connectance")+theme(axis.title=element_text(size=20))


g.no.spec|g.con|g.wcon+theme(aspect.ratio = 1)





###Correlation plots####
#wcon r50
lm<-lm(apply(w.R50th30,1,median)~apply(w.con,1,median))
plot(apply(w.con,1,median),apply(w.R50th30,1,median))
abline(lm(apply(w.R50th30,1,median)~apply(w.con,1,median)))
summary(lm)
#ggplot wcon R50
data<-data.frame(apply(w.R50th30,1,median), apply(w.con,1,median))
colnames(data)<-c("R50", "w.con")
ggplot(data, aes(w.con, R50))+
  geom_point()+
  stat_smooth(method = "lm",
                  formula = "y~x",
                  geom = "smooth",
                  col="black")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Weighted connectance")+
  ggtitle("R50 threshold 30% vs Weighted connectance")

####all in one loop wcon vs R50####
R50THplots<-list()
for (i in 1:10){
  data<-data.frame(apply(R50.list[[i]],1,median), apply(w.con,1,median))
  colnames(data)<-c("R50", "w.con")
  R50THplots[i]<-list(ggplot(data, aes(w.con, R50))+
    geom_point()+
    stat_smooth(method = "lm",
                formula = "y~x",
                geom = "smooth",
                col="black")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    xlab("Weighted connectance"))
}
R50THplots[[1]]+ggtitle("TH Binary")+R50THplots[[2]]+ggtitle("TH 10%")+R50THplots[[3]]+
  ggtitle("TH 20%")+R50THplots[[4]]+ggtitle("TH 30%")+R50THplots[[5]]+ggtitle("TH 40%")+
  R50THplots[[6]]+ggtitle("TH 50%")+R50THplots[[7]]+ggtitle("TH 60%")+
  R50THplots[[8]]+ggtitle("TH 70%")+R50THplots[[9]]+ggtitle("TH 80%")+
  R50THplots[[10]]+ggtitle("TH 90%")

####tot data Con vs R50 weighted & THS####
R50THplots<-list()
for (i in 1:10){
  data<-data.frame(as.vector(R50.list[[i]]), as.vector(w.con))
  colnames(data)<-c("R50", "w.con")
  R50THplots[i]<-list(ggplot(data, aes(w.con, R50))+
                        geom_point(alpha=1)+
                        stat_smooth(method = "lm",
                                    formula = "y~x",
                                    geom = "smooth",
                                    col="red")+
                        scale_x_continuous(labels = scales::number_format(accuracy = 0.01), 
                                           breaks = scales::pretty_breaks(n = 5))+
                        theme(axis.text.x = element_text(angle = 45,vjust = 0.5, size = 15))+
                        theme(axis.text.y = element_text(vjust = 0.5, size = 15))+
                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black"))+
                        xlab("Weighted connectance")+ylab("R50")+
                        theme(axis.title=element_text(size=15)))
}
R50THplots[[1]]+ggtitle("TH Binary")+R50THplots[[2]]+ggtitle("TH 10%")+R50THplots[[3]]+
  ggtitle("TH 20%")+R50THplots[[4]]+ggtitle("TH 30%")+R50THplots[[5]]+ggtitle("TH 40%")+
  R50THplots[[6]]+ggtitle("TH 50%")+R50THplots[[7]]+ggtitle("TH 60%")+
  R50THplots[[8]]+ggtitle("TH 70%")+R50THplots[[9]]+ggtitle("TH 80%")+
  R50THplots[[10]]+ggtitle("TH 90%")


#weighted gen
lm<-lm(apply(w.R50th30,1,median)~apply(qG.w,1,median))
plot(apply(qG.w,1,median),apply(w.R50th30,1,median))
abline(lm(apply(w.R50th30,1,median)~apply(qG.w,1,median)))
summary(lm)
#
data<-data.frame(apply(w.R50th30,1,median), apply(qG.w,1,median))
colnames(data)<-c("R50", "qG.w")
ggplot(data, aes(qG.w, R50))+
  geom_point()+
  stat_smooth(method = "lm",
              formula = "y~x",
              geom = "smooth",
              col="black")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Weighted Generality")+
  ggtitle("R50 threshold 30% vs Weighted Generality")

R50THplots<-list()
for (i in 1:10){
  data<-data.frame(apply(R50.list[[i]],1,median), apply(qG.w,1,median))
  colnames(data)<-c("R50", "qG.w")
  R50THplots[i]<-list(ggplot(data, aes(qG.w, R50))+
                        geom_point()+
                        stat_smooth(method = "lm",
                                    formula = "y~x",
                                    geom = "smooth",
                                    col="black")+
                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black"))+
                        xlab("Weighted Generality"))
}
R50THplots[[1]]+ggtitle("TH Binary")+R50THplots[[2]]+ggtitle("TH 10%")+R50THplots[[3]]+
  ggtitle("TH 20%")+R50THplots[[4]]+ggtitle("TH 30%")+R50THplots[[5]]+ggtitle("TH 40%")+
  R50THplots[[6]]+ggtitle("TH 50%")+R50THplots[[7]]+ggtitle("TH 60%")+
  R50THplots[[8]]+ggtitle("TH 70%")+R50THplots[[9]]+ggtitle("TH 80%")+
  R50THplots[[10]]+ggtitle("TH 90%")



#generality
lm<-lm(apply(w.R50th30,1,median)~apply(G,1,median))
plot(apply(G,1,median),apply(w.R50th30,1,median))
abline(lm(apply(w.R50th30,1,median)~apply(G,1,median)))
summary(lm)
#
R50THplots<-list()
for (i in 1:10){
  data<-data.frame(apply(R50.list[[i]],1,median), apply(G,1,median))
  colnames(data)<-c("R50", "G")
  R50THplots[i]<-list(ggplot(data, aes(G, R50))+
                        geom_point()+
                        stat_smooth(method = "lm",
                                    formula = "y~x",
                                    geom = "smooth",
                                    col="black")+
                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black"))+
                        xlab("Generality"))
}
R50THplots[[1]]+ggtitle("TH Binary")+R50THplots[[2]]+ggtitle("TH 10%")+R50THplots[[3]]+
  ggtitle("TH 20%")+R50THplots[[4]]+ggtitle("TH 30%")+R50THplots[[5]]+ggtitle("TH 40%")+
  R50THplots[[6]]+ggtitle("TH 50%")+R50THplots[[7]]+ggtitle("TH 60%")+
  R50THplots[[8]]+ggtitle("TH 70%")+R50THplots[[9]]+ggtitle("TH 80%")+
  R50THplots[[10]]+ggtitle("TH 90%")



#weighted vuln
lm<-lm(apply(w.R50th30,1,median)~apply(qV.w,1,median))
plot(apply(qV.w,1,median),apply(w.R50th30,1,median))
abline(lm(apply(w.R50th30,1,median)~apply(qV.w,1,median)))
summary(lm)
#
data<-data.frame(apply(w.R50th30,1,median), apply(qV.w,1,median))
colnames(data)<-c("R50", "qV.w")
ggplot(data, aes(qV.w, R50))+
  geom_point()+
  stat_smooth(method = "lm",
              formula = "y~x",
              geom = "smooth",
              col="black")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Weighted vulnerability")+
  ggtitle("R50 threshold 30% vs Weighted vulnerability")

R50THplots<-list()
for (i in 1:10){
  data<-data.frame(apply(R50.list[[i]],1,median), apply(qV.w,1,median))
  colnames(data)<-c("R50", "qV.w")
  R50THplots[i]<-list(ggplot(data, aes(qV.w, R50))+
                        geom_point()+
                        stat_smooth(method = "lm",
                                    formula = "y~x",
                                    geom = "smooth",
                                    col="black")+
                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black"))+
                        xlab("Weighted vulnerability"))
}
R50THplots[[1]]+ggtitle("TH Binary")+R50THplots[[2]]+ggtitle("TH 10%")+R50THplots[[3]]+
  ggtitle("TH 20%")+R50THplots[[4]]+ggtitle("TH 30%")+R50THplots[[5]]+ggtitle("TH 40%")+
  R50THplots[[6]]+ggtitle("TH 50%")+R50THplots[[7]]+ggtitle("TH 60%")+
  R50THplots[[8]]+ggtitle("TH 70%")+R50THplots[[9]]+ggtitle("TH 80%")+
  R50THplots[[10]]+ggtitle("TH 90%")

#vuln
lm<-lm(apply(w.R50th30,1,median)~apply(V,1,median))
plot(apply(V,1,median),apply(w.R50th30,1,median))
abline(lm(apply(w.R50th30,1,median)~apply(V,1,median)))
summary(lm)
#
R50THplots<-list()
for (i in 1:10){
  data<-data.frame(apply(R50.list[[i]],1,median), apply(V,1,median))
  colnames(data)<-c("R50", "V")
  R50THplots[i]<-list(ggplot(data, aes(V, R50))+
                        geom_point()+
                        stat_smooth(method = "lm",
                                    formula = "y~x",
                                    geom = "smooth",
                                    col="black")+
                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black"))+
                        xlab("Vulnerability"))
}
R50THplots[[1]]+ggtitle("TH Binary")+R50THplots[[2]]+ggtitle("TH 10%")+R50THplots[[3]]+
  ggtitle("TH 20%")+R50THplots[[4]]+ggtitle("TH 30%")+R50THplots[[5]]+ggtitle("TH 40%")+
  R50THplots[[6]]+ggtitle("TH 50%")+R50THplots[[7]]+ggtitle("TH 60%")+
  R50THplots[[8]]+ggtitle("TH 70%")+R50THplots[[9]]+ggtitle("TH 80%")+
  R50THplots[[10]]+ggtitle("TH 90%")


#links
lm<-lm(apply(w.R50th30,1,median)~apply(links,1,median))
plot(apply(links,1,median),apply(w.R50th30,1,median))
abline(lm(apply(w.R50th30,1,median)~apply(links,1,median)))
summary(lm)
#
R50THplots<-list()
for (i in 1:10){
  data<-data.frame(apply(R50.list[[i]],1,median), apply(links,1,median))
  colnames(data)<-c("R50", "links")
  R50THplots[i]<-list(ggplot(data, aes(links, R50))+
                        geom_point()+
                        stat_smooth(method = "lm",
                                    formula = "y~x",
                                    geom = "smooth",
                                    col="black")+
                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black"))+
                        xlab("Number of links"))
}
R50THplots[[1]]+ggtitle("TH Binary")+R50THplots[[2]]+ggtitle("TH 10%")+R50THplots[[3]]+
  ggtitle("TH 20%")+R50THplots[[4]]+ggtitle("TH 30%")+R50THplots[[5]]+ggtitle("TH 40%")+
  R50THplots[[6]]+ggtitle("TH 50%")+R50THplots[[7]]+ggtitle("TH 60%")+
  R50THplots[[8]]+ggtitle("TH 70%")+R50THplots[[9]]+ggtitle("TH 80%")+
  R50THplots[[10]]+ggtitle("TH 90%")


#qLD.w
lm<-lm(apply(w.R50th30,1,median)~apply(qLD.w,1,median))
plot(apply(qLD.w,1,median),apply(w.R50th30,1,median))
abline(lm(apply(w.R50th30,1,median)~apply(qLD.w,1,median)))
summary(lm)
#
data<-data.frame(apply(w.R50th30,1,median), apply(qLD.w,1,median))
colnames(data)<-c("R50", "qLD.w")
ggplot(data, aes(qLD.w, R50))+
  geom_point()+
  stat_smooth(method = "lm",
              formula = "y~x",
              geom = "smooth",
              col="black")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Weighted Link density")+
  ggtitle("R50 threshold 30% vs Weighted Link Density")

R50THplots<-list()
for (i in 1:10){
  data<-data.frame(apply(R50.list[[i]],1,median), apply(qLD.w,1,median))
  colnames(data)<-c("R50", "qLD.w")
  R50THplots[i]<-list(ggplot(data, aes(qLD.w, R50))+
                        geom_point()+
                        stat_smooth(method = "lm",
                                    formula = "y~x",
                                    geom = "smooth",
                                    col="black")+
                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black"))+
                        xlab("Weighted Link denisty"))
}
R50THplots[[1]]+ggtitle("TH Binary")+R50THplots[[2]]+ggtitle("TH 10%")+R50THplots[[3]]+
  ggtitle("TH 20%")+R50THplots[[4]]+ggtitle("TH 30%")+R50THplots[[5]]+ggtitle("TH 40%")+
  R50THplots[[6]]+ggtitle("TH 50%")+R50THplots[[7]]+ggtitle("TH 60%")+
  R50THplots[[8]]+ggtitle("TH 70%")+R50THplots[[9]]+ggtitle("TH 80%")+
  R50THplots[[10]]+ggtitle("TH 90%")

#link.dens
lm<-lm(apply(w.R50th30,1,median)~apply(link.dens,1,median))
plot(apply(link.dens,1,median),apply(w.R50th30,1,median))
abline(lm(apply(w.R50th30,1,median)~apply(link.dens,1,median)))
summary(lm)
#
R50THplots<-list()
for (i in 1:10){
  data<-data.frame(apply(R50.list[[i]],1,median), apply(link.dens,1,median))
  colnames(data)<-c("R50", "link.dens")
  R50THplots[i]<-list(ggplot(data, aes(link.dens, R50))+
                        geom_point()+
                        stat_smooth(method = "lm",
                                    formula = "y~x",
                                    geom = "smooth",
                                    col="black")+
                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black"))+
                        xlab("Link denisty"))
}
R50THplots[[1]]+ggtitle("TH Binary")+R50THplots[[2]]+ggtitle("TH 10%")+R50THplots[[3]]+
  ggtitle("TH 20%")+R50THplots[[4]]+ggtitle("TH 30%")+R50THplots[[5]]+ggtitle("TH 40%")+
  R50THplots[[6]]+ggtitle("TH 50%")+R50THplots[[7]]+ggtitle("TH 60%")+
  R50THplots[[8]]+ggtitle("TH 70%")+R50THplots[[9]]+ggtitle("TH 80%")+
  R50THplots[[10]]+ggtitle("TH 90%")



#qL
lm<-lm(apply(w.R50th30,1,median)~apply(qL,1,median))
plot(apply(qL,1,median),apply(w.R50th30,1,median))
abline(lm(apply(w.R50th30,1,median)~apply(qL,1,median)))
summary(lm)
#
R50THplots<-list()
for (i in 1:10){
  data<-data.frame(apply(R50.list[[i]],1,median), apply(qL,1,median))
  colnames(data)<-c("R50", "qL")
  R50THplots[i]<-list(ggplot(data, aes(qL, R50))+
                        geom_point()+
                        stat_smooth(method = "lm",
                                    formula = "y~x",
                                    geom = "smooth",
                                    col="black")+
                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black"))+
                        xlab("qualitative links"))
}
R50THplots[[1]]+ggtitle("TH Binary")+R50THplots[[2]]+ggtitle("TH 10%")+R50THplots[[3]]+
  ggtitle("TH 20%")+R50THplots[[4]]+ggtitle("TH 30%")+R50THplots[[5]]+ggtitle("TH 40%")+
  R50THplots[[6]]+ggtitle("TH 50%")+R50THplots[[7]]+ggtitle("TH 60%")+
  R50THplots[[8]]+ggtitle("TH 70%")+R50THplots[[9]]+ggtitle("TH 80%")+
  R50THplots[[10]]+ggtitle("TH 90%")

#no.spec
lm<-lm(apply(w.R50th30,1,median)~apply(no.spec,1,median))
plot(apply(no.spec,1,median),apply(w.R50th30,1,median))
abline(lm(apply(w.R50th30,1,median)~apply(no.spec,1,median)))
summary(lm)
#
R50THplots<-list()
for (i in 1:10){
  data<-data.frame(apply(R50.list[[i]],1,median), apply(no.spec,1,median))
  colnames(data)<-c("R50", "no.spec")
  R50THplots[i]<-list(ggplot(data, aes(no.spec, R50))+
                        geom_point()+
                        stat_smooth(method = "lm",
                                    formula = "y~x",
                                    geom = "smooth",
                                    col="black")+
                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black"))+
                        xlab("Number of species"))
}
R50THplots[[1]]+ggtitle("TH Binary")+R50THplots[[2]]+ggtitle("TH 10%")+R50THplots[[3]]+
  ggtitle("TH 20%")+R50THplots[[4]]+ggtitle("TH 30%")+R50THplots[[5]]+ggtitle("TH 40%")+
  R50THplots[[6]]+ggtitle("TH 50%")+R50THplots[[7]]+ggtitle("TH 60%")+
  R50THplots[[8]]+ggtitle("TH 70%")+R50THplots[[9]]+ggtitle("TH 80%")+
  R50THplots[[10]]+ggtitle("TH 90%")


#degree
lm<-lm(apply(w.R50th30,1,median)~apply(deg,1,median))
plot(apply(deg,1,median),apply(w.R50th30,1,median))
abline(lm(apply(w.R50th30,1,median)~apply(deg,1,median)))
summary(lm)
#
R50THplots<-list()
for (i in 1:10){
  data<-data.frame(apply(R50.list[[i]],1,median), apply(deg,1,median))
  colnames(data)<-c("R50", "deg")
  R50THplots[i]<-list(ggplot(data, aes(deg, R50))+
                        geom_point()+
                        stat_smooth(method = "lm",
                                    formula = "y~x",
                                    geom = "smooth",
                                    col="black")+
                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black"))+
                        xlab("Mean Degree"))
}
R50THplots[[1]]+ggtitle("TH Binary")+R50THplots[[2]]+ggtitle("TH 10%")+R50THplots[[3]]+
  ggtitle("TH 20%")+R50THplots[[4]]+ggtitle("TH 30%")+R50THplots[[5]]+ggtitle("TH 40%")+
  R50THplots[[6]]+ggtitle("TH 50%")+R50THplots[[7]]+ggtitle("TH 60%")+
  R50THplots[[8]]+ggtitle("TH 70%")+R50THplots[[9]]+ggtitle("TH 80%")+
  R50THplots[[10]]+ggtitle("TH 90%")


#con R50bin
plot(apply(Con,1,median),apply(w.R50Binary,1,median))
lm<-lm(apply(w.R50Binary,1,median)~apply(Con,1,median))
abline(lm(apply(w.R50Binary,1,median)~apply(Con,1,median)))
summary(lm)
#wcon R50bin
plot(apply(w.con,1,median),apply(w.R50Binary,1,median))
lm<-lm(apply(w.R50Binary,1,median)~apply(w.con,1,median))
abline(lm(apply(w.R50Binary,1,median)~apply(w.con,1,median)))
summary(lm)

#tot dat R50 th30 vs wcon.
dat<-data.frame(as.vector(w.R50th30), as.vector(w.con))
colnames(dat)<-c("R50", "w.con")
dat
ggplot(dat, aes(w.con, R50))+
  geom_point()+
  stat_smooth(method = "lm",
              formula = "y~x",
              geom = "smooth",
              col="red")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

summary(lm(dat$R50~dat$w.con))


#to do list:
#dataframe för metrics åt susanne for plots!!! med vilken metric i df som hennes exemple

#confidens interval för metrics plots 95%, 50%  WTF???

#Nullmodel 1000 webs/year weighted

#weighted modularity av 1000 webs/year



# mail om torsdag till Riikka, Marie check
#weighted link dens check
#regression line för metrics och R50 check
