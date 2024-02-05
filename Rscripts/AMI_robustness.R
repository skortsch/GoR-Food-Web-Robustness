library(NetIndices)

w.webs.1000<-readRDS("w.webs.1000.RDS")# weighted webs 1000/year

uw.webs.1000<-readRDS("index.arr.RDS") # unweighted webs 1000/year
df.uw.long<-readRDS("df.unweighted.RDS") #unweighted webs df, long format for plot
df.w.long<-readRDS("df.weighted.RDS")# weighted webs df, longformat for plot

df.uw.wide<-readRDS("df.uw.wide.RDS")# unweighted webs df, wide format,
df.w.wide<-readRDS("df.w.wide.RDS")# weighted webs df, wide format
#Con= connectance,binary+10%-90% values = R50 values for the thresholds
#w.con= weighted connectance, binary+10%-90% values = R50 values for the thresholds

uw.AMIind.34000<-matrix(NA,1000,34)
colnames(uw.AMIind.34000)<-c(1981:2014)
row.names(uw.AMIind.34000)<-c(1:1000)
for (i in 1:1000) {
  a<-uw.webs.1000[[i]]
  ami.list<-list()
    for (j in 1:34) {
      b<-a[[j]]
      ami.list<-UncInd(b)
      uw.AMIind.34000[i,j]<-ami.list$AMI/ami.list$Hmax
    }
}
uw.AMIind.34000

hist(uw.AMIind.34000)# looks very close to normal distr.
qqnorm(uw.AMIind.34000)
qqline(uw.AMIind.34000)# looks very close to normal distr.

plot(uw.AMIind.34000, df.uw.wide$binary)
plot(uw.AMIind.34000, df.uw.wide$`10%`)
plot(uw.AMIind.34000, df.uw.wide$`20%`)
plot(uw.AMIind.34000, df.uw.wide$`30%`)
plot(uw.AMIind.34000, df.uw.wide$`40%`)


ttt<-df.uw.wide %>% mutate(new_col=c(as.vector(uw.AMIind.34000)))
ttt<-as_tibble(as.data.frame(as.table(uw.AMIind.34000)))
colnames(ttt)<-c("iter", "year", "norm.ami")
ttt$year<-as.numeric(as.character(ttt$year))
colnames(ttt)<-c("iter", "year", "norm.ami")
ttt<-left_join(df.uw.wide,ttt)
plot.df<-pivot_longer(ttt,4:13)

ggplot(plot.df, aes(x=norm.ami, y=value, colour=year))+
  geom_point()+
  stat_smooth(method = "lm",
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
  facet_wrap(~factor(name, levels=c("binary","10%", "20%","30%","40%",
                                         "50%", "60%", "70%", "80%", "90%")),ncol=5)+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  theme(strip.text.x = element_text(size = 18))+
  ggtitle(element_text("unweighted webs normalized AMI & R50  at diffrent THs"))



w.AMIind.34000<-matrix(NA,34,1000)
rel.asc<-matrix(NA,34,1000)
for (i in 1:1000) {
  a<-w.webs.1000[[i]]
  ami.list<-list()
  for (j in 1:34) {
    b<-a[[j]]
    ami.list<-UncInd(b)
    w.AMIind.34000[j,i]<-ami.list$AMI/ami.list$Hmax
    asc<-AscInd(b)
    rel.asc[j,i]<-asc[1,4]
    #w.AC.34000[j,i]<-asc[1,4]
  }
}
w.AMIind.34000
rel.asc

plot(w.AMIind.34000, df.w.wide$binary)
plot(w.AMIind.34000, df.w.wide$`40%`)
plot(w.AMIind.34000, df.w.wide$`90%`)

plot(w.AMIind.34000, df.w.wide$w.con)

plot(rel.asc, df.w.wide$binary)
plot(rel.asc, df.w.wide$`40%`)
plot(rel.asc, df.w.wide$`90%`)

plot(rel.asc, df.uw.wide$Con)



hist(w.AMIind.34000) # not at all normal distr
qqnorm(w.AMIind.34000)# not at all normal distr
qqline(w.AMIind.34000)# not at all normal distr
#not norm distr so lm not really good here?

summary(glm(df.w.wide$`40%`~ as.vector(w.AMIind.34000)))


ttt2<-df.w.wide %>% mutate(new_col=c(as.vector(w.AMIind.34000)))
pivot_longer(ttt2, 4:13)
ggplot(pivot_longer(ttt2,4:13), aes(x=new_col, y=value, colour=year))+
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
  facet_wrap(~factor(name, levels=c("binary","10%", "20%","30%","40%",
                                    "50%", "60%", "70%", "80%", "90%")),ncol=5)+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  theme(strip.text.x = element_text(size = 18))+
  ggtitle(element_text("weighted webs normalized AMI & R50  at diffrent THs"))

   
# re-read the paper
# redo AMI check their examples how they did it.

# plot against species number w.con con etc., time CHECK but not sp

# need the maximum uncertainty and divide with AMI with that, start with normalizing!!!! CHECK

# try and calculate A/C instead CHECK

# go back and try random extinctions

# justify our deletion order and connectance as measure
