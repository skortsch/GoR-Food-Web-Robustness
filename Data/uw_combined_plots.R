#uw combined plot
uwR50th90<-readRDS("uw.R50th90.RDS")
uwR50th80<-readRDS("uw.R50th80.RDS")
uwR50th70<-readRDS("uw.R50th70.RDS")
uwR50th60<-readRDS("uw.R50th60.RDS")
uwR50th50<-readRDS("uw.R50th50.RDS")
uwR50th40<-readRDS("uw.R50th40.RDS")
uwR50th30<-readRDS("uw.R50th30.RDS")
uwR50th20<-readRDS("uw.R50th20.RDS")
uwR50th10<-readRDS("uw.R50th10.RDS")
uwR50Binary<-readRDS("uw.R50.output.1000.rds")
uwR50.list<-list(uwR50Binary, uwR50th10, uwR50th20, uwR50th30, uwR50th40,
                 uwR50th50, uwR50th60, uwR50th70, uwR50th80, uwR50th90)
saveRDS(uwR50.list, file = "unweighted_R50_list.RDS")
uwR50.list

df.allTH<-matrix(NA,34,11)
df.allTH[,11]<-c(seq(1981,2014))
colnames(df.allTH)<-c("Binary","TH 10%","TH 20%", "TH 30%", "TH 40%", "TH 50%",
                      "TH 60%", "TH 70%", "TH 80%", "TH 90%", "Year")
for (i in 1:10){
  df.allTH[,i]<-apply(uwR50.list[[i]],1, median)
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
  ggtitle("Unweighted food webs R50 with Thresholds")+theme(plot.title = element_text(hjust = 0.5, size = 25))+
  labs(color = "Thresholds")+
  theme(legend.text=element_text(size=25),
        legend.title=element_text(size=25))



uwR50THplots<-list()
for (i in 1:10){
  data<-data.frame(apply(uwR50.list[[i]],1,median), apply(Con,1,median))
  colnames(data)<-c("R50", "Con")
  uwR50THplots[i]<-list(ggplot(data, aes(Con, R50))+
                          geom_point(alpha=1)+
                          stat_smooth(method = "lm",
                                      formula = "y~x",
                                      geom = "smooth",
                                      col="black")+
                          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                panel.background = element_blank(), axis.line = element_line(colour = "black"))+
                          theme(axis.title.x = element_text(size=20))+
                          theme(axis.title.y = element_text(size=20))+
                          xlab("Connectance"))
}
uwR50THplots[[1]]+ggtitle("TH Binary")+uwR50THplots[[2]]+ggtitle("TH 10%")+uwR50THplots[[3]]+
  ggtitle("TH 20%")+uwR50THplots[[4]]+ggtitle("TH 30%")+uwR50THplots[[5]]+ggtitle("TH 40%")+
  uwR50THplots[[6]]+ggtitle("TH 50%")+uwR50THplots[[7]]+ggtitle("TH 60%")+
  uwR50THplots[[8]]+ggtitle("TH 70%")+uwR50THplots[[9]]+ggtitle("TH 80%")+
  uwR50THplots[[10]]+ggtitle("TH 90%")



# combined plot Con vs R50, tot data not only median
uwR50THplots<-list()
for (i in 1:10){
  data<-data.frame(as.vector(uwR50.list[[i]]), as.vector(Con))
  colnames(data)<-c("R50", "Con")
  uwR50THplots[i]<-list(ggplot(data, aes(Con, R50))+
                          geom_point(alpha=1)+
                          stat_smooth(method = "lm",
                                      formula = "y~x",
                                      geom = "smooth",
                                      col="Red")+
                          scale_x_continuous(labels = scales::number_format(accuracy = 0.01), 
                                             breaks = scales::pretty_breaks(n = 5))+
                          theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 15))+
                          theme(axis.text.y = element_text(vjust = 0.5, size = 15))+
                          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                panel.background = element_blank(), axis.line = element_line(colour = "black"))+
                          xlab("Connectance")+ylab("R50")+
                          theme(axis.title=element_text(size=15)))
                          
}
uwR50THplots[[1]]+ggtitle("TH Binary")+uwR50THplots[[2]]+ggtitle("TH 10%")+uwR50THplots[[3]]+
  ggtitle("TH 20%")+uwR50THplots[[4]]+ggtitle("TH 30%")+uwR50THplots[[5]]+ggtitle("TH 40%")+
  uwR50THplots[[6]]+ggtitle("TH 50%")+uwR50THplots[[7]]+ggtitle("TH 60%")+
  uwR50THplots[[8]]+ggtitle("TH 70%")+uwR50THplots[[9]]+ggtitle("TH 80%")+
  uwR50THplots[[10]]+ggtitle("TH 90%")
