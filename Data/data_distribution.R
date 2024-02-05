library(ggplot2)
library(tidyverse)

df.uw.long<-readRDS("df.unweighted.RDS") #unweighted webs df, long format for plot
df.w.long<-readRDS("df.weighted.RDS")# weighted webs df, longformat for plot

df.uw.wide<-readRDS("df.uw.wide.RDS")# unweighted webs df, wide format,
df.w.wide<-readRDS("df.w.wide.RDS")# weighted webs df, wide format
#Con= connectance,binary+10%-90% values = R50 values for the thresholds
#w.con= weighted connectance, binary+10%-90% values = R50 values for the thresholds

# lm plot at diffrent thresholds weighted webs, use "df.uw.long" for unweighted webs
ggplot(df.w.long, aes(x=w.con, y=R50, colour=year))+
  geom_point()+
  stat_smooth(method = "glm", method.args=list(family="gaussian"),
              formula = "y~x",
              geom = "smooth",
              col="Red")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 18))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  scale_color_gradient(low="yellow", high="Blue")+
  theme(plot.title = element_text(hjust = 0))+ylab("R50")+
  theme(strip.text.x = element_text(size = 15))+
  theme(axis.title=element_text(size=28))+theme(legend.text=element_text(size=25),
                                                legend.title=element_text(size=25))+
  facet_wrap(~factor(threshold, levels=c("binary","10%", "20%","30%","40%",
                                         "50%", "60%", "70%", "80%", "90%")),ncol=5)+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  theme(strip.text.x = element_text(size = 18))

#lm to get R^2, R^2 is how well the model fits the data, right?
#value from 0-1, 0= bad fit, 1= prefect fit, right?
df.w.wide
summary(lm(df.w.wide$w.con~df.w.wide$binary)) # lowest R^2, 0.001362
summary(lm(df.w.wide$w.con~df.w.wide$`10%`)) # 0.106
summary(lm(df.w.wide$w.con~df.w.wide$`20%`)) #0.2935
summary(lm(df.w.wide$w.con~df.w.wide$`30%`)) # 0.3953
summary(lm(df.w.wide$w.con~df.w.wide$`40%`)) # "best" R^2, 0.4291
summary(lm(df.w.wide$w.con~df.w.wide$`50%`)) #0.4119 
summary(lm(df.w.wide$w.con~df.w.wide$`60%`)) #0.3778
summary(lm(df.w.wide$w.con~df.w.wide$`70%`)) #0.3489
summary(lm(df.w.wide$w.con~df.w.wide$`80%`)) #0.2169
summary(lm(df.w.wide$w.con~df.w.wide$`90%`)) #0.1806



# checking data distribution of R50, normal distr or not
shapiro.test(df.w.wide$binary) # error: sample size must be between 3 and 5000
#any other official test to run?

ks.test(df.w.wide$binary, "pnorm")# tested all thresholds, all not norm.distr? Since the p-value is less than .05, we reject the null hypothesis
ks.test(df.uw.wide$binary, "pnorm")
#checking data distribution visually

# histogram plot of R50
ggplot(df.w.long, aes(x=R50))+#use df.uw.long to se for unweighted
  geom_histogram()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  facet_wrap(~factor(threshold,levels=c("binary","10%", "20%","30%","40%",
                                        "50%", "60%", "70%", "80%", "90%")),ncol=5)+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  ggtitle(element_text("R50 histplot at diffrent THs"))

ggplot(df.w.long, aes(sample=R50))+
  stat_qq()+
  stat_qq_line()+
  facet_wrap(~factor(threshold,levels=c("binary","10%", "20%","30%","40%",
                                        "50%", "60%", "70%", "80%", "90%")),ncol=5)+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  ggtitle(element_text("R50 qqplot at diffrent THs"))

#unweigheted a bit more normaly distributed for some Ths
#not really gaussian distribution for uw or weighted R50 at the diffrent thresholds, but is it close enough?
#can we still use the lm? or do we need some other model like glm or?



# is this useful?
#checking assumption of Normality of residuals
#Normality of the residuals is an assumption of running a linear model. 
#So, if your residuals are normal, it means that your assumption is 
#valid and model inference (confidence intervals, model predictions) should also be valid?
par(mfrow=c(2,2))
lmbin<-lm(df.w.wide$w.con~df.w.wide$binary)
plot(lmbin)
lm10<-lm(df.w.wide$w.con~df.w.wide$`10%`)
plot(lm10)
lm20<-lm(df.w.wide$w.con~df.w.wide$`20%`)
plot(lm20)
lm30<-lm(df.w.wide$w.con~df.w.wide$`30%`)
plot(lm30)
lm40<-lm(df.w.wide$w.con~df.w.wide$`40%`)
plot(lm40)
lm50<-lm(df.w.wide$w.con~df.w.wide$`50%`)
plot(lm50)
lm60<-lm(df.w.wide$w.con~df.w.wide$`60%`)
plot(lm60)
lm70<-lm(df.w.wide$w.con~df.w.wide$`70%`)
plot(lm70)
lm80<-lm(df.w.wide$w.con~df.w.wide$`80%`)
plot(lm80)
lm90<-lm(df.w.wide$w.con~df.w.wide$`90%`)
plot(lm90)


##qqplot w.con & hist 
ggplot(df.w.wide, aes(sample=w.con))+
  stat_qq()+
  stat_qq_line()

ggplot(df.w.wide, aes(w.con))+
  geom_histogram()
#not really gaussian distribution for w.con

#transformation ?
ggplot(df.w.long, aes(sample=log10(R50)))+
  stat_qq()+
  stat_qq_line()+
  facet_wrap(~factor(threshold,levels=c("binary","10%", "20%","30%","40%",
                                        "50%", "60%", "70%", "80%", "90%")),ncol=5)+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  ggtitle(element_text("R50 qqplot at diffrent THs log10 transformation"))

#looks better now, except for the 90% scenario

summary(lm(df.w.wide$w.con~log10(df.w.wide$binary))) # 0.001315 vs 0.001362 from non log transformed model
summary(lm(df.w.wide$w.con~log10(df.w.wide$`10%`)))# 0.09744 vs 0.106
summary(lm(df.w.wide$w.con~log10(df.w.wide$`20%`))) #0.2886 vs 0.2935
summary(lm(df.w.wide$w.con~log10(df.w.wide$`30%`))) # .406 vs 0.3953
summary(lm(df.w.wide$w.con~log10(df.w.wide$`40%`))) # "0.4503 vs 0.4291
summary(lm(df.w.wide$w.con~log10(df.w.wide$`50%`)))# 0.4221 vs 0.4119 
summary(lm(df.w.wide$w.con~log10(df.w.wide$`60%`)))#0.3693 vs 0.3778
summary(lm(df.w.wide$w.con~log10(df.w.wide$`70%`)))#0.3274 vs 0.3489
summary(lm(df.w.wide$w.con~log10(df.w.wide$`80%`))) #0.2054 vs 0.2169
summary(lm(df.w.wide$w.con~log10(df.w.wide$`90%`))) #0.1666 vs 0.1806


### glm attempt
#"For running the linear models, 
#I would use a glm with identity link for the response variable with normal dist.
#and glm with log link for the skewed ones." from Susanne
summary(glm(df.w.wide$w.con~df.w.wide$binary, family = gaussian(link = "log")))
summary(glm(df.w.wide$w.con~df.w.wide$`10%`, family = gaussian(link = "log")))
summary(glm(df.w.wide$w.con~df.w.wide$`20%`, family = gaussian(link = "log")))
summary(glm(df.w.wide$w.con~df.w.wide$`30%`, family = gaussian(link = "log")))
summary(glm(df.w.wide$w.con~df.w.wide$`40%`, family = gaussian(link = "identity")))
summary(glm(df.w.wide$w.con~df.w.wide$`50%`, family = gaussian(link = "identity")))
summary(glm(df.w.wide$w.con~df.w.wide$`60%`, family = gaussian(link = "identity")))
summary(glm(df.w.wide$w.con~df.w.wide$`70%`, family = gaussian(link = "log")))
summary(glm(df.w.wide$w.con~df.w.wide$`80%`, family = gaussian(link = "log")))
summary(glm(df.w.wide$w.con~(df.w.wide$`90%`), family =gaussian(link = "log")))



summary(glm(df.w.wide$w.con~log10(df.w.wide$binary), family = gaussian(link = "identity")))
summary(glm(df.w.wide$w.con~log10(df.w.wide$`10%`), family = gaussian(link = "identity")))
summary(glm(df.w.wide$w.con~log10(df.w.wide$`20%`), family = gaussian(link = "identity")))
summary(glm(df.w.wide$w.con~log10(df.w.wide$`30%`), family = gaussian(link = "identity")))
summary(glm(df.w.wide$w.con~log10(df.w.wide$`40%`), family = gaussian(link = "identity")))
summary(glm(df.w.wide$w.con~log10(df.w.wide$`50%`), family = gaussian(link = "identity")))
summary(glm(df.w.wide$w.con~log10(df.w.wide$`60%`), family = gaussian(link = "identity")))
summary(glm(df.w.wide$w.con~log10(df.w.wide$`70%`), family = gaussian(link = "identity")))
summary(glm(df.w.wide$w.con~df.w.wide$`80%`, family = gaussian(link = "log")))
summary(glm(df.w.wide$w.con~(df.w.wide$`90%`), family =gaussian(link = "log")))

####correct order for glm####
summary(glm(df.w.wide$binary~df.w.wide$w.con, family = gaussian(link = "identity")))
summary(glm(df.w.wide$`10%`~df.w.wide$w.con, family = gaussian(link = "identity")))
summary(glm(df.w.wide$`20%`~df.w.wide$w.con, family = gaussian(link = "identity")))
summary(glm(df.w.wide$`30%`~df.w.wide$w.con, family = gaussian(link = "identity")))
summary(glm(df.w.wide$`40%`~df.w.wide$w.con, family = gaussian(link = "identity")))
summary(glm(df.w.wide$`50%`~df.w.wide$w.con, family = gaussian(link = "identity")))
summary(glm(df.w.wide$`60%`~df.w.wide$w.con, family = gaussian(link = "identity")))
summary(glm(df.w.wide$`70%`~df.w.wide$w.con, family = gaussian(link = "log")))
summary(glm(df.w.wide$`80%`~df.w.wide$w.con, family = gaussian(link = "log")))
summary(glm(df.w.wide$`90%`~df.w.wide$w.con, family =gaussian(link = "log")))


# put p-values in captions on plots
# draft by end of november to Lai