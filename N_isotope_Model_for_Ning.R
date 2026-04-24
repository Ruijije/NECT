#################################################################N istope model
pop.mod2 <-nls(d15N.14N~ c+a/(1+b*(fR/fTemp)),start=list(a =7.04, b =2.68,c=-2.47), data=nitrogendata, trace=T)
summary(pop.mod2) 
###This model you can get the value of k and εgas
nitrogenisotope_model<--7.6265 +(10.4682/(1+1.1405*(fR/fTemp)))
fgas<-1/(1+1.1405*(fR/fTemp)) 
#####################################
fgasdata<-data.frame(cbind(fgas,ndata$MAP))
fgas.sd<-aggregate(fgasdata$fgas,by=list(fgasdata$V2),FUN=sd)
fgas_mean<-aggregate(fgasdata$fgas,by=list(fgasdata$V2),FUN=mean)
relationshipdata1<-data.frame(cbind(fgas_mean,fgas.sd$x))
relationshipdata1[relationshipdata1<0]<-0
p12<-ggplot(data=datamodel, aes(x=MAP, y=d15N.14N)) + geom_point() + geom_smooth(method = 'lm', formula = y ~ x)+
  ylab(expression('Plant'*δ^15*'N')) + xlab(expression('MAP'*~'(mm)'))+
  geom_point(aes(x=MAP,y=nitrogenisotope_model),col="red") + 
  geom_smooth(aes(y=nitrogenisotope_model,x=MAP, col='red'),method = 'lm') +
  theme_bw()+
  theme( legend.position="none",axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
         axis.text=element_text(size=rel(1.8)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))

mfas2<-ggplot(relationshipdata1, aes(x=Group.1, y=x)) + geom_point() + 
  geom_smooth(method = 'loess', formula = y ~ x, aes(color='red')) + 
  geom_errorbar(aes(ymax=x+fgas.sd$x,ymin=x-fgas.sd$x),position=position_dodge(0.9),width=1.0)+
  ylab(expression(f[gas])) +xlab(expression('MAP'*'('*'mm'*')'))+
  theme_bw()+theme(legend.position='none',axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
                   axis.text=element_text(size=rel(1.8)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))

lmnisotopef<-lm(predictedN~fR+ftemp,data=datamodel)
summary(lmnisotopef)
library(car)
crPlots(lmnisotopef,ylab = expression('Plant'*δ^15*'N'))
