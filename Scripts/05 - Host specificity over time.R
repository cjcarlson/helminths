
library(pspline)

genera <- read.csv('~/Github/helminths/Data/timeseries_NHM_v2.csv',
                   stringsAsFactors = FALSE)[,-1]



genera.sub <- genera[!genera$first=='exclude',]
#boxplot(log(hs) ~ first, data=genera.sub)
summary(aov(hs~as.factor(first)+Error(as.factor(genus)), data = genera.sub))
kruskal.test(hs~as.factor(first), data = genera.sub)
wilcox.test(hs~as.factor(first), data=genera.sub)


# genera

gen.omit <- na.omit(genera.sub)
gen.omit <- gen.omit[!(gen.omit$hs==0),]
gen.omit$loghs <- log(gen.omit$hs)


summary(aov(hs~as.factor(first)+Error(as.factor(genus)), data = gen.omit))
kruskal.test(hs~as.factor(first), data = gen.omit)
wilcox.test(hs~as.factor(first), data=gen.omit)

#plot(gen.omit$loghs ~ gen.omit$year, pch=16, xlab='year', ylab='log(host specificity)')
#lines(smooth.spline(gen.omit$year,gen.omit$loghs, spar=0.5),type='l',col='red',lwd=2)

MyGAM1 <- gam(hs ~  s(year), family=nb, data=gen.omit)
summary(MyGAM1)
#plot(MyGAM1, residuals=F, se=TRUE, pch=19, cex=0.75, scheme=0, shade=F, col='red', rug=TRUE,
#     xlab = "year", ylab = "log(host specificity)")


yofa <- read.csv('C:/Users/cjcar/Dropbox/helminths/YOFA hosts.csv')
MyGAM0 <- gam(hosts ~ s(yearone), family=nb, data=yofa)

#lines(spec, residuals=F, se=TRUE,pch=19, cex=0.75, scheme=0, shade=F, col='blue', xlab = "year", ylab = "log(host specificity)", add=TRUE)

response1 <- predict(MyGAM1, type="response", se.fit=T)
response0 <- predict(MyGAM0, type="response", se.fit=T)


layout(matrix(c(1,2,2,3), 1, 3, byrow = TRUE), 
       widths=c(1,1,1,1), heights=c(1,2))


par(mar=c(5,5,2,3))
boxplot(log(hs) ~ first, data=genera.sub, ylab='log(host specificity)', xlab='', bty="n", 
        col=c('#e9a3c9','#a1d76a'), cex.lab=1.4, cex.axis=1.4)

plot(0, type="n", bty="n", main="", xlab="Year", ylab="host specificity", lwd=3, ylim=c(0,50), xlim=c(1700,2045),
     cex.lab=1.4, cex.axis=1.4)
legend("topright", bty="n", lwd=3, col=c("#998ec3","#f1a340"), legend=c("NHM", "USNPC"),
       cex=1.4)

lines(sm.spline(MyGAM1$model$year , response1$fit) , lwd = 3 , col = "#998ec3")
lines(sm.spline(MyGAM1$model$year , response1$fit+1.96*response1$se) , lty = 2 , lwd = 2 , col = "#998ec3")
lines(sm.spline(MyGAM1$model$year , response1$fit-1.96*response1$se) , lty = 2 , lwd = 2 , col = "#998ec3")

lines(sm.spline(MyGAM0$model$yearone , response0$fit) , lwd = 3 , col = "#f1a340")
lines(sm.spline(MyGAM0$model$yearone , response0$fit + 1.96 * response0$se) , lty = 2 , lwd = 2, col = "#f1a340")
lines(sm.spline(MyGAM0$model$yearone , response0$fit - 1.96 * response0$se) , lty = 2 , lwd = 2 , col = "#f1a340")


