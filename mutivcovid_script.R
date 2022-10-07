
#### LOAD REQUIRED PACKAGES ####
library(stpca)
library(stringdist)
library(sp)
library(geonames)
library(reshape2)
library(dplyr)
library(expm)
library(rgl)
library(mgcv)
library(ggplot2)


g_un <- "" #Username for geonames

dat11 <- read.csv("latest.txt") #Read in incidence data
# reverse so in order of increasing date
dat11<-dat11[nrow(dat11):1,]
# shorten country
names(dat11)<-ifelse(names(dat11)=="countriesAndTerritories","country",names(dat11))
# shorten some names
levels(dat11$country) <- c(levels(dat11$country), gsub("_", " ", levels(dat11$country)),"United States","Bosnia","North Macedonia")
dat11$country[dat11$country=="United_States_of_America"] <- "United States"
dat11$country <- gsub("_", " ", dat11$country)
dat11$country[dat11$country=="South_Africa"] <- "South Africa"
# dat11$country[dat11$country=="Dominican_Republic"]
# dat11$country[dat11$country=="United_Arab_Emirates"]
dat11$country[dat11$country=="Bosnia_and_Herzegovina"]<-"Bosnia"
dat11$country[dat11$country=="North_Macedonia"] <- "North Macedonia"
dat11$country[dat11$country=="Lao People’s Democratic Republic"] <- "Laos"
# dat11$country[dat11$country=="Democratic_Republic_of_the_Congo"]
dat11$d2020<-difftime(as.POSIXct(strptime(as.Date(dat11$dateRep, format="%d/%m/%Y",tz="UTC"), format="%Y-%m-%d",tz="UTC")), as.POSIXct("2020-01-01",format="%Y-%m-%d", tz="UTC"), units="days")
# there is negative number for the cases on 10 days where data were adjusted
dat11[is.na(dat11$d2020) | dat11$cases<0 | dat11$deaths<0,]
# just set those to na
#dat11<-dat11[!is.na(dat11$d2020) & dat11$cases>=0 & !is.na(dat11$deaths) & dat11$deaths>=0,]
dat11$cases[dat11$cases<0]<-NA_real_
dat11$deaths[dat11$deaths<0]<-NA_real_
# and in Iran the data for the 4th April has been put on the 5th; share those out (it is an odd number)
dat11[dat11$country=="Iran" & dat11$dateRep=="04/04/2020",c("cases","deaths")]<-
  floor(dat11[dat11$country=="Iran" & dat11$dateRep=="05/04/2020",c("cases","deaths")]/2)
dat11[dat11$country=="Iran" & dat11$dateRep=="05/04/2020",c("cases","deaths")]<-
  ceiling(dat11[dat11$country=="Iran" &
                  dat11$dateRep=="05/04/2020",c("cases","deaths")]/2)
# in China a huge pile got added on the 17/04/2020; discard that day
dat11<- dat11[dat11$country!="China"|dat11$dateRep!="17/04/2020",]

# dat11sub <- dat11 %>% 
#   filter(country=="United States"|country=="Brazil"|country=="United Kingdom"|country=="Mexico") %>%
#   select(c(cases,deaths,country,d2020)) %>%
#   group_by(country)
# 
# m <- max(count(dat11sub)$n)
# nmeas <- count(dat11sub)$n
# data.cases <- array(NA,c(m,4))
# data.deaths <- array(NA,c(m,4))
# countries <- c("Brazil","Mexico","United States","United Kingdom")
# for(i in 1:4){
#   data.cases[1:nmeas[i],i] <- c(dat11sub[dat11sub$country==countries[i],"cases"])$cases
#   data.deaths[1:nmeas[i],i] <- c(dat11sub[dat11sub$country==countries[i],"deaths"])$deaths
# }
# 
# data.cov <- array(NA,c(m,4))
# for(i in 1:4){
#   data.cov[1:nmeas[i],i] <- c(dat11sub[dat11sub$country==countries[i],"d2020"])$d2020
# }



options(geonamesUsername=g_un)

apcon <- read.csv("~/Annual-data/prediction.csv")
apinf <- read.csv("~/Annual-data/AirportInfo.csv")



xx1=apply(apinf[1:1000,],1,function(x)try(GNcountryCode(lat=x[6],lng=x[5])$countryName))
xx2=apply(apinf[1001:1488,],1,function(x)try(GNcountryCode(lat=x[6],lng=x[5])$countryName))


xx <- c(xx1,xx2)

xna1 <- grepl("Error", xx1, fixed = TRUE)
xna2 <- grepl("Error", xx2, fixed = TRUE)
xna <- c(xna1,xna2)
xx[xna]<-NA
xxmiss <- sapply(which(is.na(xx)),function(x)try(GNcountryCode(apinf[x,6],lng=apinf[x,5],radius=100)$countryName))
xx[which(is.na(xx))] <- xxmiss


apinf$country <- xx

apcon$countryfrom <- apinf$country[apply(apcon,1,function(x)match(x[1],apinf$NodeName))]
apcon$countryto <- apinf$country[apply(apcon,1,function(x)match(x[2],apinf$NodeName))]

### Normalise by population size
Global_pop <- read.csv("~/OneDrive - University of Glasgow/Mepi/Wen_internship/Global_pop/Global_pop.csv", header=FALSE, comment.char="#")[-c(1,2),]
names(Global_pop) <- Global_pop[1,]
Global_pop <- Global_pop[-1,c("Country Name","2019")]
Global_pop[Global_pop$`Country Name`=="Eritrea",2]<- 3213972
Global_pop[128,1] <- "Laos"
uncs <- unique(dat11$country)


na.omit(unique(xx))[unlist(lapply(1:length(uncs),function(x)which(grepl(uncs[x],na.omit(unique(xx))))))]


predct <- na.omit(unique(Global_pop$`Country Name`))[unlist(lapply(1:length(unique(dat11$country)),function(x){zz<-which(grepl(unique(dat11$country)[x],na.omit(Global_pop$`Country Name`)));ifelse(is.null(zz),NA,zz)}))]
cntrels <- cbind(unique(dat11$country),predct)
for(i in 1:dim(cntrels)[1])if(!is.na(cntrels[i,2])){Global_pop[Global_pop[,1]==cntrels[i,2],1]<-cntrels[i,1]}

bigcnt <- na.omit(cntrels)[,1]

for(i in bigcnt)dat11$cases[dat11$country==i]<-dat11$cases[dat11$country==i]/Global_pop[Global_pop$`Country Name`==i,2]
for(i in bigcnt)dat11$deaths[dat11$country==i]<-dat11$deaths[dat11$country==i]/Global_pop[Global_pop$`Country Name`==i,2]

#### CONSTRUCT SPATIAL WEIGHT MATRIX ######
uncs <- unique(dat11$country)

na.omit(unique(xx))[unlist(lapply(1:length(uncs),function(x)which(grepl(uncs[x],na.omit(unique(xx))))))]


predct <- na.omit(unique(xx))[unlist(lapply(1:length(unique(dat11$country)),function(x){zz<-which(grepl(unique(dat11$country)[x],na.omit(unique(xx))));ifelse(is.null(zz),NA,zz)}))]
predct[69] <- na.omit(unique(xx))[144] #Niger
predct[101] <- "Lao People’s Democratic Republic" #Laos
predct[114] <- na.omit(unique(xx))[96] #Ireland
predct[154] <- NA #Dominica
cntrels <- cbind(unique(dat11$country),predct)

#bigcnt <- c("United States","Brazil","India","Russia","Peru","United Kingdom","Chile","Mexico","South Africa","Iran","Spain","Pakistan","France","Germany","China","Sweden","New Zealand","Australia", "Republic of Korea")
bigcnt2 <- na.omit(cntrels)[,1]
#bigcnt <- bigcnt[bigcnt!="Peru"] #Remove Peru as date 67 and 87 unrealistically
bigcnt <- intersect(bigcnt,bigcnt2) #Countries with pop and mobility data
spatwt = array(0,dim=rep(length(bigcnt),2))

for(k in 1:length(bigcnt)){
  cntk = match(bigcnt[k],cntrels[,1])
  for(l in 1:length(bigcnt)){
    cntl = match(bigcnt[l],cntrels[,1])
    spatwt[l,k] <- spatwt[l,k] + sum(apcon[apcon$countryfrom==cntrels[cntl,2]&apcon$countryto==cntrels[cntk,2],"PredMu"],na.rm=T)
    
  }
  
}

row.names(spatwt) <- bigcnt
colnames(spatwt) <- bigcnt

spatwtnorm <- t(sapply(1:length(bigcnt),function(x)spatwt[x,]/sum(spatwt[x,-x])))
diag(spatwtnorm) <- diag(spatwtnorm)/max(diag(spatwtnorm))
diag(spatwtnorm) <- 1
yy = sqrtm(spatwtnorm)

spatwtfin <- yy+t(yy)
row.names(spatwtfin) <- bigcnt
colnames(spatwtfin) <- bigcnt

dat11sub<-dat11[dat11$country%in%bigcnt,]

newdat <- na.omit(dcast(dat11sub,d2020~country,value.var = "cases"))[,-1]
newdatd <- na.omit(dcast(dat11sub,d2020~country,value.var = "deaths"))[,-1]


#### CONSTRUCT TEMPORAL WEIGHT MATRIX - CASES ######

casesdat <- cbind(as.numeric(row.names(newdat)),newdat)
names(casesdat)[1] <- "date"
nt <- dim(casesdat)[1]
cormc <- NULL
for(j in 1:length(bigcnt)){
  x <- casesdat[,j+1]
  gamfit <- try(gam(x~s(date),data = casesdat,method="REML"))
  if(class(gamfit)[1]=="gam"){
    cormc <- c(cormc,cor(gamfit$residuals[1:(nt-1)],gamfit$residuals[2:nt]))
  }else{
      cormc <- c(cormc,NA)
  }
  
}


IQR(cormc,na.rm = T)



#### CONSTRUCT TEMPORAL WEIGHT MATRIX - DEATHS ######

deathsdat <- cbind(as.numeric(row.names(newdatd)),newdatd)
names(deathsdat)[1] <- "date"

cormd <- NULL
for(j in 1:length(bigcnt)){
  x <- deathsdat[,j+1]
  gamfit <- try(gam(x~s(date),data = deathsdat,method="REML"))
  if(class(gamfit)[1]=="gam"){
    cormd <- c(cormd,cor(gamfit$residuals[1:(nt-1)],gamfit$residuals[2:nt]))
  }else{
    cormd <- c(cormd,NA)
  }
  
}


IQR(cormd,na.rm = T)


tempwtd <- createWeightT(dim(deathsdat)[1], abs(median(cormd,na.rm = T)), corr.form = "AR1")


tempwt <- createWeightT(dim(casesdat)[1], abs(median(cormc,na.rm = T)), corr.form = "AR1")


#### SMode PCA spatial weights ######


smods = stpca(newdat[,bigcnt],pca.mode = "Smode",pca.wt = "spatial",spatial.wt = spatwtfin)
smodu = stpca(newdat[,bigcnt],pca.mode = "Smode")

smodsd = stpca(newdatd[,bigcnt],pca.mode = "Smode",pca.wt = "spatial",spatial.wt = spatwtfin)
smodud = stpca(newdatd[,bigcnt],pca.mode = "Smode")

#### SMode PCA spatio-temporal weights ######

smodst = stpca(newdat[,bigcnt],pca.mode = "Smode",pca.wt =  "spatiotemporal",temporal.wt = tempwt,spatial.wt = spatwtfin)
smoddst = stpca(newdatd[,bigcnt],pca.mode = "Smode",pca.wt =  "spatiotemporal",temporal.wt = tempwtd,spatial.wt = spatwtfin)

#### TMode PCA spatio-temporal weights ######
tmodu = stpca(newdat[,bigcnt],pca.mode = "Tmode",pca.wt =  "unweighted")
tmods = stpca(newdat[,bigcnt],pca.mode = "Tmode",pca.wt =  "spatial",spatial.wt = spatwtfin)
tmodst = stpca(newdat[,bigcnt],pca.mode = "Tmode",pca.wt =  "spatiotemporal",temporal.wt = tempwt,spatial.wt = spatwtfin)
tmoddu = stpca(newdatd[,bigcnt],pca.mode = "Tmode",pca.wt =  "unweighted")
tmodds = stpca(newdatd[,bigcnt],pca.mode = "Tmode",pca.wt =  "spatial",spatial.wt = spatwtfin)
tmoddst = stpca(newdatd[,bigcnt],pca.mode = "Tmode",pca.wt =  "spatiotemporal",temporal.wt = tempwtd,spatial.wt = spatwtfin)



##### RESULTS AND PLOTS ######
#### SPATIAL ######
smods$spatial$loads[,1:4]

for(i in 1:4)print(smods$spatial$loads[order(abs(smods$spatial$loads[,i]),decreasing = T)[1:10],i])

cbind(bigcnt,smodu$unweighted$loads[,1:2])
smodsd$spatial$loads[,1:2]
cbind(bigcnt,smodud$unweighted$loads[,1:2])


par(mfrow=c(1,2))
plot(smodu$unweighted$loads[,1],smods$spatial$loads[,1],pch=16,xlab="Unweighted PC1",ylab="S-mode PC1")
abline(0,1,col=2)
plot(smodu$unweighted$loads[,2],smods$spatial$loads[,2],pch=16,xlab="Unweighted PC2",ylab="S-mode PC")
abline(0,1,col=2)

par(mfrow=c(3,3))

for(i in 1:9){
  plot(as.Date("2020-01-01")+as.numeric(row.names(smods$spatial$scores)),smods$spatial$scores[,i],xlab=paste0("Date"),ylab=paste0("PC",i),pch=16)
}

for(i in 1:9){
  plot(as.Date("2020-01-01")+as.numeric(row.names(smodsd$spatial$scores)),smodsd$spatial$scores[,i],xlab=paste0("Date"),ylab=paste0("PC",i),pch=16)
}



ticks <- as.Date("2020-01-01")+as.numeric(row.names(smods$spatial$scores))
plot3d(x=as.Date("2020-01-01")+as.numeric(row.names(smods$spatial$scores)),y=smods$spatial$scores[,1],z=smods$spatial$scores[,2], size=10, type='p',xaxt = "n", col=1:4, axes = FALSE)
box3d()
axes3d(c("y", "z")) # default axis labels here
axis3d("x", at = ticks, 
       labels = format(ticks, format = "%b")) 


ticks <- as.Date("2020-01-01")+as.numeric(row.names(smodsd$spatial$scores))
plot3d(x=as.Date("2020-01-01")+as.numeric(row.names(smodsd$spatial$scores)),y=smodsd$spatial$scores[,1],z=smodsd$spatial$scores[,2], size=10, type='p',xaxt = "n", col=1:4, axes = FALSE)
box3d()
axes3d(c("y", "z")) # default axis labels here
axis3d("x", at = ticks, 
       labels = format(ticks, format = "%b")) 


#### SPATIOTEMPORAL ######


sapply(1:4,function(j)names(unlist(lapply(1:10,function(i)smodst$spatiotemporal$loads[order(abs(smodst$spatiotemporal$loads[,j]),decreasing = T)[i],j]))))
smoddst$spatiotemporal$loads[order(abs(smoddst$spatiotemporal$loads[,1]),decreasing = T)[1:10],1]


par(mfrow=c(3,3))

for(i in 1:9){
  plot(as.Date("2020-01-01")+as.numeric(row.names(smodst$spatiotemporal$scores)),smodst$spatiotemporal$scores[,i],xlab=paste0("Date"),ylab=paste0("PC",i),pch=16)

}


for(i in 1:9){
  plot(as.Date("2020-01-01")+as.numeric(row.names(smoddst$spatiotemporal$scores)),smoddst$spatiotemporal$scores[,i],xlab=paste0("Date"),ylab=paste0("PC",i),pch=16)
}



ticks <- as.Date("2020-01-01")+as.numeric(row.names(smodst$spatiotemporal$scores))
plot3d(x=as.Date("2020-01-01")+as.numeric(row.names(smodst$spatiotemporal$scores)),y=smodst$spatiotemporal$scores[,1],z=smodst$spatiotemporal$scores[,2], size=10, type='p',xaxt = "n", col=1:4, axes = FALSE)
box3d()
axes3d(c("y", "z")) # default axis labels here
axis3d("x", at = ticks, 
       labels = format(ticks, format = "%b")) 


ticks <- as.Date("2020-01-01")+as.numeric(row.names(smoddst$spatiotemporal$scores))
plot3d(x=as.Date("2020-01-01")+as.numeric(row.names(smoddst$spatiotemporal$scores)),y=smoddst$spatiotemporal$scores[,1],z=smoddst$spatiotemporal$scores[,2], size=10, type='p',xaxt = "n", col=1:4, axes = FALSE)
box3d()
axes3d(c("y", "z")) # default axis labels here
axis3d("x", at = ticks, 
       labels = format(ticks, format = "%b")) 


library(gridExtra)
library(egg)

p1 = ggplot(mapping=aes(x=as.Date("2020-01-01")+as.numeric(row.names(smodu$unweighted$scores)),y=smodu$unweighted$scores[,1]))+ geom_point(size=2, shape=19)+
  ggtitle("Unwweighted SMode") + xlab("Date") + ylab("PC1")
p2 = ggplot(mapping=aes(x=as.Date("2020-01-01")+as.numeric(row.names(smods$spatial$scores)),smods$spatial$scores[,1]))+ geom_point(size=2, shape=19)+
  ggtitle("Spatial SMode") +  xlab("Date") + ylab("PC1")
p3 = ggplot(mapping=aes(x=as.Date("2020-01-01")+as.numeric(row.names(smodst$spatiotemporal$scores)),smodst$spatiotemporal$scores[,1]))+ geom_point(size=2, shape=19)+
  ggtitle("Spatiotemporal SMode") + xlab("Date") + ylab("PC1")
p4 = ggplot(mapping=aes(x=as.Date("2020-01-01")+as.numeric(row.names(smodu$unweighted$scores)),y=smodu$unweighted$scores[,2]))+ geom_point(size=2, shape=19)+
  xlab("Date") + ylab("PC2")
p5 = ggplot(mapping=aes(x=as.Date("2020-01-01")+as.numeric(row.names(smods$spatial$scores)),smods$spatial$scores[,2]))+ geom_point(size=2, shape=19)+
  xlab("Date") + ylab("PC2")
p6 = ggplot(mapping=aes(x=as.Date("2020-01-01")+as.numeric(row.names(smodst$spatiotemporal$scores)),smodst$spatiotemporal$scores[,2]))+ geom_point(size=2, shape=19)+
  xlab("Date") + ylab("PC2")
p7 = ggplot(mapping=aes(x=as.Date("2020-01-01")+as.numeric(row.names(smodu$unweighted$scores)),y=smodu$unweighted$scores[,3]))+ geom_point(size=2, shape=19)+
  xlab("Date") + ylab("PC3")
p8 = ggplot(mapping=aes(x=as.Date("2020-01-01")+as.numeric(row.names(smods$spatial$scores)),smods$spatial$scores[,3]))+ geom_point(size=2, shape=19)+
  xlab("Date") + ylab("PC3")
p9 = ggplot(mapping=aes(x=as.Date("2020-01-01")+as.numeric(row.names(smodst$spatiotemporal$scores)),smodst$spatiotemporal$scores[,3]))+ geom_point(size=2, shape=19)+
  xlab("Date") + ylab("PC3")

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, ncol = 3)

ggsave("PCA_cases.pdf",device = 'pdf')

p1 = ggplot(mapping=aes(x=as.Date("2020-01-01")+as.numeric(row.names(smodud$unweighted$scores)),y=smodud$unweighted$scores[,1]))+ geom_point(size=2, shape=19)+
  ggtitle("Unwweighted SMode") + xlab("Date") + ylab("PC1")
p2 = ggplot(mapping=aes(x=as.Date("2020-01-01")+as.numeric(row.names(smodsd$spatial$scores)),smodsd$spatial$scores[,1]))+ geom_point(size=2, shape=19)+
  ggtitle("Spatial SMode") +  xlab("Date") + ylab("PC1")
p3 = ggplot(mapping=aes(x=as.Date("2020-01-01")+as.numeric(row.names(smoddst$spatiotemporal$scores)),smoddst$spatiotemporal$scores[,1]))+ geom_point(size=2, shape=19)+
  ggtitle("Spatiotemporal SMode") + xlab("Date") + ylab("PC1")
p4 = ggplot(mapping=aes(x=as.Date("2020-01-01")+as.numeric(row.names(smodud$unweighted$scores)),y=smodud$unweighted$scores[,2]))+ geom_point(size=2, shape=19)+
  xlab("Date") + ylab("PC2")
p5 = ggplot(mapping=aes(x=as.Date("2020-01-01")+as.numeric(row.names(smodsd$spatial$scores)),smodsd$spatial$scores[,2]))+ geom_point(size=2, shape=19)+
  xlab("Date") + ylab("PC2")
p6 = ggplot(mapping=aes(x=as.Date("2020-01-01")+as.numeric(row.names(smoddst$spatiotemporal$scores)),smoddst$spatiotemporal$scores[,2]))+ geom_point(size=2, shape=19)+
  xlab("Date") + ylab("PC2")
p7 = ggplot(mapping=aes(x=as.Date("2020-01-01")+as.numeric(row.names(smodud$unweighted$scores)),y=smodud$unweighted$scores[,3]))+ geom_point(size=2, shape=19)+
  xlab("Date") + ylab("PC3")
p8 = ggplot(mapping=aes(x=as.Date("2020-01-01")+as.numeric(row.names(smodsd$spatial$scores)),smodsd$spatial$scores[,3]))+ geom_point(size=2, shape=19)+
  xlab("Date") + ylab("PC3")
p9 = ggplot(mapping=aes(x=as.Date("2020-01-01")+as.numeric(row.names(smoddst$spatiotemporal$scores)),smoddst$spatiotemporal$scores[,3]))+ geom_point(size=2, shape=19)+
  xlab("Date") + ylab("PC3")

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, ncol = 3)

ggsave("PCA_deaths.pdf",device = 'pdf')




p1 = ggplot(mapping=aes(x=names(newdat),y=tmodu$unweighted$scores[,1]))+ geom_point(size=2, shape=19)+
  ggtitle("Unwweighted Tmode") + xlab("Country") + ylab("PC1")
p2 = ggplot(mapping=aes(x=names(newdat),y=tmods$spatial$scores[,1]))+ geom_point(size=2, shape=19)+
  ggtitle("Spatial Tmode") +  xlab("Country") + ylab("PC1")
p3 = ggplot(mapping=aes(x=names(newdat),y=tmodst$spatiotemporal$scores[,1]))+ geom_point(size=2, shape=19)+
  ggtitle("Spatiotemporal Tmode") + xlab("Country") + ylab("PC1")
p4 = ggplot(mapping=aes(x=names(newdat),y=tmodu$unweighted$scores[,2]))+ geom_point(size=2, shape=19)+
  xlab("Country") + ylab("PC2")
p5 = ggplot(mapping=aes(x=names(newdat),y=tmods$spatial$scores[,2]))+ geom_point(size=2, shape=19)+
  xlab("Country") + ylab("PC2")
p6 = ggplot(mapping=aes(x=names(newdat),tmodst$spatiotemporal$scores[,2]))+ geom_point(size=2, shape=19)+
  xlab("Country") + ylab("PC2")
p7 = ggplot(mapping=aes(x=names(newdat),y=tmodu$unweighted$scores[,3]))+ geom_point(size=2, shape=19)+
  xlab("Country") + ylab("PC3")
p8 = ggplot(mapping=aes(x=names(newdat),tmods$spatial$scores[,3]))+ geom_point(size=2, shape=19)+
  xlab("Country") + ylab("PC3")
p9 = ggplot(mapping=aes(x=names(newdat),tmodst$spatiotemporal$scores[,3]))+ geom_point(size=2, shape=19)+
  xlab("Country") + ylab("PC3")

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, ncol = 3)

ggsave("PCA_cases_tmode.pdf",device = 'pdf')




p1 = ggplot(mapping=aes(x=names(newdatd),y=tmoddu$unweighted$scores[,1]))+ geom_point(size=2, shape=19)+
  ggtitle("Unwweighted Tmode") + xlab("Country") + ylab("PC1")
p2 = ggplot(mapping=aes(x=names(newdatd),y=tmodds$spatial$scores[,1]))+ geom_point(size=2, shape=19)+
  ggtitle("Spatial Tmode") +  xlab("Country") + ylab("PC1")
p3 = ggplot(mapping=aes(x=names(newdatd),y=tmoddst$spatiotemporal$scores[,1]))+ geom_point(size=2, shape=19)+
  ggtitle("Spatiotemporal Tmode") + xlab("Country") + ylab("PC1")
p4 = ggplot(mapping=aes(x=names(newdatd),y=tmoddu$unweighted$scores[,2]))+ geom_point(size=2, shape=19)+
  xlab("Country") + ylab("PC2")
p5 = ggplot(mapping=aes(x=names(newdatd),y=tmodds$spatial$scores[,2]))+ geom_point(size=2, shape=19)+
  xlab("Country") + ylab("PC2")
p6 = ggplot(mapping=aes(x=names(newdatd),tmoddst$spatiotemporal$scores[,2]))+ geom_point(size=2, shape=19)+
  xlab("Country") + ylab("PC2")
p7 = ggplot(mapping=aes(x=names(newdatd),y=tmoddu$unweighted$scores[,3]))+ geom_point(size=2, shape=19)+
  xlab("Country") + ylab("PC3")
p8 = ggplot(mapping=aes(x=names(newdatd),tmodds$spatial$scores[,3]))+ geom_point(size=2, shape=19)+
  xlab("Country") + ylab("PC3")
p9 = ggplot(mapping=aes(x=names(newdatd),tmoddst$spatiotemporal$scores[,3]))+ geom_point(size=2, shape=19)+
  xlab("Country") + ylab("PC3")

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, ncol = 3)

ggsave("PCA_deaths_tmode.pdf",device = 'pdf')



############
#CLUSTERING#
############

library(dtwclust)
library(dplyr)

caselst <- lapply(unique(dat11$country),function(x)dat11$cases[dat11$country==x])
deathlst <- lapply(unique(dat11$country),function(x)dat11$deaths[dat11$country==x])

datelst <- lapply(unique(dat11$country),function(x)dat11$d2020[dat11$country==x])

names(caselst) <- unique(dat11$country)
names(deathlst) <- unique(dat11$country)
names(datelst) <- unique(dat11$country)


acf_fun <- function(series, ...) {
  lapply(series, function(x) {
    y=as.numeric(acf(na.omit(x), lag.max = 5, plot = FALSE,na.action = na.omit)$acf)#;ifelse(is.na(y)|is.na(y),1,y)
  })
  
}
# Fuzzy c-means
fc <- tsclust(tslist(caselst), type = "f", k = 6L,
              preproc = acf_fun, distance = "L2",
              seed = 42)

plot(fc)

lapply(1:6,function(i)names(which(apply(fc@fcluster,1,which.max)==i)))
xx<-data.frame(lapply(lapply(1:6,function(i)names(which(apply(fc@fcluster,1,which.max)==i))), "length<-", max(lengths(lapply(1:6,function(i)names(which(apply(fc@fcluster,1,which.max)==i)))))))
names(xx)<-sapply(1:6,function(i)paste0("Cluster_",i))

table(apply(fc@fcluster,1,which.max))

##############
#            #
# CLUSTERING #
#            #
##############


newdatl <- lapply(seq_len(ncol(newdat)), function(i) newdat[,i])

clusta_12<-tsclust(newdatl,type="partitional",k=12L,distance = "dtw",clustering = "pam")
plot(clusta_12,type="series")
plot(clusta_12,type="sc")
plot(clusta_12,type="centroids")
plot(clusta_12,type="centroids",clus=1L)

#hierachical


hier.12<-tsclust(t(newdat),type = "h",k=12L,distance = "dtw")
plot(hier.12,type="series")
plot(hier.12)
cvi(hier.12)
cvi(clusta_12)
cls <- data.frame(cutree(hier.12,k=6L))
row.names((cls))[cls==6]
table(cutree(hier.12,k=6L))
hier.12$labels

lapply(1:6,function(i)names(which(cls==i)))

par(mfrow=c(2,1))

plot(rev(complete_data[,15]))
plot(rev(complete_data[,16]))


acf_fun <- function(series, ...) {
  lapply(series, function(x) {
    y=as.numeric(acf(x, lag.max = 50, plot = FALSE,na.action = na.omit)$acf);ifelse(is.na(y),1,y)
  })
  
}
# Fuzzy c-means
fc <- tsclust(t(newdat), type = "f", k = 6L,
              preproc = acf_fun, distance = "L2",
              seed = 42)

plot(fc)

lapply(1:6,function(i)names(which(apply(fc@fcluster,1,which.max)==i)))

knitr::kable(xx)

table(apply(fc@fcluster,1,which.max))