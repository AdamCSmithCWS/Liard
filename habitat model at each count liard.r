### script to run BBS-style model on Liard data
#####
# install.packages("R2admb")
# install.packages("glmmADMB", 
#                  repos=c("http://glmmadmb.r-forge.r-project.org/repos",
#                          getOption("repos")),
#                  type="source")

setwd("c:/Liard")

dir.data = paste0(getwd(),"/data")
dir.model = paste0(getwd(),"/model")

Y = 2017 
library(rjags)
library(jagsUI)
library(loo)
library(glmmADMB)
library(MCMCvis)
library(mgcv)
library(ggplot2)
library(ggrepel)


################ habitat data cleaning and interpolation
# generate linear and squared predictors for
### tree height
### canopy cover
### shrub density
### % conifer

veg = read.csv(paste0(dir.data,"/Master habitat.csv"),
                  stringsAsFactors = F)

############ interpolate veg data
veg$ss.uni = paste(veg$StandNumber,veg$Station,sep = "_")
fyear = tapply(veg$YEAR,veg$ss.uni,min)
lyear = tapply(veg$YEAR,veg$ss.uni,max)

veg.int = expand.grid(ss.uni = unique(veg$ss.uni),
                      YEAR = c(1998:2017))
vpreds = names(veg)[c(6:13)]
for(si in unique(veg$ss.uni)){
  tmp = veg[which(veg$ss.uni == si),]
  if(nrow(tmp)== 1){next}
  for(v in vpreds){
    res = tmp[,v]
    y = tmp$YEAR
    for(j in 2:length(res)){
      y1 = y[j-1]
      y2 = y[j]
      res1 = res[j-1]
      res2 = res[j]
      dres = res2-res1
      dresy = dres/(y2-y1)
      
      if(res1 == res2 | any(is.na(c(res1,res2)))){
        resint = rep(res1,times = length(y1:y2))
      }else{
        resint = seq(from = res1,to = res2,by = dresy)
      }
      
      names(resint) <- c(y1:y2)
      for(yy in c(y1:y2)){
        veg.int[which(veg.int$ss.uni == si & veg.int$YEAR == yy),v] = resint[as.character(yy)]
      }
    }
    
    
  }
}



for(v in vpreds){
  ncl = paste0(v,"_2")
  veg.int[,ncl] = veg.int[,v]^2
  
}

## next line assumes that the habitat measured in 1999 can be 
## substituted in to the same sites in 1998
for(s in names(fyear)[which(fyear == 1999)]){
  veg.int[which(veg.int$ss.uni == s & veg.int$YEAR == 1998),c(3:ncol(veg.int))] = veg.int[which(veg.int$ss.uni == s & veg.int$YEAR == 1999),c(3:ncol(veg.int))]
}
veg.int$ss.uni = as.character(veg.int$ss.uni)


#simplified table defining the species-specific models (global models only)
models = read.csv(paste0(dir.model,"/Model descriptions.csv"),stringsAsFactors = F)

#table linking the predictor numbers with their names in the habitat file
predlist = read.csv(paste0(dir.model,"/predictor list.csv"),stringsAsFactors = F)

mods = models[,c("English_Name","SpeciesID","Proposed.hab.predictors.NUM")]


pretest = F### default of TRUE runs a preliminary (non Bayes) regression to 
### remove non-signficant predictors:
### change pretest to FALSE to retain full, species level model 


for(sp in c("MOWA","PAWA")){  # 
 
  
  ##### if we want to apply the budworm to other species, this next line can be modified
  ##### if we want to apply the budworm to other species, this next line can be modified
  ##### if we want to apply the budworm to other species, this next line can be modified
  ##### if we want to apply the budworm to other species, this next line can be modified
  ##### if we want to apply the budworm to other species, this next line can be modified
if(sp %in% c("TEWA","CMWA","BBWA")){addbud = T}else{addbud = F}


wpreds = mods[which(mods$SpeciesID == sp),3]  
wpreds = as.integer(unlist(strsplit(x = wpreds,split = "[[:punct:]]")))

sppreds = vpreds[which(vpreds %in% predlist[which(predlist$habN %in% wpreds),"habitat"])]  



if(addbud){
budw = read.csv(paste0(dir.data,"/SBW_Master.csv"),stringsAsFactors = F)



budw$ss.uni = paste(budw$StandNumber,budw$Station,sep = "_")

budw = unique(budw[,c("ss.uni","YEAR","NEAR_DIST","SEVERITY")]) #removes the "round" replication

budb = expand.grid(ss.uni = unique(budw$ss.uni),
                   YEAR = c(min(veg.int$YEAR):max(veg.int$YEAR)))

budw = merge(budb,budw,by = c("ss.uni","YEAR"),all.x = T)

cats = c(-1,0,1000,10000,50000,1000000)
budw$budcat = cut(budw$NEAR_DIST,
                  breaks = cats,
                  ordered_result = T)
budw$budw = as.integer(budw$budcat)

yrswbud = unique(budw[which(!is.na(budw$budw)),"YEAR"])

### add reasonable approximate data for budworm in teh missing years

### alternatively, drop 2000 data?


veg.int = merge(veg.int,budw[,c("ss.uni","YEAR","budw")],
                by = c("ss.uni","YEAR"),
                all.x = T)


 
  veg.int = veg.int[which(veg.int$YEAR %in% yrswbud),]
 
vcols = c(sppreds,"budw")
}else{
  vcols = c(sppreds)
  
}


veg.sp = veg.int[,c("ss.uni","YEAR",vcols)]

Npreds = 500
veg.predict = data.frame(step = 1:Npreds)

#veg.sp.p = veg.int[which(veg.int$ss.uni %in% c(predsites))]



for(v in sppreds){

  #ncl = paste0(v,"_2")
  ncls = paste0(v,"_2_s")
  nvs = paste0(v,"_s")
  veg.sp[,nvs] = (veg.sp[,v]-mean(veg.sp[,v],na.rm = T))/sd(veg.sp[,v],na.rm = T)
  veg.sp[,ncls] = veg.sp[,nvs]^2
  veg.predict[,v] = seq(min(veg.sp[,v],na.rm = T),max(veg.sp[,v],na.rm = T),length.out = Npreds)
  #veg.predict[,ncl] = veg.predict[,v]^2
  veg.predict[,nvs] = (veg.predict[,v]-mean(veg.sp[,v],na.rm = T))/sd(veg.sp[,v],na.rm = T)
  veg.predict[,ncls] = veg.predict[,nvs]^2
  
  #   veg.sp[,nvs] = (veg.sp[,v]-mean(veg.sp[,v],na.rm = T))/sd(veg.sp[,v],na.rm = T)
  # veg.sp[,ncls] = (veg.sp[,ncl]-mean(veg.sp[,v],na.rm = T))/sd(veg.sp[,ncl],na.rm = T)
  # veg.predict[,v] = seq(min(veg.sp[,v],na.rm = T),max(veg.sp[,v],na.rm = T),length.out = Npreds)
  # veg.predict[,ncl] = veg.predict[,v]^2
  # veg.predict[,nvs] = (veg.predict[,v]-mean(veg.sp[,v],na.rm = T))/sd(veg.sp[,v],na.rm = T)
  # veg.predict[,ncls] = (veg.predict[,ncl]-mean(veg.sp[,ncl],na.rm = T))/sd(veg.sp[,ncl],na.rm = T)
  # 
  
}



dat = read.csv(paste0(dir.data,"/",sp,".csv"),
               stringsAsFactors = F)

names(dat)[which(names(dat) == "Abundance.x")] = "count"

dat$ss.uni = paste(dat$StandNumber,dat$Station,sep = "_")


dmerg = merge(dat,veg.sp,
              by = c("ss.uni","YEAR"),
              all.x = T)


dnoveg = dmerg[which(complete.cases(dmerg) == F),]
dmerg = dmerg[complete.cases(dmerg),]

write.csv(dnoveg,paste0("counts w missing veg ",sp,".csv"))



# maxperstat = tapply(dat$Count,dat$StationID,max)
# statkeep = as.integer(names(maxperstat)[which(maxperstat > 0)])
# dat = dat[which(dat$StationID %in% statkeep),]




dat = dmerg



dat$yr = dat$YEAR-1997 #year variable starting at 1, if species data include 1998

midyear = floor(max(dat$yr)/2) #middle of time-series, used to center the year effects in the model

yrs = unique(dat$yr) #vector of years with bird data

# dat$obser = as.integer(factor(dat$ObserverID)) #new unique observer ID number with no missing values

# nobservers = max(dat$obser) #number of observers


doy1 = gsub(pattern = "[,].*", replacement = "",gsub(pattern = ".*y[,] ", replacement = "", dat$Date))
doy2 = as.Date(doy1, format = "%B %d")
dat$doy = as.integer(doy2-as.Date("June 1", format = "%B %d"))
dat$doy_s = (dat$doy-15)/10
dat$doy_2_s = dat$doy_s^2


#summary(dat$Time)
tthr = floor(dat$Time/100)
ttmin = dat$Time-(tthr*100)
minpost = tthr*60 + ttmin
dat$time = minpost
dat$time_s = (dat$time-mean(dat$time))/sd(dat$time)
dat$time_2_s = dat$time_s^2



# nvisits = max(dat$Visit) #maximum number of visits, perhaps uncecessary

#dat$sit = as.integer(factor(dat$StandNumber,ordered = T,
#                            levels = as.character(sort(unique(dat$StandNumber))))) 
# generates a unique site ID variable with no missing values, but ordered in the same numeric
# ordering as dat$Site
#dat$sitf = factor(dat$ss.uni)
lyear = tapply(dat$YEAR,dat$ss.uni,max)

siteswveg = as.character(unique(dat[,"ss.uni"]))
#sites that have complete veg data in 3 years
# these are the sites used to estimate the species trajectories including the observed change in habitat
sitesw3veg = siteswveg[which(siteswveg %in% names(lyear)[which(lyear == 2017)])]

#veg.sp
missvegsp = unique(veg.sp[which(complete.cases(veg.sp) == F),"ss.uni"])

if(any(missvegsp %in% sitesw3veg)){
sitesw3veg = sitesw3veg[-which(sitesw3veg %in% missvegsp)]
}
#sites that have complete veg data in 2 years (i.e. missing 2017 veg surveys)
sitesw2veg = siteswveg[-which(siteswveg %in% sitesw3veg)]
# sitesw3veg = sitesw3veg[which(sitesw3veg %in% unique(dat$ss.uni))]



dat$sitf = (factor(dat$ss.uni,ordered = T,levels = c(sitesw3veg,sitesw2veg))) 


dat$sit = as.integer(dat$sitf) 
sitord = unique(dat[,c("sit","ss.uni")])
sitord = sitord[order(sitord$sit),]
nsits = max(dat$sit)

veg.sp = merge(veg.sp,sitord,by = "ss.uni")
veg.sp = veg.sp[order(veg.sp$sit,veg.sp$YEAR),]




dat = dat[order(dat$sit,dat$YEAR),]
#above generates a unique Station ID variable with no missing values, but ordered in the same numeric

apreds = c("doy","time",sppreds)

predcols = paste0(rep(apreds,each = 2),c("_s","_2_s"))

if(addbud){
predcols = c(predcols,"budw")
}

spform = as.formula(paste0("count~YEAR+",paste(predcols,collapse = "+")))

spformg = as.formula(paste0("count~s(YEAR)+",paste(predcols,collapse = "+")))
spformnaive = as.formula("count~YEAR")

######### running preliminary, simplified analysis to remove non-signficant predictors
### this step could be removed

if(pretest){

  dattest = dat
  if(addbud){
  dattest$budw = factor(dattest$budw)
  }
test = glmmadmb(spform,data = dattest,
                family = "poisson",
                random = ~ (1 | sitf))

testnaive = glmmadmb(formula = spformnaive,data = dattest,
                family = "poisson",
                random = ~ (1 | sitf))
st = summary(test)
stnaive = summary(testnaive)

coefout = as.data.frame(st$coefficients)
names(coefout) = paste0("slope",c("",".SE",".Z",".pval"))
coefout[,paste0("slope.naive",c("",".SE",".Z",".pval"))] = NA
coefout["YEAR",paste0("slope.naive",c("",".SE",".Z",".pval"))] = as.numeric(coefficients(stnaive)[2,])

write.csv(coefout,paste0(sp,"all test model coefficients Jan2019.csv"))


pvals = st$coefficients[,4]
sigpreds = names(pvals)[which(pvals < 0.1)]## using a rather liberal "significance" threshold

if(any(grepl(sigpreds,pattern = "budw"))){
  
  sigpreds[grep(sigpreds,pattern = "budw")] = "budw"
  sigpreds = unique(sigpreds)
  
}

if(addbud){  
  if(any(c("(Intercept)","YEAR","budw") %in% sigpreds)){
  sigpreds = sigpreds[-which(sigpreds %in% c("(Intercept)","YEAR","budw"))]
}
}else{
if(any(c("(Intercept)","YEAR") %in% sigpreds)){
  sigpreds = sigpreds[-which(sigpreds %in% c("(Intercept)","YEAR"))]
}
}

sppreds2 = c(apreds[1:2])
for(pp in sigpreds){
  ppk = strsplit(pp,split = "_",fixed = T)[[1]][1]
  sppreds2 = c(sppreds2,ppk)
}
sppreds2 = unique(sppreds2)
predcols2 = paste0(rep(sppreds2,each = 2),c("_s","_2_s"))
 
if(addbud){predcols2 = c(predcols2,"budw")}

spform2 = as.formula(paste0("count~YEAR+",paste(predcols2,collapse = "+")))




}else{

  
  
    sppreds2 = apreds
  predcols2 = predcols
  
  
  spform2 = spform
}





npredictors = length(sppreds2)


hab = dat[,paste0(sppreds2,"_s")]
hab2 = dat[,paste0(sppreds2,"_2_s")]


nyears = 20
nvegpredictors = npredictors-2
nsitespredict = length(sitesw3veg)
habp = array(dim = c(nsitespredict,nyears,nvegpredictors))
hab2p = habp
sppreds2veg = sppreds2[-c(1:2)]

if(addbud){
  habbudp = habp[,,1]

for(s in 1:nsitespredict){
  tmp = veg.sp[which(veg.sp$sit == s),]
  tmp[c((nrow(tmp)+1):(nrow(tmp)+3)),"YEAR"] = c(2000,2006,2007)
  tmp = tmp[order(tmp$YEAR),]
  #if(any(is.na(tmp))){break}
  habp[s,,] = as.matrix(tmp[,paste0(sppreds2veg,"_s")])
  hab2p[s,,] = as.matrix(tmp[,paste0(sppreds2veg,"_2_s")])
  habbudp[s,] = tmp$budw
}#s
predyrs = yrswbud-1997
jagsdat = list(count = dat$count,
               #statn = dat$statn,
               sit = dat$sit,
               yr = dat$yr,
               budw = dat$budw,
               nbudclass = max(dat$budw),
               #visit = dat$Visit,
               #desc = dat$desc,
               #ndescs = ndescs,
               #nobservers = nobservers,
               #yrs = yrs,
               predyrs = predyrs,
               nyrs = max(dat$yr),
               midyear = midyear,
               #nstatns = nstatns,
               nsits = nsits,
               #nvisits = nvisits,
               ncounts = nrow(dat),
               qzeros = c(0,0),
               npredictors = npredictors,
               nvegpredictors = nvegpredictors,
               nsitespredict = nsitespredict,
               diagon = matrix(c(1,0,0,1),
                               nrow = 2,byrow = T),
               hab = hab,#observed habitat values
               hab2 = hab2,#observed habitat values (squared)
               habp = habp,#reduced matrix of habitat values for sites with habitat in all years
               hab2p = hab2p,
               habbudp = habbudp)



sp.params = c("beta",
              "alpha",
              "sdnoise",
              "n",
              "nh",
              "trend",
              "trendh",
              "trenddiff",
              "sdsite",
              #"sdstation",
              "site",
              #"station",
              "bhab",
              "bbud")

mod <- (paste0(dir.model,"/Liard model w bud.txt"))  

}else{
  
  for(s in 1:nsitespredict){
    tmp = veg.sp[which(veg.sp$sit == s),]
    tmp = tmp[order(tmp$YEAR),]
    #if(any(is.na(tmp))){break}
    habp[s,,] = as.matrix(tmp[,paste0(sppreds2veg,"_s")])
    hab2p[s,,] = as.matrix(tmp[,paste0(sppreds2veg,"_2_s")])
  }#s

  predyrs = 1:nyears

jagsdat = list(count = dat$count,
               #statn = dat$statn,
               sit = dat$sit,
               yr = dat$yr,
               #visit = dat$Visit,
               #desc = dat$desc,
               #ndescs = ndescs,
               #nobservers = nobservers,
               #yrs = yrs,
               predyrs = predyrs,
               nyrs = max(dat$yr),
               midyear = midyear,
               #nstatns = nstatns,
               nsits = nsits,
               #nvisits = nvisits,
               ncounts = nrow(dat),
               qzeros = c(0,0),
               npredictors = npredictors,
               nvegpredictors = nvegpredictors,
               nsitespredict = nsitespredict,
               diagon = matrix(c(1,0,0,1),
                               nrow = 2,byrow = T),
               hab = hab,#observed habitat values
               hab2 = hab2,#observed habitat values (squared)
               habp = habp,#reduced matrix of habitat values for sites with habitat in all years
               hab2p = hab2p)



sp.params = c("beta",
              "alpha",
              "sdnoise",
              "n",
              "nh",
              "trend",
              "trendh",
              "trenddiff",
              "sdsite",
              #"sdstation",
              "site",
              #"station",
              "bhab")

mod <- (paste0(dir.model,"/Liard model no bud.txt"))  

  
  
  
  
}



adaptSteps = 1000              # Number of steps to "tune" the samplers.
burnInSteps = 70000            # Number of steps to "burn-in" the samplers.
nChains = 3                   # Number of chains to run.
numSavedSteps=10000           # Total number of steps to save.
thinSteps=300                   # Number of steps to "thin" (1=keep every step).
nIter = burnInSteps + ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.

t1 = Sys.time()


jagsMod = jags( model.file = mod, 
                      data= jagsdat ,  
                      #inits= sp.inits,  
                      n.chains= nChains , 
                      #n.adapt= NULL,
                n.iter=nIter,
                n.burnin = burnInSteps,
                n.thin = thinSteps,
                parallel = T,
                parameters.to.save = sp.params,
                store.data = T)


t2 = Sys.time()
t2-t1


# The saved MCMC chain:
pdf(paste(sp,"post jags diagnostic Jan2019.pdf"))
plot(jagsMod$samples)
dev.off()


save(list = c("jagsMod","sp","t2","t1","dat","jagsdat"),
     file = paste(sp, "modeling output Jan2019.RData"))

sums = as.data.frame(jagsMod$summary)


### calculating and visualizing the trajectories based on trends



### calculating and visualizing the trajectories based on trends


### calculating and visualizing the trajectories based on trends


### calculating and visualizing the trajectories based on trends

#str(jagsdat)
names(sums)[3:7] <- c("lci","lquart","med","uquart","uci")
indices = as.data.frame(sums[paste0("n[",1:nyears,"]"),])
indices$year = (1:nyears)+1997
indices$version = "Accounting for habitat"

indicesnohab = as.data.frame(sums[paste0("nh[",1:nyears,"]"),])
indicesnohab$year = (1:nyears)+1997
indicesnohab$version = "Including habitat"

indices = rbind(indices,indicesnohab)
#######
#indices is a dataframe of two sets of annual indices
#Including habitat = "nh[y]" = predicted annual indices, including the effects of habitat change
#Accounting for habitat = "n[y]" = predicted annual indices, after removning the effects of any changes in habitat




beta = sums["beta",]
beta[1:7] <- (exp(beta[1:7])-1)*100
slopplot = paste0(sp," trend = ",signif(beta$med,2),"%/year [",
                 signif(beta$lci,2), " : ",signif(beta$uci,2),"]")

trend = sums["trend",]
trendh = sums["trendh",]
trenddiff = sums["trenddiff",]

tout = paste0(signif(trend[,"med"],2)," [",signif(trend[,"lci"],2),":",signif(trend[,"uci"],2),"]")
thout = paste0(signif(trendh[,"med"],2)," [",signif(trendh[,"lci"],2),":",signif(trendh[,"uci"],2),"]")
tdiffout = paste0(sp,"difference in trends = ", signif(trenddiff[,"med"],2)," [",signif(trenddiff[,"lci"],2),":",signif(trenddiff[,"uci"],2),"]")

trendsout = rbind(trend,
                  trendh,
                  trenddiff,
                  beta)
trendsout$trend = c("endpoint habitat removed",
                    "endpoint with habitat",
                    "difference in endpoints",
                    "slope parameter")
write.csv(trendsout,paste0(sp,"trend estimates.csv"))

kp = which(dat$ss.uni %in% sitesw3veg)
obsmeans = tapply(dat[kp,"count"],dat[kp,"YEAR"],mean)
obsmeans = data.frame(YEAR = as.numeric(as.character(names(obsmeans))),
                      obs = as.numeric(obsmeans))
### plotting of the trajectories
indend = indices[which(indices$year == 2017),]
typ = indend[which(indend$version == "Accounting for habitat"),"med"]
typh = indend[which(indend$version == "Including habitat"),"med"]

indiceso = indices
indices = indices[which(!is.na(indices$mean)),]

pdf(paste0(sp,"indices Jan2019.pdf"))
  
  p = ggplot(data = indices, aes(x = year,y = med))+
    geom_ribbon(data = indices,aes(x = year,ymin = lci,ymax = uci,fill = version),
                alpha = 0.25)+
    geom_line(data = indices, aes(x = year,y = med,color = version))+
    geom_point(data = obsmeans,aes(x = YEAR,y = obs),col = "black",alpha = 0.5)+
    labs(x = "",y = "index of abundance",title = tdiffout)+
    annotate(geom = "text",x = 2017,y = typ,label = tout)+
    annotate(geom = "text",x = 2017,y = typh,label = thout)+
    coord_cartesian(xlim = c(1998,2019),ylim = c(0,max(max(indices$uci),max(obsmeans[,2],na.rm = T))))
  print(p)
  

dev.off()




#### plotting of the veg and survey-timing effects



#### plotting of the veg and survey-timing effects



#### plotting of the veg and survey-timing effects



#### plotting of the veg and survey-timing effects



bhabsn = paste0("bhab[",rep(1:npredictors,each = 2),",",rep(c(1,2),times = npredictors),"]")

bhabs = as.data.frame(sums[bhabsn,])
bhabs$predictor = c(paste0(rep(sppreds2,each = 2),rep(c("_s","_2_s"),times = npredictors)))

### bhabs = dataframe of the coefficient estimates for the two
###         timing predictors (doy and time), and the other habitat predictors

if(addbud){
bbudsn = paste0("bbud[",1:5,"]")
bbud = as.data.frame(sums[bbudsn,])
 
bbud$predictor = paste0("budworm defol ",c("=","<","<","<",">"),cats[c(2,3,4,5,5)]/1000,"km")

bhabsout = rbind(bhabs,bbud)
}else{
  bhabsout = bhabs
}
write.csv(bhabsout,paste0(sp," coefficient estimates Jan2019.csv"))



#### create predicted values to visulize the effects of the vege and survey timing predictors
posterior = MCMCchains(jagsMod$samples,params = "bhab")
vegeffects = veg.predict
vegeffects$doy = seq(min(dat$doy),max(dat$doy),length.out = Npreds) 
vegeffects$doy_s = (vegeffects$doy-15)/10
vegeffects$doy_2_s = vegeffects$doy_s^2
vegeffects$time = seq(min(dat$time),max(dat$time),length.out = Npreds)
vegeffects$time_s = (vegeffects$time-mean(dat$time))/sd(dat$time)
vegeffects$time_2_s = vegeffects$time_s^2

interc = MCMCchains(jagsMod$samples,params = "alpha")

#############

for(v in seq(1,(nrow(bhabs)-1),by = 2)){
b = bhabsn[v]
b2 = bhabsn[v+1]
bn = bhabs[v,"predictor"]
bn2 = bhabs[v+1,"predictor"]

postvec = posterior[,b]
postvec2 = posterior[,b2]

for(i in 1:nrow(vegeffects)){
tmp = exp(interc + postvec*vegeffects[i,bn]+postvec2*vegeffects[i,bn2])

vegeffects[i,paste0("pred.",bn)] = median(tmp)
vegeffects[i,paste0("pred.",bn,".lci")] = quantile(tmp,0.025)
vegeffects[i,paste0("pred.",bn,".uci")] = quantile(tmp,0.975)

}



}#v
 
if(addbud){
bbudpost = MCMCchains(jagsMod$samples,params = "bbud")
for(j in 1:nrow(bbud)){
  
  tmp = exp(interc + bbudpost[,j])
  
  
  bbud[j,"n.birds"] = median(tmp)
  bbud[j,"n.birds.lci"] = quantile(tmp,0.025)
  bbud[j,"n.birds.uci"] = quantile(tmp,0.975)
  
  
}
bbud$pred = gsub(bbud$predictor,pattern = " defol",replacement = "")
bbud$pred = factor(bbud$pred,ordered = T,levels = bbud$pred)


pdf(paste0(sp,"predictor effects with budworm.pdf"))
}else{
  pdf(paste0(sp,"predictor effects no budworm.pdf"))
  
}
for(v in seq(1,(nrow(bhabs)-1),by = 2)){
  bn = bhabs[v,"predictor"]
  bnn = gsub(bn,pattern = "_s",replacement = "")
  
  tmp = vegeffects[,c(bnn,bn,paste0("pred.",bn),
                      paste0("pred.",bn,".lci"),
                      paste0("pred.",bn,".uci"))]
  names(tmp) = c("pred","preds","med","lci","uci")
  datmp = dat
  names(datmp)[which(names(datmp) == bnn)] <- "pred"
  brk = quantile(datmp$pred,probs = seq(0,1,length = 20),names = F)
  brk = unique(brk)
  mbrk = brk[-length(brk)]+(diff(brk)/2)
  
  datmp$predc = cut(datmp$pred,breaks = brk)
  
  obsm = tapply(datmp$count,datmp$predc,mean,na.rm = T)
  obslci = tapply(datmp$count,datmp$predc,quantile,probs = 0.25,na.rm = T)
  obsuci = tapply(datmp$count,datmp$predc,quantile,probs = 0.75,na.rm = T)
  datgr = data.frame(pred = mbrk,
                     mn = obsm,
                     lci = obslci,
                     uci = obsuci)
  datlg = datgr[which.min(datgr$mn),]
  # datlg$mn = max(tmp$uci)
  # datlg$pred = as.numeric(quantile(tmp$pred,0.9))
  
  p = ggplot(data = tmp, aes(x = pred,y = med))+
    geom_pointrange(data = datgr,aes(x = pred,y = mn,ymin = lci,ymax = uci),
                    colour = grey(0.8))+
    #geom_point(data = datgr,aes(x = pred,y = mn),shape = 20,size = 2)+
    geom_ribbon(data = tmp,aes(x = pred,ymin = lci,ymax = uci),
                fill = grey(0.5),alpha = 0.5)+
    geom_line(data = tmp, aes(x = pred,y = med))+
    #geom_point(data = datmp,aes(x = pred,y = count),shape = 1,alpha = 0.3)+
    labs(x = bnn,y = "number of birds",title = paste(sp,bnn))+
    geom_point(data = datlg,aes(x = pred,y = mn),shape = 20,size = 2,colour = grey(0.8))+
    geom_text_repel(data = datlg,aes(x = datlg$pred,y = datlg$mn),label = "Observed means",colour = grey(0.4))
  print(p)
  
}#v

if(addbud){
p = ggplot(data = bbud, aes(x = pred,y = n.birds))+
  geom_pointrange(data = bbud,aes(x = pred,ymin = n.birds.lci,ymax = n.birds.uci),
              colour = grey(0.5))+
   geom_point(data = bbud,aes(x = pred,y = n.birds),shape = 20,size = 5)+
  labs(x = "budworm proximity",y = "predicted mean number of birds",title = "Budworm effects")
print(p)
}


dev.off()









print(paste(sp,"completed",t2-t1))


}### end of species loop






















# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# out = summary(cd13)
# 
# indices = out$quantiles
# indices = indices[paste0("n[",1:max(yrs),"]"),]
# 
# 
# qs = out$quantiles
# 
# ### species trend
# slopeqs = qs["beta",]
# 
# trendqs = round((exp(slopeqs)-1)*100,1)
# 
# print(paste0(sp," trend estimate = ",trendqs["50%"],"%/year [",trendqs["2.5%"]," - ",trendqs["97.5%"],"]"))
# 
# x11()
# source("transparency function.r")
# 
# 
# obsmean = tapply(dat$Count,dat$YEAR,mean)
# yup = max(c(quantile(dat$Count,0.95),max(indices[,"97.5%"])))
# xplot = c(min(dat$YEAR):max(dat$YEAR))
# plot(y = indices[,"50%"],
#      x = xplot,
#      ylim = c(0,yup),
#      type = "l")#estimated trajectory
# polygon(y = c(indices[,"2.5%"],rev(indices[,"97.5%"])),
#         x = c(xplot,rev(xplot)),
#         border = NA,
#         col = transp.func(grey(0.5),0.2))#95% credible interval
# points(y = obsmean,
#        x = as.integer(names(obsmean)),
#        pch = 19,
#        cex = 1.5,
#        col = "red")#observed mean counts by year
# points(y = dat$Count+rnorm(1:nrow(dat),0,0.05),
#        x = dat$YEAR+rnorm(1:nrow(dat),0,0.05),
#        cex = 0.5,
#        col = transp.func("darkgreen",0.3))#jittered counts - raw data
# 
# 








