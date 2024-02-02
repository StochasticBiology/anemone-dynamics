library(ggplot2)
library(ggpubr)
library(writexl)

############ code layout:
# helper functions
# data input
# the set of tables which are standlone analyses (mainly LMs) of specific data
# bootstrap analysis for multiphase behaviour
# the set of tables which pull together information from more than one of the above

############ helper functions

# pull in likelihood functions for multi-phase models
source("multiphasemodel.R")

# extract various statistics from lm summaries
pull.stats = function(lm.fit, nvar = 2, half.life = FALSE) {
  s = summary(lm.fit)
  if(half.life == TRUE) {
    if(s$coefficients[2,1] > 0) {
      return(data.frame(ci.lo = log(2)/(s$coefficients[2,1]+1.96*s$coefficients[2,2]),
                        ci.hi = log(2)/(s$coefficients[2,1]-1.96*s$coefficients[2,2])))
    } else {
      return( data.frame(ci.lo = log(2)/abs(s$coefficients[2,1]-1.96*s$coefficients[2,2]),
                         ci.hi = log(2)/abs(s$coefficients[2,1]+1.96*s$coefficients[2,2])))
    }
  }
  if(nvar == 2) {
    return(data.frame( slope=s$coefficients[2,1], slope.sd=s$coefficients[2,2], 
                       intercept=s$coefficients[1,1], intercept.sd=s$coefficients[1,2], 
                       R2=s$r.squared, AIC=AIC(lm.fit) ))
  } else {
    return(data.frame( slope.1=s$coefficients[2,1], slope.1.sd=s$coefficients[2,2], 
                       slope.2=s$coefficients[3,1], slope.2.sd=s$coefficients[3,2], 
                       intercept=s$coefficients[1,1], intercept.sd=s$coefficients[1,2], 
                       R2=s$r.squared, AIC=AIC(lm.fit) )) 
  }
}

# return mean and sem-based ci size
mean.sem = function(x) {
  return(data.frame(mean=mean(x), sem=1.96*sd(x)/sqrt(length(x))))
}

# return mean and sd of estimator from "pulled stats"
mean.sd.lm = function(x) {
  return(data.frame(mean=x[1,1], sd=x[1,2]))
}

# return mean and sd
mean.sd = function(x) {
  return(data.frame(mean=mean(x), sd=sd(x)))
}

# number of bootstrap resamples to run
total.nboot = 100

theme_set(theme_light())

############ read in data

# read new data on body size, cell size, cell number for growth, shrinkage_new, and shrinkage_old experiments
df.new = read.csv2("Data/NEWdata_shrinkage_new_shrinkage_old_growth.csv")
df.new = df.new[df.new$day != "20RF" & df.new$day != 22,]
df.new$day = as.numeric(as.character(df.new$day))

# read data on body size across different time course experiments
df.all = read.csv2("Data/ALLDATA_bodysize.csv")
df.all = df.all[df.all$day != "20RF",]
df.all$day = as.numeric(as.character(df.all$day))

# read data on Aiptasia
df.aip = read.csv2("Data/longer-starved_aiptasia.txt", sep="\t")
df.aip = df.aip[df.aip$day != "20RF",]
df.aip$day = as.numeric(as.character(df.aip$day))

# summary plot new data
ggplot(df.new, aes(x=day,y=log(mm2))) + geom_point() + facet_wrap(~experiment, scales="free_x")
# summary plot all other data
ggplot(df.all, aes(x=day,y=log(mm2))) + geom_point() + facet_wrap(~experiment, scales="free_x")
# summary plot Aiptasia
ggplot(df.aip, aes(x=day,y=log(mm2))) + geom_point() + facet_wrap(~symbiosis, scales="free_x")

############ new post-review section: growth rate with body size

# assign new labels to individuals that pulls the final element from the current label
# e.g. S_1_01 -> 01; S_21d_24 -> 24
df.new$i.label = ""
for(i in 1:nrow(df.new)) {
  ss = strsplit(df.new$individual[i], split="_")[[1]]
  df.new$i.label[i] = ss[length(ss)]
}

df.new$lin.rate = NA
df.new$exp.rate = NA
expts = unique(df.new$experiment)
# loop over different experiments (growth, shrinkage old/new)
for(expt in expts) {
  for(i in 1:nrow(df.new)) {
    if(df.new$experiment[i] == expt) {
      # get the set of rows in this experiment that correspond to the time series for this individual
      this.one = df.new$i.label[i]
      this.one.refs = which(df.new$experiment == expt & df.new$i.label == this.one)
      # get the set of future times for this individual
      future.refs = this.one.refs[which(this.one.refs > i)]
      if(length(future.refs) > 0) {
        # estimate growth rate by comparing with the closest future datapoint
        future.ref = future.refs[1]
        timediff = df.new$day[future.ref]-df.new$day[i]
        
        # linear difference -- simple s2-s1
        lin.sizediff = df.new$mm2[future.ref]-df.new$mm2[i]
        # exponential difference:
        # s2 = s0 exp(a t2), s1 = s0 exp(a t1)
        # s2/s1 = exp(a(t2-t1))
        # log s2/s1 = a (t2-t1)
        exp.sizediff = log(df.new$mm2[future.ref]/df.new$mm2[i])
        
        df.new$lin.rate[i] = lin.sizediff/timediff
        df.new$exp.rate[i] = exp.sizediff/timediff
      }
    }
  }
}

# compare absolute (less relevant) and exponential (more relevant) behaviours
ggarrange(ggplot(df.new, aes(x=mm2, y=lin.rate)) + geom_point() + geom_smooth(method="lm") + facet_wrap(~ experiment),
          ggplot(df.new, aes(x=mm2, y=exp.rate)) + geom_point() + geom_smooth(method="lm") + facet_wrap(~ experiment),
          nrow=2)

# look at relationships after different time periods post-treatment
ggarrange( 
  ggplot(df.new[df.new$experiment != "shrinkage_old",], aes(x=mm2, y=exp.rate, color=factor(day), fill=factor(day))) + geom_point() + geom_smooth(method="lm", alpha=0.25) + facet_wrap(~ experiment) + xlim(0,8) + labs(x="Body size / mm2", y="Rate of exponential dynamics"),
  ggplot(df.new[df.new$experiment != "shrinkage_old",], aes(x=mm2, y=exp.rate)) + geom_point() + geom_smooth(method="lm", alpha=0.25) + facet_wrap(~ experiment) + xlim(0,8) + labs(x="Body size / mm2", y="Rate of exponential dynamics"),
  nrow=2
)

############ standalone tables

###### Table 1 -- linear model stats for growth behaviours
t.1.df = data.frame()
### Table 1, first three rows
sub = df.new[df.new$experiment == "growth",]
ggplot(sub, aes(x=day, y=log(mm2))) + geom_point() + geom_smooth(method="lm")
t.1.df = rbind(t.1.df, data.frame(label="body size", 
                                  pull.stats(lm(log(mm2) ~ day, data=sub))))
ggplot(sub, aes(x=day, y=log(cellnumber))) + geom_point() + geom_smooth(method="lm")
t.1.df = rbind(t.1.df, data.frame(label="cell number", 
                                  pull.stats(lm(log(cellnumber) ~ day, data=sub))))
ggplot(sub, aes(x=day, y=log(cellsizemedian))) + geom_point() + geom_smooth(method="lm")
t.1.df = rbind(t.1.df, data.frame(label="cell size", 
                                  pull.stats(lm(log(cellsizemedian) ~ day, data=sub))))

### Table 1, last three rows (refeeding, partitioned by period or not partitioned)
sub = df.all[df.all$experiment=="170d_200d_starved_refed",]
ggplot(sub, aes(x=day, y=log(mm2), color=condition)) + geom_point() + geom_smooth(method="lm")
t.1.df = rbind(t.1.df, data.frame(label="body size, 170d", 
                                  pull.stats(lm(log(mm2) ~ day, data=sub[sub$condition=="170d",]))))
t.1.df = rbind(t.1.df, data.frame(label="body size, 200d", 
                                  pull.stats(lm(log(mm2) ~ day, data=sub[sub$condition=="200d",]))))
t.1.df = rbind(t.1.df, data.frame(label="body size, 170+200d", 
                                  pull.stats(lm(log(mm2) ~ day, data=sub))))

####### Table 2+3 -- relationships between body size, cell number, cell size
t.2.df = data.frame()
sub = df.new[df.new$experiment == "growth",]
sub$day = factor(sub$day)
ggplot(sub, aes(x=log(cellnumber), y=log(mm2), color=day)) + geom_point()
ggplot(sub, aes(x=log(cellsizemedian), y=log(mm2), color=day)) + geom_point()

t.2.df = rbind(t.2.df, data.frame(label="body size ~ cell size + cell number", 
                                  pull.stats(lm(log(mm2) ~ log(cellsizemedian) + log(cellnumber), data=sub), nvar=3)))
t.3.df = data.frame()
t.3.df = rbind(t.3.df, data.frame(label="body size ~ cell size", 
                                  pull.stats(lm(log(mm2) ~ log(cellsizemedian), data=sub))))
t.3.df = rbind(t.3.df, data.frame(label="body size ~ cell number", 
                                  pull.stats(lm(log(mm2) ~ log(cellnumber), data=sub))))

####### Tables 6 and 7 -- relationships between body size, cell number, cell size partitioned by time period
sub = df.new[df.new$experiment == "shrinkage_new",]
sub$Day = factor(sub$day)
ggplot(sub, aes(x=log(cellnumber), y=log(mm2), color=Day)) + geom_point()
ggplot(sub, aes(x=log(cellsizemedian), y=log(mm2), color=Day)) + geom_point()

### Table 6
t.6.df = data.frame()
t.6.df = rbind(t.6.df, data.frame(label="body size ~ cell size + cell number, t <= 2",
                                  pull.stats(lm(log(mm2) ~ log(cellsizemedian) + log(cellnumber), data=sub[sub$day <= 2,]), nvar=3) ))
t.6.df = rbind(t.6.df, data.frame(label="body size ~ cell size + cell number, t >= 2",
                                  pull.stats(lm(log(mm2) ~ log(cellsizemedian) + log(cellnumber), data=sub[sub$day >= 2,]), nvar=3)))

### Table 7
t.7.df = data.frame()
t.7.df = rbind(t.7.df, data.frame(label="body size ~ cell size, t <= 2",
                                  pull.stats(lm(log(mm2) ~ log(cellsizemedian), data=sub[sub$day <= 2,])) ))
t.7.df = rbind(t.7.df, data.frame(label="body size ~ cell number, t <= 2",
                                  pull.stats(lm(log(mm2) ~ log(cellnumber), data=sub[sub$day <= 2,])) ))
t.7.df = rbind(t.7.df, data.frame(label="body size ~ cell size, t >= 2",
                                  pull.stats(lm(log(mm2) ~ log(cellsizemedian), data=sub[sub$day >= 2,])) ))
t.7.df = rbind(t.7.df, data.frame(label="body size ~ cell number, t >= 2",
                                  pull.stats(lm(log(mm2) ~ log(cellnumber), data=sub[sub$day >= 2,])) ))

####### Table 10 -- linear models for time behaviour of body size, cell number, cell size, partitioned by time interval
t.10.df = data.frame()
sub = df.new[df.new$experiment == "shrinkage_new",]
t.10.df = rbind(t.10.df, data.frame(label="body size, 2 <= t <= 5",
                                    pull.stats(lm(log(mm2) ~ day, data=sub[sub$day >= 2 & sub$day <= 5,])) ))
t.10.df = rbind(t.10.df, data.frame(label="cell number, 2 <= t <= 5",
                                    pull.stats(lm(log(cellnumber) ~ day, data=sub[sub$day >= 2 & sub$day <= 5,])) ))
t.10.df = rbind(t.10.df, data.frame(label="cell size, 2 <= t <= 5",
                                    pull.stats(lm(log(cellsizemedian) ~ day, data=sub[sub$day >= 2 & sub$day <= 5,])) ))
t.10.df = rbind(t.10.df, data.frame(label="body size, 5 <= t <= 12",
                                    pull.stats(lm(log(mm2) ~ day, data=sub[sub$day >= 5 & sub$day <= 12,])) ))
t.10.df = rbind(t.10.df, data.frame(label="cell number, 5 <= t <= 12",
                                    pull.stats(lm(log(cellnumber) ~ day, data=sub[sub$day >= 5 & sub$day <= 12,])) ))
t.10.df = rbind(t.10.df, data.frame(label="cell size, 5 <= t <= 12",
                                    pull.stats(lm(log(cellsizemedian) ~ day, data=sub[sub$day >= 5 & sub$day <= 12,])) ))
t.10.df = rbind(t.10.df, data.frame(label="body size, 12 <= t <= 21",
                                    pull.stats(lm(log(mm2) ~ day, data=sub[sub$day >= 12 & sub$day <= 21,])) ))
t.10.df = rbind(t.10.df, data.frame(label="cell number, 12 <= t <= 21",
                                    pull.stats(lm(log(cellnumber) ~ day, data=sub[sub$day >= 12 & sub$day <= 21,])) ))
t.10.df = rbind(t.10.df, data.frame(label="cell size, 12 <= t <= 21",
                                    pull.stats(lm(log(cellsizemedian) ~ day, data=sub[sub$day >= 12 & sub$day <= 21,])) ))

###### fitting for the various multiphase models. this section will take a few minutes

##
# test different optimisation protocols on an awkward dataset
# consider changing the "expt" parameter in initial.guess and "method" in optim if we're not converging well
data.set = "RES"
sub = df.all[df.all$experiment == "Cyclicfeeding_starving" & !is.na(log(df.all$mm2)),]
bx = sub$day[sub$condition == data.set]
by = log(sub$mm2[sub$condition == data.set]) 
model.type = 4
theta = initial.guess(expt=0, model.type, bx, by)
opt = optim(theta, neg.log.likelihood, model=model.type, y=by, x=bx, method="BFGS")
2*opt$val + 2*(model.type*2+1)
##

# initial structures for bootstrap output
boot.df = data.frame()
plot.df = data.frame()

# loop over different experiments
for(data.set in c("AL", "RES", "size21", "cellnum21", "cellsize21", "140d_starved", "200d_starved", "apo", "sym")) {
  # pull predictor-response data for this experiment
  if(data.set == "AL" | data.set == "RES") {
    sub = df.all[df.all$experiment == "Cyclicfeeding_starving" & !is.na(log(df.all$mm2)),]
    x = sub$day[sub$condition == data.set]
    y = log(sub$mm2[sub$condition == data.set]) 
  } else if(data.set == "140d_starved" | data.set == "200d_starved") {
    sub = df.all[df.all$experiment == data.set & !is.na(log(df.all$mm2)),]
    x = sub$day
    y = log(sub$mm2)
  } else if(data.set == "apo" | data.set == "sym") {
    sub = df.apt[df.apt$symbiosis == data.set,]
    x = sub$day
    y = log(sub$mm2)
  } else {
    sub = df.new[df.new$experiment == "shrinkage_new",]
    x = sub$day
    if(data.set == "size21") { 
      y = log(sub$mm2)
    } else if(data.set == "cellnum21") {
      y = log(sub$cellnumber)
    } else if(data.set == "cellsize21") {
      y = log(sub$cellsizemedian)
    }
  }
  
  # loop through changepoint numbers
  for(model.type in 1:4) { 
    # loop over bootstrap resamples
    for(nboot in 1:total.nboot) {
      print(paste(c(data.set, " model ", model.type, " bootstrap ", nboot, " of ", total.nboot), collapse=""))
      # for the first sample, retain original data, otherwise sample with replacement
      if(nboot == 1) {  boot = 1:length(y) }
      else {  boot = sample(length(y), length(y), replace=T) }
      bx = x[boot]
      by = y[boot]
      # initial guess for parameters -- reflecting a flat profile with changepoints based on prelim knowledge      
      theta = initial.guess(0, model.type, bx, by)
      
      # attempt to optimise
      opt = optim(theta, neg.log.likelihood, model=model.type, y=by, x=bx, method="BFGS")
      
      # append slopes to lists, this only for saving and plotting the model outputs of distributions of slopes
      boot.df = rbind(boot.df, data.frame(data.set = data.set, model=model.type,boot=nboot,
                                          b1=opt$par[1],b2=opt$par[4],b3=opt$par[6],b4=opt$par[8],
                                          tau1=opt$par[5], tau2=opt$par[7], tau3=opt$par[9],
                                          AIC = 2*opt$val + 2*(model.type*2+1) ))
      # add this model prediction to plotting data frame
      xmodel = floor(min(x)):ceiling(max(x))
      ymodel = unlist(lapply(xmodel, model.val, theta=opt$par, model=model.type)) 
      plot.df = rbind(plot.df, data.frame(data.set = data.set, model=model.type, boot=nboot, x=xmodel, y=ymodel))
    }
  }
}

## this awkward code just prepares a set of plots of each model's bootstrap resamples along with its original data, labelled by AIC value
expt.ref = "AL"
mod.labs = c("1" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 1 & boot.df$boot == 1]),
             "2" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 2 & boot.df$boot == 1]),
             "3" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 3 & boot.df$boot == 1]),
             "4" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 4 & boot.df$boot == 1]))
g.mp.1 = ggplot() + geom_point(data= df.all[df.all$experiment == "Cyclicfeeding_starving" & df.all$condition == "AL" & !is.na(log(df.all$mm2)),],
                      aes(x=day,y=log(mm2))) + geom_line(data = plot.df[plot.df$data.set==expt.ref,], aes(x=x,y=y,color=factor(boot))) + 
  theme(legend.position = "none") + facet_wrap(~ model, labeller = labeller(model = mod.labs)) 

expt.ref = "RES"
mod.labs = c("1" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 1 & boot.df$boot == 1]),
             "2" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 2 & boot.df$boot == 1]),
             "3" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 3 & boot.df$boot == 1]),
             "4" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 4 & boot.df$boot == 1]))
g.mp.2 = ggplot() + geom_point(data= df.all[df.all$experiment == "Cyclicfeeding_starving" & df.all$condition == "RES" & !is.na(log(df.all$mm2)),],
                      aes(x=day,y=log(mm2))) + geom_line(data = plot.df[plot.df$data.set==expt.ref,], aes(x=x,y=y,color=factor(boot))) + 
  theme(legend.position = "none") + facet_wrap(~ model, labeller = labeller(model = mod.labs)) 

expt.ref = "140d_starved"
mod.labs = c("1" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 1 & boot.df$boot == 1]),
             "2" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 2 & boot.df$boot == 1]),
             "3" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 3 & boot.df$boot == 1]),
             "4" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 4 & boot.df$boot == 1]))
g.mp.3 = ggplot() + geom_point(data= df.all[df.all$experiment == "140d_starved" & !is.na(log(df.all$mm2)),],
                      aes(x=day,y=log(mm2))) + geom_line(data = plot.df[plot.df$data.set==expt.ref,], aes(x=x,y=y,color=factor(boot))) + 
  theme(legend.position = "none") + facet_wrap(~ model, labeller = labeller(model = mod.labs))

expt.ref = "200d_starved"
mod.labs = c("1" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 1 & boot.df$boot == 1]),
             "2" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 2 & boot.df$boot == 1]),
             "3" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 3 & boot.df$boot == 1]),
             "4" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 4 & boot.df$boot == 1]))
g.mp.4 = ggplot() + geom_point(data= df.all[df.all$experiment == "200d_starved" & !is.na(log(df.all$mm2)),],
                      aes(x=day,y=log(mm2))) + geom_line(data = plot.df[plot.df$data.set==expt.ref,], aes(x=x,y=y,color=factor(boot))) + 
  theme(legend.position = "none") + facet_wrap(~ model, labeller = labeller(model = mod.labs))

expt.ref = "size21"
mod.labs = c("1" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 1 & boot.df$boot == 1]),
             "2" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 2 & boot.df$boot == 1]),
             "3" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 3 & boot.df$boot == 1]),
             "4" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 4 & boot.df$boot == 1]))
g.mp.5 = ggplot() + geom_point(data= df.new[df.new$experiment == "shrinkage_new",],
                      aes(x=day, y=log(mm2))) + geom_line(data = plot.df[plot.df$data.set==expt.ref,], aes(x=x,y=y,color=factor(boot))) + 
  theme(legend.position = "none") + facet_wrap(~ model, labeller = labeller(model = mod.labs))

expt.ref = "cellnum21"
mod.labs = c("1" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 1 & boot.df$boot == 1]),
             "2" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 2 & boot.df$boot == 1]),
             "3" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 3 & boot.df$boot == 1]),
             "4" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 4 & boot.df$boot == 1]))
g.mp.6 = ggplot() + geom_point(data= df.new[df.new$experiment == "shrinkage_new",],
                      aes(x=day, y=log(cellnumber))) + geom_line(data = plot.df[plot.df$data.set==expt.ref,], aes(x=x,y=y,color=factor(boot))) + 
  theme(legend.position = "none") + facet_wrap(~ model, labeller = labeller(model = mod.labs))

expt.ref = "cellsize21"
mod.labs = c("1" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 1 & boot.df$boot == 1]),
             "2" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 2 & boot.df$boot == 1]),
             "3" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 3 & boot.df$boot == 1]),
             "4" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 4 & boot.df$boot == 1]))
g.mp.7 = ggplot() + geom_point(data= df.new[df.new$experiment == "shrinkage_new",],
                      aes(x=day, y=log(cellsizemedian))) + geom_line(data = plot.df[plot.df$data.set==expt.ref,], aes(x=x,y=y,color=factor(boot))) + 
  theme(legend.position = "none") + facet_wrap(~ model, labeller = labeller(model = mod.labs))

expt.ref = "apo"
mod.labs = c("1" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 1 & boot.df$boot == 1]),
             "2" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 2 & boot.df$boot == 1]),
             "3" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 3 & boot.df$boot == 1]),
             "4" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 4 & boot.df$boot == 1]))
g.mp.8 = ggplot() + geom_point(data= df.aip[df.aip$symbiosis == "apo",],
                               aes(x=day, y=log(mm2))) + geom_line(data = plot.df[plot.df$data.set==expt.ref,], aes(x=x,y=y,color=factor(boot))) + 
  theme(legend.position = "none") + facet_wrap(~ model, labeller = labeller(model = mod.labs))

expt.ref = "sym"
mod.labs = c("1" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 1 & boot.df$boot == 1]),
             "2" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 2 & boot.df$boot == 1]),
             "3" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 3 & boot.df$boot == 1]),
             "4" = as.character(boot.df$AIC[boot.df$data.set == expt.ref & boot.df$model == 4 & boot.df$boot == 1]))
g.mp.9 = ggplot() + geom_point(data= df.aip[df.aip$symbiosis == "sym",],
                               aes(x=day, y=log(mm2))) + geom_line(data = plot.df[plot.df$data.set==expt.ref,], aes(x=x,y=y,color=factor(boot))) + 
  theme(legend.position = "none") + facet_wrap(~ model, labeller = labeller(model = mod.labs))

# pull all these plots together into a multi-panel figure
sf = 2
png("multiphase-models.png", width=800*sf, height=800*sf, res=72*sf)
ggarrange(g.mp.1, g.mp.2, g.mp.3, g.mp.4, g.mp.5, g.mp.6, g.mp.7, g.mp.8, g.mp.9,
          labels = c("AL", "RES", "140d", "200d", "bs", "cn", "cs", "apo", "sym"),
          #labels = c("AL", "RES", "140d_starved", "200d_starved", "size21", "cellnum21", "cellsize21"),
          nrow=3, ncol=3)
dev.off()

############ outputs that pull various subsets of the above together

###### Table 4 -- estimates of rates from different experiments
# gather slope estimates from bootstrapping. the model that is selected for each experiment is on us, based on AIC etc
# for some of these, the changepoints may be inferred to lie outside the range of the timescale. we'll lazily subset these out -- a better solution would be to constraint the changepoint intervals
t.4.df = rbind(data.frame(pheno="AL", phase=1, mean.sd(boot.df$b1[boot.df$model==4 & boot.df$data.set=="AL"])),
               data.frame(pheno="AL", phase=2, mean.sd(boot.df$b2[boot.df$model==4 & boot.df$data.set=="AL"])),
               data.frame(pheno="AL", phase=3, mean.sd(boot.df$b3[boot.df$model==4 & boot.df$data.set=="AL"])),
               data.frame(pheno="AL", phase=4, mean.sd(boot.df$b4[boot.df$model==4 & boot.df$data.set=="AL"])),
               data.frame(pheno="RES", phase=1, mean.sd(boot.df$b1[boot.df$model==4 & boot.df$data.set=="RES"])),
               data.frame(pheno="RES", phase=2, mean.sd(boot.df$b2[boot.df$model==4 & boot.df$data.set=="RES"])),
               data.frame(pheno="RES", phase=3, mean.sd(boot.df$b3[boot.df$model==4 & boot.df$data.set=="RES"])),
               data.frame(pheno="RES", phase=4, mean.sd(boot.df$b4[boot.df$model==4 & boot.df$data.set=="RES"])),
               data.frame(pheno="Body size", phase=1, mean.sd(boot.df$b1[boot.df$model==1 & boot.df$data.set=="size21"])),
               data.frame(pheno="Cell size", phase=1, mean.sd(boot.df$b1[boot.df$model==1 & boot.df$data.set=="cellsize21"])),
               data.frame(pheno="Cell number", phase=1, mean.sd(boot.df$b1[boot.df$model==2 & boot.df$data.set=="cellnum21" & boot.df$tau1 > 2 & boot.df$tau1 < 10])),
               data.frame(pheno="Cell number", phase=2, mean.sd(boot.df$b2[boot.df$model==2 & boot.df$data.set=="cellnum21" & boot.df$tau1 > 2 & boot.df$tau1 < 10])), 
               data.frame(pheno="140", phase=1, mean.sd(boot.df$b1[boot.df$model==3 & boot.df$data.set=="140d_starved"])),
               data.frame(pheno="140", phase=2, mean.sd(boot.df$b2[boot.df$model==3 & boot.df$data.set=="140d_starved"])),
               data.frame(pheno="140", phase=3, mean.sd(boot.df$b3[boot.df$model==3 & boot.df$data.set=="140d_starved"])),
               data.frame(pheno="200", phase=1, mean.sd(boot.df$b1[boot.df$model==4 & boot.df$data.set=="200d_starved" & boot.df$tau1 > 2])),
               data.frame(pheno="200", phase=2, mean.sd(boot.df$b2[boot.df$model==4 & boot.df$data.set=="200d_starved" & boot.df$tau1 > 2])),
               data.frame(pheno="200", phase=3, mean.sd(boot.df$b3[boot.df$model==4 & boot.df$data.set=="200d_starved" & boot.df$tau1 > 2])),
               data.frame(pheno="200", phase=4, mean.sd(boot.df$b4[boot.df$model==4 & boot.df$data.set=="200d_starved" & boot.df$tau1 > 2]))
)


# summary plot of slope estimates
ggplot(t.4.df, aes(x = mean-1.96*sd, width=2*1.96*sd,
                   y=paste(pheno,phase), height=0.5,
                   fill=factor(phase), color=factor(phase))) + geom_tile() + 
  xlab("95% c.i. of estimated slope") + ylab("Experiment") + theme(legend.position = "none")

###### Table 5 -- estimates of changepoints from different experiments
t.5.df = rbind(data.frame(pheno="AL", phase=1, mean.sd(boot.df$tau1[boot.df$model==4 & boot.df$data.set=="AL"])),
               data.frame(pheno="AL", phase=2, mean.sd(boot.df$tau2[boot.df$model==4 & boot.df$data.set=="AL"])),
               data.frame(pheno="AL", phase=3, mean.sd(boot.df$tau3[boot.df$model==4 & boot.df$data.set=="AL"])),
               data.frame(pheno="RES", phase=1, mean.sd(boot.df$tau1[boot.df$model==4 & boot.df$data.set=="RES"])),
               data.frame(pheno="RES", phase=2, mean.sd(boot.df$tau2[boot.df$model==4 & boot.df$data.set=="RES"])),
               data.frame(pheno="RES", phase=3, mean.sd(boot.df$tau3[boot.df$model==4 & boot.df$data.set=="RES"])),
               data.frame(pheno="Cell size", phase=1, mean.sd(boot.df$tau1[boot.df$model==2 & boot.df$data.set=="cellsize21"])),
               data.frame(pheno="Cell number", phase=1, mean.sd(boot.df$tau1[boot.df$model==2 & boot.df$data.set=="cellnum21"])),
               data.frame(pheno="140", phase=1, mean.sd(boot.df$tau1[boot.df$model==3 & boot.df$data.set=="140d_starved"])),
               data.frame(pheno="140", phase=2, mean.sd(boot.df$tau2[boot.df$model==3 & boot.df$data.set=="140d_starved"])),
               data.frame(pheno="200", phase=1, mean.sd(boot.df$tau1[boot.df$model==4 & boot.df$data.set=="200d_starved"])),
               data.frame(pheno="200", phase=2, mean.sd(boot.df$tau2[boot.df$model==4 & boot.df$data.set=="200d_starved"])),
               data.frame(pheno="200", phase=3, mean.sd(boot.df$tau3[boot.df$model==4 & boot.df$data.set=="200d_starved"])) )

ggplot(t.5.df, aes(x = mean-1.96*sd, width=2*1.96*sd,
                   y=paste(pheno,phase), height=0.5,
                   fill=factor(phase), color=factor(phase))) + geom_tile() +
   xlab("95% c.i. of estimated changepoint") + ylab("Experiment") + theme(legend.position = "none")

####### Table 8 -- estimates of double times from various experiments
# removed "Body growth 0-21d phase 1"

t.8.df = rbind(data.frame(label="growth 10d",
                          mean.sd.lm(pull.stats(lm(log(mm2) ~ day, 
                                                   data=df.all[df.all$experiment=="10d_feeding",])))),
               
               data.frame(label="cell proliferation",
                          mean.sd.lm(pull.stats(lm(log(cellnumber) ~ day, 
                                                   data=df.new[df.new$experiment=="growth",])))),
               data.frame(label="cell growth",
                          mean.sd.lm(pull.stats(lm(log(cellsizemedian) ~ day, 
                                                   data=df.new[df.new$experiment=="growth",])))),
               data.frame(label="regrowth 170d",
                          mean.sd.lm(pull.stats(lm(log(mm2) ~ day, 
                                                   data=df.all[df.all$experiment=="170d_200d_starved_refed" & df.all$condition=="170d",])))),
               data.frame(label="regrowth 200d",
                          mean.sd.lm(pull.stats(lm(log(mm2) ~ day, 
                                                   data=df.all[df.all$experiment=="170d_200d_starved_refed" & df.all$condition=="200d",])))),
               data.frame(label="AL phase 1", 
                          mean.sd(boot.df$b1[boot.df$model==4 & boot.df$data.set=="AL"])),
               data.frame(label="AL phase 3", 
                          mean.sd(boot.df$b3[boot.df$model==4 & boot.df$data.set=="AL"])),
               data.frame(label="RES phase 1", 
                          mean.sd(boot.df$b1[boot.df$model==4 & boot.df$data.set=="RES"])),
               data.frame(label="RES phase 3", 
                          mean.sd(boot.df$b3[boot.df$model==4 & boot.df$data.set=="RES"])) )

# convert to doubling times
t.8.df$Td.lo = log(2)/(t.8.df$mean+1.96*t.8.df$sd)
t.8.df$Td.hi = log(2)/(t.8.df$mean-1.96*t.8.df$sd)

####### Table 9 -- estimates of half-lives. from various experiments
mean.140 = t.4.df$mean[(which(t.4.df$pheno=="140" & t.4.df$phase==2))]
sd.140 = t.4.df$sd[(which(t.4.df$pheno=="140" & t.4.df$phase==2))]
mean.200 = t.4.df$mean[(which(t.4.df$pheno=="200" & t.4.df$phase==2))]
sd.200 = t.4.df$sd[(which(t.4.df$pheno=="200" & t.4.df$phase==2))]

t.9.df = rbind(data.frame(label="degrowth 21d phase 2",
                          mean.sd.lm(pull.stats(lm(log(mm2) ~ day, 
                                                   data=df.new[df.new$experiment=="shrinkage_new",])))),
               
               data.frame(label="cell death 21d phase 2",
                          mean.sd.lm(pull.stats(lm(log(cellnumber) ~ day, 
                                                   data=df.new[df.new$experiment=="shrinkage_new",])))),
               
               data.frame(label="cell degrowth 21d phase 2",
                          mean.sd.lm(pull.stats(lm(log(cellsizemedian) ~ day, 
                                                   data=df.new[df.new$experiment=="shrinkage_new",])))),
               
               data.frame(label="body degrowth 140d", mean=mean.140, sd=sd.140),
               data.frame(label="body degrowth 200d", mean=mean.200, sd=sd.200),
               data.frame(label="AL phase 2", 
                          mean.sd(boot.df$b2[boot.df$model==4 & boot.df$data.set=="AL"])),
               data.frame(label="AL phase 4", 
                          mean.sd(boot.df$b4[boot.df$model==4 & boot.df$data.set=="AL"])),
               data.frame(label="RES phase 2", 
                          mean.sd(boot.df$b2[boot.df$model==4 & boot.df$data.set=="RES"])),
               data.frame(label="RES phase 4", 
                          mean.sd(boot.df$b4[boot.df$model==4 & boot.df$data.set=="RES"])) )

# convert to half-lives
t.9.df$Td.lo = log(2)/abs(t.9.df$mean-1.96*t.9.df$sd)
t.9.df$Td.hi = log(2)/abs(t.9.df$mean+1.96*t.9.df$sd)

####### Tables 11-13 -- Aiptasia data

t.aip.b.df = rbind(data.frame(pheno="Apo", phase=1, mean.sd(boot.df$b1[boot.df$model==2 & boot.df$data.set=="apo"])),
                 data.frame(pheno="Apo", phase=2, mean.sd(boot.df$b2[boot.df$model==2 & boot.df$data.set=="apo"])),
                 data.frame(pheno="Sym", phase=1, mean.sd(boot.df$b1[boot.df$model==2 & boot.df$data.set=="sym"])),
                 data.frame(pheno="Sym", phase=2, mean.sd(boot.df$b2[boot.df$model==2 & boot.df$data.set=="sym"])))
t.aip.tau.df = rbind(data.frame(pheno="Apo", phase=1, mean.sd(boot.df$tau1[boot.df$model==2 & boot.df$data.set=="apo"])),
                   data.frame(pheno="Sym", phase=1, mean.sd(boot.df$tau1[boot.df$model==2 & boot.df$data.set=="sym"])))

t.11.df = t.aip.b.df[t.aip.b.df$phase == 1,]
t.11.df$Td.lo = log(2)/t.11.df$mean-1.96*t.11.df$sd
t.11.df$Td.hi = log(2)/t.11.df$mean+1.96*t.11.df$sd

t.12.df = t.aip.b.df[t.aip.b.df$phase == 2,]
t.12.df$Td.lo = log(2)/abs(t.12.df$mean-1.96*t.12.df$sd)
t.12.df$Td.hi = log(2)/abs(t.12.df$mean+1.96*t.12.df$sd)

t.13.df = t.aip.tau.df
t.13.df$lo = t.13.df$mean-1.96*t.13.df$sd
t.13.df$hi = t.13.df$mean+1.96*t.13.df$sd
               
############ pull everything together to an Excel sheet
# build a dataframe for the contents page
contents.df = data.frame()
contents.df = rbind(contents.df, 
                    data.frame(sheet="Contents", description="this contents page"),
                    data.frame(sheet="Table.1", description="summary statistics of linear models with time for feeding and refeeding experiments"),
                    data.frame(sheet="Table.2", description="summary statistics of multiple linear regression for body size, cell size, cell number"),
                    data.frame(sheet="Table.3", description="summary statistics for simple linear regression for body size, cell size, cell number"),
                    data.frame(sheet="Table.4", description="summary statistics of slopes with time in multi-phase linear models"),
                    data.frame(sheet="Table.5", description="summary statistics of changepoint times in multi-phase linear models"),
                    data.frame(sheet="Table.6", description="summary statistics of multiple linear regression for body size, cell size, cell number, partitioned by time period"),
                    data.frame(sheet="Table.7", description="summary statistics of simple linear regression for body size, cell size, cell number, partitioned by time period"),
                    data.frame(sheet="Table.8", description="doubling times for growth models in different experiments"),
                    data.frame(sheet="Table.9", description="half-lives for decay models in different experiments"),
                    data.frame(sheet="Table.10", description="loss rates for decay models, partitioned by time period"),
                    data.frame(sheet="Table.11", description="growth rates for Aiptasia"),
                    data.frame(sheet="Table.12", description="decay rates for Aiptasia"),
                    data.frame(sheet="Table.13", description="changepoints for Aiptasia") )

# write the contents and all the data tables
write_xlsx(list(Contents = contents.df,
                Table.1 = t.1.df,
                Table.2 = t.2.df,
                Table.3 = t.3.df,
                Table.4 = t.4.df,
                Table.5 = t.5.df,
                Table.6 = t.6.df,
                Table.7 = t.7.df,
                Table.8 = t.8.df,
                Table.9 = t.9.df,
                Table.10 = t.10.df,
                Table.11 = t.11.df,
                Table.12 = t.12.df,
                Table.13 = t.13.df), "outputs.xlsx")

############ reproduce various summary plots for checking
###### summary plot in Fig S2D
focus = c("growth 10d", "AL phase 1", 
          "RES phase 1", "AL phase 3", 
          "RES phase 3", "regrowth 200d", 
          "regrowth 170d")
plot.s2d.df = t.8.df[which(t.8.df$label %in% focus),]
ggplot(plot.s2d.df, aes(x = mean-1.96*sd, width=2*1.96*sd,
                        y=label, height=0.5,
                        fill=label, color=label)) + geom_tile() +
  xlab("95% c.i. of estimated slope") + ylab("Experiment") + theme(legend.position="none")

###### summary plot in Fig S3H
focus = c("degrowth 21d phase 2", "body degrowth 140d", 
          "body degrowth 200d", "AL phase 2", 
          "RES phase 2", "AL phase 4", 
          "RES phase 4")
plot.s3h.df = t.9.df[which(t.9.df$label %in% focus),]
ggplot(plot.s3h.df, aes(x = mean-1.96*sd, width=2*1.96*sd,
                        y=label, height=0.5,
                        fill=label, color=label)) + geom_tile() +
  xlab("95% c.i. of estimated slope") + ylab("Experiment")+ theme(legend.position="none")

###### summary plot in Fig S4O
# big change since first instance of paper
plot.s4o.df = rbind(data.frame(pheno="Body size", phase=1, mean.sd(boot.df$b1[boot.df$model==1 & boot.df$data.set=="size21"])),
                                   data.frame(pheno="Cell size", phase=1, mean.sd(boot.df$b1[boot.df$model==1 & boot.df$data.set=="cellsize21"])),
                                   data.frame(pheno="Cell number", phase=1, mean.sd(boot.df$b1[boot.df$model==2 & boot.df$data.set=="cellnum21" & boot.df$tau1 > 2 & boot.df$tau1 < 10])),
                                   data.frame(pheno="Cell number", phase=2, mean.sd(boot.df$b2[boot.df$model==2 & boot.df$data.set=="cellnum21" & boot.df$tau1 > 2 & boot.df$tau1 < 10]))
)

plot.s4o.df$phase = factor(plot.s4o.df$phase)
ggplot(plot.s4o.df, aes(x = mean-1.96*sd, width=2*1.96*sd,
                        y=paste(pheno, phase), height=0.5,
                        fill=phase, color=phase)) + geom_tile() +
  xlab("95% c.i. of estimated slope") + ylab("Experiment")+ theme(legend.position="none")

