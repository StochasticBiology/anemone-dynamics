library(ggplot2)
library(gridExtra)

df = read.csv2("NEWdata_shrinkage_new_shrinkage_old_growth.csv")

### relationships between cell size, cell number, body size

lm.csize.cnumber.old = lm(log(df$cellnumber[df$experiment=="shrinkage_old"]) ~ 
                        log(df$cellsizemedian[df$experiment=="shrinkage_old"]))
lm.csize.cnumber.new = lm(log(df$cellnumber[df$experiment=="shrinkage_new"]) ~ 
                          log(df$cellsizemedian[df$experiment=="shrinkage_new"]))

lm.bsize.cnumber.old = lm(log(df$mm2[df$experiment=="shrinkage_old"]) ~ 
                            log(df$cellnumber[df$experiment=="shrinkage_old"]))
lm.bsize.cnumber.new = lm(log(df$mm2[df$experiment=="shrinkage_new"]) ~ 
                            log(df$cellnumber[df$experiment=="shrinkage_new"]))

lm.bsize.csize.old = lm(log(df$mm2[df$experiment=="shrinkage_old"]) ~ 
                            log(df$cellsizemedian[df$experiment=="shrinkage_old"]))
lm.bsize.csize.new = lm(log(df$mm2[df$experiment=="shrinkage_new"]) ~ 
                            log(df$cellsizemedian[df$experiment=="shrinkage_new"]))

lm.bsize.cnumber.growth = lm(log(df$mm2[df$experiment=="growth"]) ~ 
                            log(df$cellnumber[df$experiment=="growth"]))
lm.bsize.csize.growth = lm(log(df$mm2[df$experiment=="growth"]) ~ 
                          log(df$cellsizemedian[df$experiment=="growth"]))

summary(lm.bsize.csize.growth) # for growth, log body size ~ 3.25 log cell size (***)
summary(lm.bsize.cnumber.growth) # for growth, log body size ~ 0.73 log cell number (***)

summary(lm.bsize.csize.old) # for old shrinkage, log body size ~ 0.38 log cell size (NS)
summary(lm.bsize.csize.new) # for new shrinkage, log body size ~ 0.83 log cell size (***)

summary(lm.bsize.cnumber.old) # for old shrinkage, log body size ~ 0.14 log cell number (***)
summary(lm.bsize.cnumber.new) # for new shrinkage, log body size ~ 0.23 log cell number (***)

g.number = ggplot(df, aes(x=log(cellnumber), y=log(mm2), color=factor(day))) + geom_point() + facet_wrap(~experiment)
g.size = ggplot(df, aes(x=log(cellsizemedian), y=log(mm2), color=factor(day))) + geom_point() + facet_wrap(~experiment)

#### difference between old and new experiments?

shrinkage.df = df[df$experiment != "growth",]

lm.bsize.cnumber.both = lm(log(shrinkage.df$mm2) ~ log(shrinkage.df$cellnumber)*shrinkage.df$experiment)
summary(lm.bsize.cnumber.both)
# intercept difference p < 0.05, slope difference p > 0.05

lm.bsize.csize.both = lm(log(shrinkage.df$mm2) ~ log(shrinkage.df$cellsizemedian)*shrinkage.df$experiment)
summary(lm.bsize.csize.both)
# intercept and slope p > 0.05

#### time series

lm.bsize.t.growth = lm(log(df$mm2[df$experiment=="growth"]) ~ 
                         df$day[df$experiment=="growth"])
lm.bsize.t.old = lm(log(df$mm2[df$experiment=="shrinkage_new"]) ~ 
                         df$day[df$experiment=="shrinkage_new"])
lm.bsize.t.new = lm(log(df$mm2[df$experiment=="shrinkage_old"]) ~ 
                         df$day[df$experiment=="shrinkage_old"])

summary(lm.bsize.t.growth) # for growth, log body size ~ 0.27 * time (***)
g.timeseries =  ggplot(df, aes(x=day, y=log(mm2), color=factor(day))) + geom_point() + facet_wrap(~experiment)

df$igjlabel = ""
for(i in 1:nrow(df)) {
  ss = strsplit(df$individual[i], split="_")[[1]]
  df$igjlabel[i] = ss[length(ss)]
}
g.timeserieslines = ggplot(df, aes(x=day,y=log(mm2),color=igjlabel)) + geom_line(alpha=0.3) + geom_point() + facet_wrap(~experiment)

# take a look at growth rates as function of body size
df$rate = NA
expts = unique(df$experiment)
for(expt in expts) {
  for(i in 1:nrow(df)) {
  if(df$experiment[i] == expt)
  {
    this.one = df$igjlabel[i]
    this.one.refs = which(df$experiment == expt & df$igjlabel == this.one)
    future.refs = this.one.refs[which(this.one.refs > i)]
    if(length(future.refs) > 0) {
       future.ref = future.refs[1]
       sizediff = df$mm2[future.ref]-df$mm2[i]
       timediff = df$day[future.ref]-df$day[i]
       df$rate[i] = sizediff/timediff
    }
  }
  }
}

g.rates = ggplot(df, aes(x=mm2,y=rate,color=factor(day))) + geom_point() + facet_wrap(~experiment)

lm.rate.bsize.growth = lm(df$rate[df$experiment=="growth"] ~ 
                         df$mm2[df$experiment=="growth"])
lm.rate.bsize.old = lm(df$rate[df$experiment=="shrinkage_new"] ~ 
                      df$mm2[df$experiment=="shrinkage_new"])
lm.rate.bsize.new = lm(df$rate[df$experiment=="shrinkage_old"] ~ 
                      df$mm2[df$experiment=="shrinkage_old"])

summary(lm.rate.bsize.growth) # rate ~ 0.34 body size (***)
summary(lm.rate.bsize.old) # rate ~ -0.33 body size (***)
summary(lm.rate.bsize.new) # rate ~ -0.51 body size (***)

sf = 2
png("summary-plots.png", width=800*sf, height=1000*sf, res=72*sf)
grid.arrange(g.timeserieslines, g.number, g.size, g.rates, nrow=4)
dev.off()



