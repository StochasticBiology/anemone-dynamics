#imports
library(ggplot2)
library(gridExtra)

#FEEDING
#open data files .csv for feeding periods 0-10d
#df.growth contains body size, cell number and cell size
#df.growth.body contains only body size (with more data points than the data frame above)
df.growth = read.csv("Data/nematostella_growth_bgp.csv")
df.growth.body=read.csv("Data/body_size_growth_bgp.csv")

#fit simple linear regressions for log body size, log cell number and log cell size  with respect to day
my.lm.bodysize = lm(log.body.size ~ Day, data=df.growth.body)
summary(my.lm.bodysize) #summary statistics of the linear model
AIC(my.lm.bodysize)#compute the AIC value for future model selection

my.lm.cellnumber = lm(log.cell.number ~ Day_y, data=df.growth)
summary(my.lm.cellnumber)
AIC(my.lm.cellnumber)

my.lm.cellsize = lm(log.cell.size ~ Day, data=df.growth)
summary(my.lm.cellsize)
AIC(my.lm.cellsize)

#plots of the simple linear regressions

ggplot(df.growth.body, aes(x = Day, y = log.body.size)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "darkblue")+geom_point(color='cadetblue3')+
  labs(title ="(A) Feeding 0-10d", y='Body size log')+
  theme(plot.title = element_text(face="bold", size=16))

ggplot(df.growth, aes(x = Day_y, y = log.cell.number)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "darkblue")+geom_point(color='cadetblue3')+
  labs(y='Cell number log', x='Day')

ggplot(df.growth, aes(x = Day, y = log.cell.size)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "darkblue")+geom_point(color='cadetblue3')+
  labs(y='Cell size log')


#Fit a multiple linear regression model to find correlations between log body size and log cell number and size
my.lm.growth = lm(log.body.size ~ log.cell.size + log.cell.number, data=df.growth)
summary(my.lm.growth)
AIC(my.lm.growth)

#Fit two simple linear regressions to find correlations between body size and cell size, and body size
#and cell number
my.lm.growth.size = lm(log.body.size ~ log.cell.size, data=df.growth)
summary(my.lm.growth.size)
AIC(my.lm.growth.size)

my.lm.growth.numb = lm(log.body.size ~ log.cell.number, data=df.growth)
summary(my.lm.growth.numb)
AIC(my.lm.growth.numb)

#plots of the independent correlations
ggplot(df.growth, aes(x = log.cell.size, y = log.body.size, col=factor(Day))) + 
  geom_point()+
  labs(title = paste("Adj R2 = ",signif(summary(my.lm.growth.size)$adj.r.squared, 5),
                     "Intercept =",signif(my.lm.growth.size$coef[[1]],5 ),
                     " Slope =",signif(my.lm.growth.size$coef[[2]], 5),
                     " P =",signif(summary(my.lm.growth.size)$coef[2,4], 5)), y='Body size log',
       x='Cell size log')

ggplot(df.growth, aes(x = log.cell.number, y = log.body.size, col=factor(Day))) + 
  geom_point()+
  labs(title = paste("Adj R2 = ",signif(summary(my.lm.growth.numb)$adj.r.squared, 5),
                     "Intercept =",signif(my.lm.growth.numb$coef[[1]],5 ),
                     " Slope =",signif(my.lm.growth.numb$coef[[2]], 5),
                     " P =",signif(summary(my.lm.growth.numb)$coef[2,4], 5)), y='Body size log',
       x='Cell number log')

#STARVATION
#open data file containing body size, cell number and cell size for starvation periods 0-21d
df.degrowth = read.csv("Data/nematostella_degrowth_bgp.csv")
df.degrowth02=subset(df.degrowth, Day <=2) #get first part where growth happens
df.degrowth221 = subset(df.degrowth, Day>= 2) #get second part where degrowth happens

#Fit a multiple linear regression model to find correlations between log body size and log cell number and size
#for growth and degrowth phases in starvation
my.lm.degrowth02 = lm(log.body.size~ log.cell.size + log.cell.number, data=df.degrowth02)
summary(my.lm.degrowth02)
AIC(my.lm.degrowth02)

my.lm.degrowth221 = lm(log.body.size ~ log.cell.size + log.cell.number, data=df.degrowth221)
summary(my.lm.degrowth221)
AIC(my.lm.degrowth221)


#Fit simple linear regressions to find correlations between body size and cell size, and body size
#and cell number, separately

#growth phase
my.lm.degrowth02.size = lm(log.body.size ~ log.cell.size, data=df.degrowth02)
summary(my.lm.degrowth02.size)
AIC(my.lm.degrowth02.size)

my.lm.degrowth02.numb = lm(log.body.size ~ log.cell.number, data=df.degrowth02)
summary(my.lm.degrowth02.numb)
AIC(my.lm.degrowth02.numb)

#plots
ggplot(df.degrowth02, aes(x = log.cell.size, y = log.body.size, col=factor(Day))) + 
  geom_point()+
  labs(title = paste("Adj R2 = ",signif(summary(my.lm.degrowth02.size)$adj.r.squared, 5),
                     "Intercept =",signif(my.lm.degrowth02.size$coef[[1]],5 ),
                     " Slope =",signif(my.lm.degrowth02.size$coef[[2]], 5),
                     " P =",signif(summary(my.lm.degrowth02.size)$coef[2,4], 5)), y='Body size log',
       x='Cell size log')



ggplot(df.degrowth02, aes(x = log.cell.number, y = log.body.size, col=factor(Day))) + 
  geom_point()+
  labs(title = paste("Adj R2 = ",signif(summary(my.lm.degrowth02.numb)$adj.r.squared, 5),
                     "Intercept =",signif(my.lm.degrowth02.numb$coef[[1]],5 ),
                     " Slope =",signif(my.lm.degrowth02.numb$coef[[2]], 5),
                     " P =",signif(summary(my.lm.degrowth02.numb)$coef[2,4], 5)), y='Body size log',
       x='Cell number log')

#degrowth phase
my.lm.degrowth221.size = lm(log.body.size ~ log.cell.size, data=df.degrowth221)
summary(my.lm.degrowth221.size)
AIC(my.lm.degrowth221.size)


my.lm.degrowth221.numb = lm(log.body.size ~ log.cell.number, data=df.degrowth221)
summary(my.lm.degrowth221.numb)
AIC(my.lm.degrowth221.numb)

#plots
ggplot(df.degrowth221, aes(x = log.cell.size, y = log.body.size, col=factor(Day))) + 
  geom_point()+
  labs(title = paste("Adj R2 = ",signif(summary(my.lm.degrowth221.size)$adj.r.squared, 5),
                     "Intercept =",signif(my.lm.degrowth221.size$coef[[1]],5 ),
                     " Slope =",signif(my.lm.degrowth221.size$coef[[2]], 5),
                     " P =",signif(summary(my.lm.degrowth221.size)$coef[2,4], 5)), y='Body size log',
       x='Cell size log')

ggplot(df.degrowth221, aes(x = log.cell.number, y = log.body.size, col=factor(Day))) + 
  geom_point()+
  labs(title = paste("Adj R2 = ",signif(summary(my.lm.degrowth221.numb)$adj.r.squared, 5),
                     "Intercept =",signif(my.lm.degrowth221.numb$coef[[1]],5 ),
                     " Slope =",signif(my.lm.degrowth221.numb$coef[[2]], 5),
                     " P =",signif(summary(my.lm.degrowth221.numb)$coef[2,4], 5)), y='Body size log',
       x='Cell number log')


#REFEEDING
#Open data files for the two different treatments, 170d and 200d
df.treatments = read.csv("Data/treatments.csv") #data points for 170d and 200d merged together
df.treatment_170d = read.csv("Data/treatment_170d.csv") #data points for 170d
df.treatment_200d = read.csv("Data/treatment_200d.csv") #data points for 200d

#Fit simple linear regression models
my.lm.treatments = lm(log.body.size ~ Day, data=df.treatments)
summary(my.lm.treatments)
AIC(my.lm.treatments)

my.lm.treatment_170d = lm(log.body.size ~ Day, data=df.treatment_170d)
summary(my.lm.treatment_170d)
AIC(my.lm.treatment_170d)

my.lm.treatment_200d = lm(log.body.size ~ Day, data=df.treatment_200d)
summary(my.lm.treatment_200d)
AIC(my.lm.treatment_200d)

#plot
p=ggplot(df.treatment_170d, aes(x = Day, y = log.body.size, color="170d")) + 
  geom_point(col='orange') +
  stat_smooth(method = "lm", aes(color="170d"))+
  labs(title='Re-feeding', y='Body size log', x='Day') 
p=p+geom_point(data = df.treatment_200d, aes(color="200d"))+
  stat_smooth(data = df.treatment_200d, method = "lm", aes(color="200d"))
p=p+
  stat_smooth(data = df.treatments, method = "lm", aes(color="both"))+
  scale_color_manual(name = "", values = c("170d" = "orange", 
                                           "200d" = "purple", "both"="black"))
print(p)

#plot removing the linear regression model for the merged dataset 
p=ggplot(df.treatment_170d, aes(x = Day, y = log.body.size, color="170 days")) + 
  geom_point(col='orange') +
  stat_smooth(method = "lm", aes(color="170 days"))+
  labs(title='(B) Refeeding', y='Body size log', x='Day')+
  theme(plot.title = element_text(face="bold", size=16)) 
p=p+geom_point(data = df.treatment_200d, aes(color="200 days"))+
  stat_smooth(data = df.treatment_200d, method = "lm", aes(color="200 days"))+
  scale_color_manual(name = "After starved for:", values = c("170 days" = "orange", 
                                                             "200 days" = "purple"))
print(p)

#long starvation data
df.longB = read.csv("Data/long_starvation_B.csv")
df.longanB = read.csv("Data/long_starvation_analogous_B.csv")

my.lm.longB = lm(log.body.size ~ Day, data=df.longB)
summary(my.lm.longB)
sigma(my.lm.longB)
AIC(my.lm.longB)

my.lm.longanB = lm(log.body.size ~ Day, data=df.longanB)
summary(my.lm.longanB)
AIC(my.lm.longanB)

#plot of long B
ggplot(df.longB, aes(x = Day, y = log.body.size)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "darkblue")+geom_point(color='cadetblue3')+
  labs(title ="Long starvation (i)", y='Body size log')+
  theme(plot.title = element_text(face="bold", size=16))

#Comparison of growth and degrowth data by merging dataframes,respectevly, and analysing if a nested
#model for the whole data is better than several models, one for each dataset
df.all_growth_02 = read.csv("Data/all_growth_02.csv") #all growth lines including simple regression starvation 0-2d
df.all_degrowth_221=read.csv("Data/all_degrowth_221.csv") #all degrowth lines including simple regression starvation 2-21d

my.lm.g_02 = lm(log.body.size ~ Day, data=df.all_growth_02)
summary(my.lm.g_02)
AIC(my.lm.g_02)

my.lm.dg_221 = lm(log.body.size ~ Day, data=df.all_degrowth_221)
summary(my.lm.dg_221)
AIC(my.lm.dg_221)

#perform likelihood ratio test
library(lmtest)

#growth
nested1=update(my.lm.bodysize, data=my.lm.g_02$model) 
lrtest(my.lm.g_02, nested1) #nested model vs body size 0-10d
nested2=update(my.lm.treatment_170d, data=my.lm.g_02$model)
lrtest(my.lm.g_02, nested2) #nested model vs 170d
nested3=update(my.lm.treatment_200d, data=my.lm.g_02$model)
lrtest(my.lm.g_02, nested3) #nested model vs 200d
nested4=update(my.lm.bodysize_02, data=my.lm.g_02$model)
lrtest(my.lm.g_02, nested4) #nested model vs 0-2d

#degrowth 2-21

#first find the linear regression models missing for degrowth 2-21d body~Day 
my.lm.bodysize_221 = lm(log.body.size ~ Day, data=df.degrowth221)
summary(my.lm.degrowth221.size)
AIC(my.lm.degrowth221.size)

nest1=update(my.lm.dg_221, data=my.lm.bodysize_221$model)
lrtest(nest1, my.lm.bodysize_221) #nested model vs body size 2-21d
nest2=update(my.lm.dg_221, data=my.lm.longB$model)
lrtest(nest2, my.lm.longB) #nested model vs long starvation (i)
nest3=update(my.lm.dg_221, data=my.lm.longanB$model)
lrtest(nest3, my.lm.longanB) #nested model vs long starvation (ii)


#Computation of the log likelihood functions of each of the models
logLik(my.lm.g_02)
logLik(my.lm.bodysize)
logLik(my.lm.treatment_170d)
logLik(my.lm.treatment_200d)
logLik(my.lm.bodysize)
logLik(my.lm.longB)
logLik(my.lm.longanB)
logLik(my.lm.dg_221)

