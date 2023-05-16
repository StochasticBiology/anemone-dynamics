# model = (number of phases)
# parameters theta = { b1, c, sigma, b2, tau1, b3, tau2, b4,tau3 }
# note that c is the intercept b0 in the text documents

#definition of the negative log likelihood function
neg.log.likelihood = function(theta, model, y, x)
{
  # parse argument to parameters
  b1 = theta[1]
  c = theta[2]
  sigma = theta[3]
  if(model > 1) {#more than 1 phase
    b2 = theta[4]
    tau1 = theta[5]
    if(model > 2) {#more than 2 phases
      b3 = theta[6]
      tau2 = theta[7]
      if(model > 3) {#more than 3 phases
        b4 = theta[8]
        tau3 = theta[9]
      } else { tau3 = Inf }#two changepoints, three phases
    } else { tau2 =tau3= Inf }#just one changepoint, 2 phases
  } else { tau1 = tau2 = tau3= Inf }#no changepoints, just one phase
  
  # initialise log lik
  llik = 0
  for(i in 1:length(x)) {
    if(x[i] < tau1) {
      llik = llik + dnorm(y[i], mean=b1*x[i]+c, sd = sigma, log=T)
    } else if(x[i] < tau2) {
      llik = llik + dnorm(y[i], mean=b2*(x[i]-tau1)+(b1*tau1+c), sd = sigma, log=T) 
    } else if(x[i] < tau3) {
      llik = llik + dnorm(y[i], mean=b3*(x[i]-tau2)+(b2*(tau2-tau1) + b1*tau1+c), sd = sigma, log=T) 
    } else {
      llik = llik + dnorm(y[i], mean=b4*(x[i]-tau3)+(b3*(tau3-tau2) +b2*(tau2-tau1)+ b1*tau1+c), sd = sigma, log=T) 
    }
  }
  return(-llik)
}


# retrieve the mean value from a parameterised model
# arguments as above
model.val = function(theta, model, x) {
  # parse argument to parameters
  b1 = theta[1]
  c = theta[2]
  sigma = theta[3]
  if(model > 1) {#more than 1 phase
    b2 = theta[4]
    tau1 = theta[5]
    if(model > 2) {#more than 2 phases
      b3 = theta[6]
      tau2 = theta[7]
      if(model > 3) {#more than 3 phases
        b4 = theta[8]
        tau3 = theta[9]
      } else { tau3 = Inf }#two changepoints, three phases
    } else { tau2 =tau3= Inf }#just one changepoint, 2 phases
  } else { tau1 = tau2 = tau3= Inf }#no changepoints, just one phase
  
  
  # simply read out mean prediction from changepoint model
  if(x < tau1) { return(b1*x+c) }
  else if(x < tau2) { return(b2*(x-tau1) + b1*tau1 +c) }
  else if(x < tau3) { return(b3*(x-tau2) + b2*(tau2-tau1) + b1*tau1 +c) }
  else { return(b4*(x-tau3) + b3*(tau3-tau2) + b2*(tau2-tau1)+ b1*tau1+c)}
}



# plot data and/or model fit
# new.plot = 1 (start new plot, draw data and AIC); -1 (redraw points on existing plot); 0 (draw lines on existing plot)
# (use for bootstrapping output)
plotfn = function(opt, model, y, x, new.plot=1) {
  # pull parameters from opt output
  theta = opt$par
  # create new set of x values and compute corresponding ys
  xmodel = floor(min(x)):ceiling(max(x))
  ymodel = unlist(lapply(xmodel, model.val, theta=theta, model=model))  
  # DOFs for AIC calculation
  if(model == 1) { dof = 3}
  if(model == 2) { dof = 5}
  if(model == 3) { dof = 7}
  if(model == 4) { dof = 9}
  # plot
  if(new.plot == 1) { plot(x,y, col = "cadetblue3", xlab="Day",
                           ylab="log(y)", pch = 16,
                           main=paste0("AIC=", round(2*opt$val+2*dof, digits = 2)),
                           cex.main=1, cex.lab=0.75, cex.axis=0.9, cex.sub=0.75,
                           abline(v=x,col="gray93",lwd=300), las = 1)}
  else if(new.plot == -1) { points(x,y) }
  else {lines(xmodel,ymodel, col="darkblue", las = 1, lwd=1)} #change lwd=3 for just 1 boot, thicker line
}

# read data
df.1 = read.csv("Data/long_starvation_analogous_B.csv")
df.2 = read.csv("Data/nematostella_degrowth_bgp.csv")
df.AL=read.csv("Data/updated_AL.csv")
df.RES=read.csv("Data/updated_RES.csv")
df.apo=read.csv("Data/apo.csv")
df.symb=read.csv("Data/symb.csv")

# set up trellis plot
par(mfrow=c(3,4))
par(mar=c(3, 2, 2, 1), mgp=c(3, 1, 0), las=1)


#lists of slopes, get them for a given model and a given data frame 
#WARNING: Run the code for just one expt and one model.type, then merge the data of each slopes together (done in the .ipynb file)
#comment # b4_l and/or b3_l for 3-phases or 2-phases models
b1_l=list()
b2_l=list()
b3_l=list()
b4_l=list()

# loop through experiment
for(expt in 1:8) {
  # noise to put on top of x values
  jitter = 0.
  # pull data out of frames depending on experiment
  if(expt == 1) { x = df.1$Day + rnorm(nrow(df.1),sd=jitter); y = df.1$log.body.size }
  if(expt == 2) { x = df.2$Day+ rnorm(nrow(df.2),sd=jitter); y = df.2$log.body.size }
  if(expt == 3) { x = df.2$Day+ rnorm(nrow(df.2),sd=jitter); y = df.2$log.cell.number}
  if(expt == 4) { x = df.2$Day+ rnorm(nrow(df.2),sd=jitter); y = df.2$log.cell.size}
  if(expt == 5) { x = df.AL$day+ rnorm(nrow(df.AL),sd=jitter); y = df.AL$log.body.size}
  if(expt == 6) { x = df.RES$day+ rnorm(nrow(df.RES),sd=jitter); y = df.RES$log.body.size}
  if(expt == 7) { x = df.apo$day+ rnorm(nrow(df.apo),sd=jitter); y = df.apo$log.body.size}
  if(expt == 8) { x = df.symb$day+ rnorm(nrow(df.symb),sd=jitter); y = df.symb$log.body.size}
  
  # loop through changepoint numbers
  for(model.type in 1:4)
  { 
    # bootstrap
    for(nboot in 1:3) #1:100 for plotting less messy figures
    {
      # for the first sample, retain original data, otherwise sample with replacement
      if(nboot == 1) {  boot = 1:length(y) }
      else {  boot = sample(length(y), length(y), replace=T) }
      bx = x[boot]
      by = y[boot]
      
      # initial guess -- flat, sigma = 1, changepoints at 1/4, 1/2 and 3/4 through x
      theta = c(0, mean(by[bx==min(bx)]), 1, 0, max(bx)/3, 0, 2*max(bx)/3, 0, 8*max(bx)/9) 
      if(expt == 1 & model.type != 3) { theta = c(0.5, mean(by[bx==min(bx)]), 1, 0,75, 0, 80, 0,180) }
      if(expt == 1 & model.type == 3) { theta = c(0.5, mean(by[bx==min(bx)]), 1, 0, 5, 0, 80, 0,180) }
      if(expt == 2) { theta = c(0, mean(by[bx==min(bx)]), 1, 0, 3, 0, 11, 0, 18) }
      if(expt == 3) { theta = c(0, mean(by[bx==min(bx)]), 1, 0, 3, 0, 15, 0, 18) }
      if(expt == 4) { theta = c(0, mean(by[bx==min(bx)]), 1, 0, 3, 0, 15, 0, 18) }
      if(expt == 5) { theta = c(0, mean(by[bx==min(bx)]), 1, 0, 5, 0, 90, 0, 120) }
      if(expt == 6) { theta = c(0, mean(by[bx==min(bx)]), 1, 0, 10, 0, 70, 0, 120) }
      if(expt == 7) { theta = c(0, mean(by[bx==min(bx)]), 1, 0, 14, 0, 40, 0, 100) }
      if(expt == 8) { theta = c(0, mean(by[bx==min(bx)]), 1, 0, 14, 0, 40, 0, 100) }
      
      # attempt to optimise
      opt = optim(theta, neg.log.likelihood, model=model.type, y=by, x=bx)
      
      #append slopes to lists, this only for saving and plotting the model outputs of distributions of slopes
      #NOTE that working with these lists bi_l is only for one expt and one model.type, 
      #for running the entire algorithm, comment them #
      b1_l=append(b1_l, opt$par[1])
      b2_l=append(b2_l, opt$par[4])
      b3_l=append(b3_l, opt$par[6])
      b4_l=append(b4_l, opt$par[8])
      
      # plot data and AIC for first "sample", otherwise fitted line
      if(nboot == 1) { new.plot = 1}
      else {new.plot = 0}
      plotfn(opt, model.type, by, bx, new.plot)
    }
  }
}


#put together into a .csv file the slope estimates from bootstrapping. One .csv file for each model
#remark that this is only done when we have run the loop above for a single expt and model.type
#unlist the lists to create a data frame
lst_b1=unlist(b1_l, use.names = FALSE)
lst_b2 =unlist(b2_l, use.names = FALSE)
lst_b3=unlist(b3_l, use.names = FALSE)
lst_b4 =unlist(b4_l, use.names = FALSE)

#name the data frame by dfEXPT_mMODEL.TYPE, for the choice of EXPT from 1 to 6 and MODEL.TYPE from 1 to 4
#for example:
future_csv_df6_m4=data.frame("b_1 RES"=c(lst_b1), "b_2 RES"=c(lst_b2), 
                             "b_3 RES"=c(lst_b3),"b_4 RES"=c(lst_b4))
#save as .csv file
write.csv(future_csv_df6_m4, "Data/df6_m4.csv")

#Open data file with all the slope estimates for all models
#this data file is constructed using Python code in the .ipynb file
df = read.csv("Data/slopes_boots_newformat_updated.csv")

#Select all columns for body size, removing like this the data for cell measurements
library(dplyr)
library(tidyverse)

df_body= df%>% filter(Legend == c( "b_2.0.21d","b_2.long.starvation..ii.","b_2.AL",
                                   "b_4.AL",  "b_2.RES", "b_4.RES","b_1.long.starvation..ii.",
                                   "b_3.long.starvation..ii.","b_1.AL", "b_3.AL", "b_1.RES",
                                   "b_3.RES","b_1.0.21d"))

#order the columns so that in the plot the negative and positive slopes appear clustered
df_ordered_clust=df_body%>% 
  mutate(Legend = factor(Legend, levels = c( "b_2.0.21d","b_2.long.starvation..ii.","b_2.AL",
                                             "b_4.AL",  "b_2.RES", "b_4.RES","b_1.long.starvation..ii.",
                                             "b_3.long.starvation..ii.","b_1.AL", "b_3.AL", "b_1.RES","b_3.RES")))

#plot of the slopes distributions, in the same x-axis, Nematostella data
library(ggplot2)
library(ggridges)


ggplot(df_ordered_clust, aes(x=Boots.value, y =Legend, fill =Legend)) +
  geom_density_ridges() +
  theme_ridges(font_size = 13, grid = TRUE) +
  labs(title = 'Distribution of slopes from bootstrapping',x='Value') +
  theme(legend.position = "none", axis.title.y = element_blank())


#Select columns for cell level measurements
df_cell= df%>% filter(Legend == c( "b_2.0.21d.cell.number", "b_2.0.21d.cell.size","b_3.0.21d.cell.size",
                                   "b_1.0.21d.cell.number", "b_1.0.21d.cell.size"))

#order the columns so that in the plot the negative and positive slopes appear clustered
df_ordered_clust=df_cell %>% 
  mutate(Legend = factor(Legend, levels = c( "b_2.0.21d.cell.number", "b_2.0.21d.cell.size","b_3.0.21d.cell.size",
                                             "b_1.0.21d.cell.number", "b_1.0.21d.cell.size"))) %>% arrange(Legend)

#plot of the slopes distributions, in the same x-axis
ggplot(df_ordered_clust, aes(x=Boots.value, y =Legend, fill =Legend)) +
  geom_density_ridges() +
  theme_ridges(font_size = 13, grid = TRUE) +
  xlim(-0.25, 0.75)+
  labs(
    title = 'Distribution of slopes from bootstrapping',x='Value') +
  theme(legend.position = "none", axis.title.y = element_blank())



#Construct confidence intervals from bootstrapping with a 95% confidence

#read data file with all slopes estimates together, the construction of this data file is in the .ipynb file
df = read.csv("Data/slopes_boots_updated.csv")
df_aiptasia= read.csv("Data/slopes_apo_symb.csv")
  
#select each column, which corresponds to a slope, and compute the values at the 2.5th and 97th percentile
#these two values are the endpoints of the confidence interval of that slope
n_body_b1=quantile(df$b_1.0.21d, na.rm = T,probs = c(0.025,0.975))
n_body_b2=quantile(df$b_2.0.21d, na.rm = T,probs = c(0.025,0.975))
n_cell_numb_b1=quantile(df$b_1.0.21d.cell.number, na.rm = T,probs = c(0.025,0.975))
n_cell_numb_b2=quantile(df$b_2.0.21d.cell.number, na.rm = T,probs = c(0.025,0.975))
n_cell_size_b1=quantile(df$b_1.0.21d.cell.size, na.rm = T,probs = c(0.025,0.975))
n_cell_size_b2=quantile(df$b_2.0.21d.cell.size, na.rm = T,probs = c(0.025,0.975))
n_cell_size_b3=quantile(df$b_3.0.21d.cell.size, na.rm = T,probs = c(0.025,0.975))
n_body_long_b1=quantile(df$b_1.long.starvation..ii., na.rm = T,probs = c(0.025,0.975))
n_body_long_b2=quantile(df$b_2.long.starvation..ii., na.rm = T,probs = c(0.025,0.975))
n_body_long_b3=quantile(df$b_3.long.starvation..ii., na.rm = T,probs = c(0.025,0.975))
n_AL_b1=quantile(df$b_1.AL, na.rm = T,probs = c(0.025,0.975))
n_AL_b2=quantile(df$b_2.AL, na.rm = T,probs = c(0.025,0.975))
n_AL_b3=quantile(df$b_3.AL, na.rm = T,probs = c(0.025,0.975))
n_AL_b4=quantile(df$b_4.AL, na.rm = T,probs = c(0.025,0.975))
n_RES_b1=quantile(df$b_1.RES, na.rm = T,probs = c(0.025,0.975))
n_RES_b2=quantile(df$b_2.RES, na.rm = T,probs = c(0.025,0.975))
n_RES_b3=quantile(df$b_3.RES, na.rm = T,probs = c(0.025,0.975))
n_RES_b4=quantile(df$b_4.RES, na.rm = T,probs = c(0.025,0.975))

a_apo_b1=quantile(df_aiptasia$b_1.apo, na.rm = T,probs = c(0.025,0.975))
a_apo_b2=quantile(df_aiptasia$b_2.apo, na.rm = T,probs = c(0.025,0.975))
a_symb_b1=quantile(df_aiptasia$b_1.symb, na.rm = T,probs = c(0.025,0.975))
a_symb_b2=quantile(df_aiptasia$b_2.symb, na.rm = T,probs = c(0.025,0.975))

#organise values as a dataframe
df_CI_boots=data.frame(n_body_b1, n_body_b2, n_cell_numb_b1, n_cell_numb_b2,
                       n_cell_size_b1, n_cell_size_b2, n_cell_size_b3, n_body_long_b1,
                       n_body_long_b2, n_body_long_b3, n_AL_b1, n_AL_b2, n_AL_b3,
                       n_AL_b4, n_RES_b1, n_RES_b2, n_RES_b3, n_RES_b4, a_apo_b1, a_apo_b2,
                       a_symb_b1, a_symb_b2)
#save as .csv file
write.csv(df_CI_boots, "Data/CI_boots.csv")
#note that the CI's for AL and RES for b1,b2,b3 and b4 are replaced in the Python code for 
#plotting and in the Suplemental Information as we select the outputs from the custom-written 
#simulated annealing method for optimising the log likelihood function in C++ 

#Compute the coefficient of determination R^2 of each multiphase model:

#log body size starvation 0-21d, 2 phases
#define the model from the parameter estimates obtained by maximising the log likelihood function
model_mp=function(x){
  if (x<=2){y=0.232*x+1.21}
  else {y=-0.024*(x-2)+1.67}
  return (y)
}

x_obs=df.2$Day
predictions=unlist(lapply(x_obs, model_mp))
SSE=sum((df.2$log.body.size-predictions)**2)
SST=sum((df.2$log.body.size-mean(df.2$log.body.size))**2)
R2=1-SSE/SST
print(R2)

#log cell number 0-21d, 2 phases
model_mp=function(x){
  if (x<=2){y=0.4*x+13.94}
  else {y=-0.1*(x-2)+14.74}
  return (y)
}
x_obs=df.2$Day
predictions=unlist(lapply(x_obs, model_mp))
SSE=sum((df.2$log.cell.number-predictions)**2)
SST=sum((df.2$log.cell.number-mean(df.2$log.cell.number))**2)
R2=1-SSE/SST
print(R2)

#log cell size 0-21d, 3-phases
model_mp=function(x){
  if (x<=2){y=0.032*x+10.9}
  else if (x>1 & x<=17){y=-0.0487*(x-2)+10.95}
  else {y=0.055*(x-17)+10.4}
  return (y)
}

x_obs=df.2$Day
predictions=unlist(lapply(x_obs, model_mp))
SSE=sum((df.2$log.cell.size-predictions)**2)
SST=sum((df.2$log.cell.size-mean(df.2$log.cell.size))**2)
R2=1-SSE/SST
print(R2)
plot(x_obs, predictions)

#long starvation (ii), 3-phases
model_mp=function(x){
  if (x<=5){y=0.0488*x+0.0147}
  else if (x>5 & x<=65){y=-0.0279*(x-5)+0.258}
  else {y=0.00213*(x-65)-1.137}
  return (y)
}

x_obs=df.1$Day
predictions=unlist(lapply(x_obs, model_mp))
SSE=sum((df.1$log.body.size-predictions)**2)
SST=sum((df.1$log.body.size-mean(df.1$log.body.size))**2)
R2=1-SSE/SST
print(R2)


#AL 4-phases
model_mp=function(x){
  if (x<=10){y=0.26*x}
  else if (x>10 & x<=100){y=-0.0219*(x-10)+3}
  else if (x>100 & x<=120){y=0.0537*(x-100)+1.4}
  else {y=-0.0227*(x-120)+3}
  return (y)
}
x_obs=df.AL$day
predictions=unlist(lapply(x_obs, model_mp))
SSE=sum((df.AL$log.body.size-predictions)**2)
SST=sum((df.AL$log.body.size-mean(df.AL$log.body.size))**2)
R2=1-SSE/SST
print(R2)

#RES 4-phases
model_mp=function(x){
  if (x<=10){y=0.168*x+0.8}
  else if (x>10 & x<=100){y=-0.0399*(x-10)+2.4}
  else if (x>100 & x<=120){y=0.0498*(x-100)+0.1}
  else {y=-0.0216*(x-120)+1.8}
  return (y)
}

x_obs=df.RES$day
predictions=unlist(lapply(x_obs, model_mp))
SSE=sum((df.RES$log.body.size-predictions)**2)
SST=sum((df.RES$log.body.size-mean(df.RES$log.body.size))**2)
R2=1-SSE/SST
print(R2)

