# model = (number of phases)
# parameters theta = { b1, c, sigma, b2, tau1, b3, tau2, b4,tau3 }
# note that c is the intercept b0 in the text documents

#definition of the negative log likelihood function
neg.log.likelihood = function(theta, model, y, x)
{
  # parse argument to parameters
  b1 = theta[1]
  c = theta[2]
  sigma = abs(theta[3])
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

initial.guess = function(expt=0, model.type, bx, by) {
  # parameters theta = { b1, c, sigma, b2, tau1, b3, tau2, b4,tau3 }
  # initial guess -- flat, sigma = 1, changepoints at 1/4, 1/2 and 3/4 through x
  theta = c(0, mean(by[bx==min(bx)]), 1, 0, max(bx)/3, 0, 2*max(bx)/3, 0, 8*max(bx)/9)
  if(expt == 0) {
  theta = c(0, 2, 2, 0, 10, 0, 100, 0, 110)
  } else {  
    if(expt == 1 & model.type != 3) { theta = c(0.5, mean(by[bx==min(bx)]), 1, 0,75, 0, 80, 0,180) }
  if(expt == 1 & model.type == 3) { theta = c(0.5, mean(by[bx==min(bx)]), 1, 0, 5, 0, 80, 0,180) }
  if(expt == 2) { theta = c(0, mean(by[bx==min(bx)]), 1, 0, 3, 0, 11, 0, 18) }
  if(expt == 3) { theta = c(0, mean(by[bx==min(bx)]), 1, 0, 3, 0, 15, 0, 18) }
  if(expt == 4) { theta = c(0, mean(by[bx==min(bx)]), 1, 0, 3, 0, 15, 0, 18) }
  if(expt == 5) { theta = c(0, mean(by[bx==min(bx)]), 1, 0, 5, 0, 90, 0, 120) }
  if(expt == 6) { theta = c(0, mean(by[bx==min(bx)]), 1, 0, 10, 0, 70, 0, 120) }
  if(expt == 7) { theta = c(0, mean(by[bx==min(bx)]), 1, 0, 14, 0, 40, 0, 100) }
  if(expt == 8) { theta = c(0, mean(by[bx==min(bx)]), 1, 0, 14, 0, 40, 0, 100) }
  }
  return(theta)
}




