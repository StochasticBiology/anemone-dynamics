library(ggplot2)
library(svglite)

# filenames and labels for model outputs -- extensible to any set of output files from the optimiser
input.set = c("cycle-AL.txt", "cycle-RES.txt")
label.set = c("Cyclic ad libitum", "Cyclic restricted")

# initial summary dataframe and plots lists
selected.stats = data.frame()
traces.plot = stats.plot = list()

# loop through experiments
for (input.index in 1:length(input.set)) {
  # labelling
  this.input = input.set[input.index]
  this.label = label.set[input.index]
  summary.fname = paste(c("summary-", this.input, ".csv"), collapse = "")
  traces.fname = paste(c("output-", this.input, ".csv"), collapse = "")
  
  # file I/O
  data.df = read.table(this.input, header = FALSE)
  colnames(data.df) = c("t", "y")
  stats.df = read.csv(summary.fname)
  traces.df = read.csv(traces.fname)
  
  # work through the statistics dataframe
  summary = data.frame()
  converged = data.frame()
  # look at each model in turn
  for (this.model in unique(stats.df$model)) {
    # the 9th column is "b1", the next are the later bs. there are only as many bs of interest as the model label
    min.param = 9
    max.param = 9 + this.model - 1
    # for each b (slope) in this model, pull the set of corresponding estimates from the bootstrap set
    for (i in min.param:max.param) {
      vals = stats.df[stats.df$model == this.model, i]
      # calculate s.d. range for each
      lo = mean(vals) - sd(vals)
      hi = mean(vals) + sd(vals)
      summary = rbind(summary,
                      data.frame(model = this.model, slope = i-8, lo = lo, hi = hi))
    }
    # for use later: identify log-likelihood outliers that correspond to non-converged numerics
    lls = stats.df$log.lik[stats.df$model==this.model]
    converged.boots = stats.df$boot.label[stats.df$model==this.model & stats.df$log.lik > mean(lls)-2*sd(lls)]
    converged = rbind(converged, traces.df[traces.df$model==this.model & traces.df$boot.label %in% converged.boots,])
  }
  
  # produce stats summary plot for all models for this experiment
  stats.plot[[input.index]] = ggplot(summary) +
    geom_rect(aes(xmin = lo, xmax = hi, ymin = slope - 0.25,ymax = slope + 0.25)) +
    facet_wrap( ~ model) + xlab("slope") + ylab("phase")
  
  # append the selected statistics (estimates of the four slopes for the four-phase model) to the selected dataframe
  for (i in 1:4) {
    selected.stats = rbind(selected.stats,
      data.frame(index = (input.index - 1) * 4 + i,
                 label = paste(c(this.label, i), collapse = " "),
                 lo = summary$lo[summary$model == 4 & summary$slope == i],
                 hi = summary$hi[summary$model == 4 & summary$slope == i]))
  }
  
  # time series plots -- overlay the bootstrap set of model fits (within 2s.d.s of the mean log lik) on the datapoints
  traces.plot[[input.index]] = ggplot() +
    geom_point(data = data.df, aes(x = t, y = y), size = 0.2, alpha = 0.5) +
    geom_line(data = converged, aes(x = t, y = y, color = factor(boot.label))) +
    facet_wrap( ~ model) +
    theme_classic() + theme(legend.position = "none")

}


par(mfrow=c(4,2))
stats.df = read.csv("summary-cycle-AL.txt.csv")
sub = stats.df[stats.df$model==4,]
sub = sub[sub$log.lik > mean(sub$log.lik)-1.5*sd(sub$log.lik),]
hist(sub$b1); hist(sub$b2); hist(sub$b3); hist(sub$b4)

stats.df = read.csv("summary-cycle-RES.txt.csv")
sub = stats.df[stats.df$model==4,]
sub = sub[sub$log.lik > mean(sub$log.lik)-1.5*sd(sub$log.lik),]
hist(sub$b1); hist(sub$b2); hist(sub$b3); hist(sub$b4)

# plot all four-phase model stats together
selected.plot = ggplot(selected.stats) +
  geom_rect(aes(xmin = lo, xmax = hi, ymin = index - 0.25, ymax = index + 0.25, fill = factor(label))) +
  xlab("Slope (+- s.d.)") + ylab("Phase") +
  scale_y_continuous(breaks = 1:nrow(selected.stats), labels = selected.stats$label) +
  theme_classic() + theme(legend.position = "none")

# output whatever specific plots we want to SVG

svglite("selected-stats-300.svg", width=4, height=4)
selected.plot
dev.off()

svglite("AL-traces-300.svg", width=4, height=4)
traces.plot[[1]]
dev.off()

svglite("RES-traces-300.svg", width=4, height=4)
traces.plot[[2]]
dev.off()
