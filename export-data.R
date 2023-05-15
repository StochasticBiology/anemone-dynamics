library(ggplot2)

# read data
df = read.csv2("ALLDATA_bodysize.csv")

# summary plot
ggplot(df, aes(x=day,y=log(mm2))) + geom_point() + facet_wrap(~experiment)

# subset and write ad libitum feeding data
sub = df[df$experiment=="Cyclicfeeding_starving" & df$condition=="AL" & !is.na(df$mm2),]
write.table(data.frame(x=sub$day, y=log(sub$mm2)), "cycle-AL.txt", row.names=FALSE, col.names=FALSE)

# subset and write restricted feeding data
sub = df[df$experiment=="Cyclicfeeding_starving" & df$condition=="RES" & !is.na(df$mm2),]
write.table(data.frame(x=sub$day, y=log(sub$mm2)), "cycle-RES.txt", row.names=FALSE, col.names=FALSE)

# subset and write others
sub = df[df$experiment=="10d_feeding" & !is.na(df$mm2),]
write.table(data.frame(x=sub$day, y=log(sub$mm2)), "data-10d.txt", row.names=FALSE, col.names=FALSE)

sub = df[df$experiment=="140d_starved" & !is.na(df$mm2),]
write.table(data.frame(x=sub$day, y=log(sub$mm2)), "data-140d.txt", row.names=FALSE, col.names=FALSE)

sub = df[df$experiment=="170d_200d_starved_refed" & df$condition == "170d" & !is.na(df$mm2),]
write.table(data.frame(x=sub$day, y=log(sub$mm2)), "data-170d-170.txt", row.names=FALSE, col.names=FALSE)

sub = df[df$experiment=="170d_200d_starved_refed" & df$condition == "200d" & !is.na(df$mm2),]
write.table(data.frame(x=sub$day, y=log(sub$mm2)), "data-170d-200.txt", row.names=FALSE, col.names=FALSE)

sub = df[df$experiment=="200d_starved" & !is.na(df$mm2),]
write.table(data.frame(x=sub$day, y=log(sub$mm2)), "data-200d.txt", row.names=FALSE, col.names=FALSE)

sub = df[df$experiment=="21d_starving" & !is.na(df$mm2),]
write.table(data.frame(x=sub$day, y=log(sub$mm2)), "data-21d.txt", row.names=FALSE, col.names=FALSE)


