df = read.table("longer-starved_aiptasia.txt", header=TRUE)
df.1 = df[df$symbiosis=="apo",]
write.table(data.frame(x=df.1$day, y=log(as.numeric(gsub(",",".",df.1$mm2)))), "data-apo.txt", row.names=FALSE, col.names=FALSE)
df.2 = df[df$symbiosis=="sym",]
write.table(data.frame(x=df.2$day, y=log(as.numeric(gsub(",",".",df.2$mm2)))), "data-sym.txt", row.names=FALSE, col.names=FALSE)

