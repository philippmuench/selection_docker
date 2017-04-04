library(ggplot2)
library(stringr)

# import codeml data
dat <- read.csv2("normal_dnds.txt", header=T)
df <- as.data.frame(dat)
df$dnds <- as.numeric(as.matrix(df$dnds))
df$dN <- as.numeric(as.matrix(df$dN))
df$dS <- as.numeric(as.matrix(df$dS))

# process eden output
listeden <- list.files("dnds/", full.names = T)
df$eden_dnds <- -1
df$eden_dN <- -1
df$eden_dS <- -1
for (file in listeden){
  dat <- read.table(file, header=T)
  dN <- (sum(dat$Nd) + sum(dat$Nd_tip)) / sum(dat$N)
  dS <- (sum(dat$Sd) + sum(dat$Sd_tip)) / sum(dat$S)
  dN_dS <- dN/dS
  filename <- str_extract(file,  perl('(?<=[/])([^/]+)(?=\\.[^.]+)'))
  name <- sub('\\.txt.DnDsRatio$', '', filename)
  df[which(df$file == name),]$eden_dnds <- dN_dS
  df[which(df$file == name),]$eden_dN <- dN
  df[which(df$file == name),]$eden_dS <- dS
}

# remove if no eden value
df <- df[-which(df$eden_dN == -1),]

# plot dN
p <- ggplot(df, aes(x=dN, y=eden_dN))
p <- p + geom_point() + scale_x_log10() #+ scale_y_log10()
p <- p + theme_minimal()
p

p <- ggplot(df, aes(x=dS, y=eden_dS))
p <- p + geom_point() #+ scale_y_log10() #+ scale_y_log10()
p <- p + theme_minimal()
p
