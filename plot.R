library(ggplot2)
library(stringr)

# import codeml data
#dat <- read.csv2("normal_dnds.txt", header=T)
dat <- read.csv2("pair4.txt", header=T)

df <- as.data.frame(dat)
df$dN <- as.numeric(as.matrix(df$dN))
df$dS <- as.numeric(as.matrix(df$dS))
df$dnds <- as.numeric(as.matrix(df$dN))/ as.numeric(as.matrix(df$dS))

# process eden output
listeden <- list.files("dnds/", full.names = T)
df$eden_dnds <- -1
df$eden_dN <- -1
df$eden_dS <- -1
df$eden_N <- -1
df$eden_S <- -1
for (file in listeden){
  dat <- read.table(file, header=T)
  N <- sum(dat$N)
  S <- sum(dat$S)
  dN <- (sum(dat$Nd) + sum(dat$Nd_tip)) / sum(dat$N)
  dS <- (sum(dat$Sd) + sum(dat$Sd_tip)) / sum(dat$S)
  dN_dS <- dN/dS
  filename <- str_extract(file,  perl('(?<=[/])([^/]+)(?=\\.[^.]+)'))
  name <- sub('\\.txt.DnDsRatio$', '', filename)
  if (length(which(df$file == name))>0){
    df[which(df$file == name),]$eden_dnds <- dN_dS
    df[which(df$file == name),]$eden_dN <- dN
    df[which(df$file == name),]$eden_dS <- dS
    df[which(df$file == name),]$eden_S <- S
    df[which(df$file == name),]$eden_N <- N
  }
}

# remove if no eden value
df <- df[-which(df$eden_dN == -1),]

# plot dN
p <- ggplot(df, aes(x=dN, y=eden_N))
p <- p + geom_point() #+ scale_x_log10() 
p <- p + theme_minimal() 
p

p <- ggplot(df, aes(x=dN, y=eden_S))
p <- p + geom_point() #+ scale_y_log10() #+ scale_y_log10()
p <- p + theme_minimal()
p
