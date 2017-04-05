library(ggplot2)
library(stringr)

names <- c("TIGR00016", "TIGR00064", "TIGR00070", "TIGR00101", "TIGR00117",
           "TIGR00152", "TIGR00170", "TIGR00173", "TIGR00183", "TIGR00185",
           "TIGR00186", "TIGR00196", "TIGR00218", "TIGR00219", "TIGR00385",
           "TIGR00390", "TIGR00392", "TIGR00408", "TIGR00409", "TIGR00422", 
           "TIGR00446", "TIGR00457", "TIGR00458", "TIGR00462", "TIGR00468",
           "TIGR00469", "TIGR00471", "TIGR00476", "TIGR00487", "TIGR00488",
           "TIGR00490", "TIGR00492", "TIGR00498", "TIGR00499", "TIGR00501",
           "TIGR00513", "TIGR00589", "TIGR00593", "TIGR00595", "TIGR00614",
           "TIGR00972", "TIGR01188", "TIGR01306", "TIGR01316", "TIGR01497"
           )


nums <- c(4, 4, 4, 4, 4,
          4, 5, 4, 4, 5,
          4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 
          6, 5, 4, 5, 4,
          4, 4, 4, 4, 4,
          4, 6, 4, 4, 4,
          4,4, 4, 5,  4,
          6, 6, 6,6,6
          )

hypy_dnds <- c(0.351, 0.324, 0.301, 0.410, 0.230,
               0.337, 0.226, 0.433, 0.357, 0.323,
               0.279, 0.375, 0.266, 0.447, 0.460,
               0.321052, 0.305032, 0.315611, 0.430501, 0.324813,
               0.374881, 0.321806, 0.345806, 0.26118,  0.218748,
               0.323695, 0.324827, 0.309928, 0.239135, 0.373966,
               0.355002, 0.2442,  0.360657, 0.194837, 0.270743,
               0.24374, 0.262539, 0.266614, 0.333748 , 0.370157,
               0.420503, 0.385379, 0.263919, 0.264632, 0.298434)



df <- data.frame(names=names, hypy = hypy_dnds, nums=nums)

# process eden output
listeden <- list.files("dnds", full.names = T)
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
  cat(name,'/n')
  if (length(which(df$names == name))>0){
    df[which(df$names == name),]$eden_dnds <- dN_dS
    df[which(df$names == name),]$eden_dN <- dN
    df[which(df$names == name),]$eden_dS <- dS
    df[which(df$names == name),]$eden_S <- S
    df[which(df$names == name),]$eden_N <- N
  }
}

df$gap <- -1
listgap <- list.files("gap", full.names = T)
for (file in listgap){
  dat <- read.table(file, header=F)
  filename <- str_extract(file,  perl('(?<=[/])([^/]+)(?=\\.[^.]+)'))
  name <- sub('\\.gap$', '', filename)
  cat(name,'/n')
  if (length(which(df$names == name))>0){
    df[which(df$names == name),]$gap <- mean(dat$V1)
  }
}


# remove if no eden value
#df <- df[-which(df$eden_dnds == -1),]

df$completness <- 1- df$gap
# plot dN
p <- ggplot(df, aes(x=hypy, y=eden_dnds, color=gap))
p <- p + geom_point() #+ scale_x_log10() 
p <- p + theme_minimal() + coord_equal() + scale_colour_gradient(low = "#4db898", high = "#ea5420")
p

cor.test(df$hypy, df$eden_dnds)
