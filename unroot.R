# code unroots tree
library(ape)
library(ggplot2)
trees <- list.files("tree", full.names = T)

for (tree in trees){
  name <- basename(tree)
  tr <- read.tree(tree)
  if(tr$Nnode>1){
    unrooted_tr <- unroot(tr)
    write.tree(unrooted_tr, paste("unrooted/",name, sep=""))
 }
}
