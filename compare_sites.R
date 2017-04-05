source("plot.functions.R")



tigrs <- c("TIGR00016", "TIGR00064", "TIGR00070", "TIGR00101", "TIGR00117",
           "TIGR00152", "TIGR00170", "TIGR00173", "TIGR00183", "TIGR00185",
           "TIGR00186", "TIGR00196", "TIGR00218", "TIGR00219", "TIGR00385",
           "TIGR00390", "TIGR00392", "TIGR00408", "TIGR00409", "TIGR00422", 
           "TIGR00446", "TIGR00457", "TIGR00458", "TIGR00462", "TIGR00468",
           "TIGR00469", "TIGR00471", "TIGR00476", "TIGR00487", "TIGR00488",
           "TIGR00490", "TIGR00492", "TIGR00498", "TIGR00499", "TIGR00501",
           "TIGR00513", "TIGR00589", "TIGR00593", "TIGR00595", "TIGR00614",
           "TIGR00972", "TIGR01188", "TIGR01306", "TIGR01316", "TIGR01497"
)

for(i in 1:length(tigrs)){
tigr <- tigrs[i]
cat(tigr)


#p1 <- create_msa_plot(dnds_path ="dnds/TIGR00016.txt.DnDsRatio.txt",
#                gap_path="gap/TIGR00016.gap.txt", points=T)
#p1

#p2 <- create_msa_plot_slac(dnds_path ="X/TIGR00016.csv",
#                      gap_path="gap/TIGR00016.gap.txt", points=T)
#p2

p3 <- create_joined(dnds_path_eden =paste("dnds/",tigr,".txt.DnDsRatio.txt", sep=''),
  dnds_path_alt =paste("X/", tigr,".csv", sep=''),
  gap_path=paste("gap/", tigr,".gap.txt"), points=T)
pdf(file=paste(tigr,'.pdf', sep=''))
print(p3)
dev.off()
}