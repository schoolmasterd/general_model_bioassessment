#This script creates figure 2 of the main text

# get data
d1<-read.csv("Data/derive_model_queries_summary.csv",check.names = F)

# transpose
n <- d1$operation # store names
d1 <- as.data.frame(t(d1[,-c(1:4)])) # transpose excluding first column (names)
colnames(d1) <- n # recover names

# check name order
ifelse(names(table(rownames(d1)==c("indicator.species","ecotoxicology","oe","ambi","mmi",'dscore','Upsilon')))==T,'leaf names are correct.','incorrect leaf names. edit tlab')

# set names for display
tlab=c('IS','ET','O/E','AMBI','MMI','D-score','MVC')

# get distance matrix
dist(d1,method="binary")

# write the plot to a file
pdf('Output/figure_2.pdf', width=4, height=4.5)
plot(hclust(dist(d1,method="binary")),
     axes = T, frame.plot = F,
     main='',ylab='Jaccard Distance',sub='',
     labels=tlab,
     xlab='',)
dev.off()
