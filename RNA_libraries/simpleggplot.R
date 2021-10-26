args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  #args[2] = "piRNAlibrary"
}


library(ggplot2)
library(stringr)
library(gridExtra)

tab=read.table(args[1],sep="\t",header=F)

test=data.frame(tab)

colnames(test)=c("Locus","PairBase","Position","Reads")

test$Reads=log10(test$Reads)

plot_list_199 <- list()

plot_list_1224 <- list()

for(locus in unique(as.character(test[,1]))){
temp=test[which(as.character(test[,1]) == locus),]

if(str_detect(locus,"21ur-199")){

coords=c(21,42,59,80,322,343,420,441,490,511,766,787)

p = ggplot(data = temp, aes(x= Position, y = Reads)) + geom_point() + facet_wrap(.~Locus) +
geom_vline(xintercept = coords, color = "red", linetype = "dotted") + ylab(expression(paste("Read count in ",log[10]))) + theme(plot.title = element_text(hjust = 0.5)) +
ylim(0,max(test$Reads)+.2)

plot_list_199 <- c(plot_list_199, list(p))

}


if(str_detect(locus,"21ur-1224")){

coords=c(151,172,332,353,501,522,824,845,1051,1072,1174,1195)

p = ggplot(data = temp, aes(x= Position, y = Reads)) + geom_point() + facet_wrap(.~Locus) +
geom_vline(xintercept = coords, color = "red", linetype = "dotted") + ylab(expression(paste("Read count in ",log[10]))) + theme(plot.title = element_text(hjust = 0.5)) +
ylim(0,max(test$Reads)+.2)

plot_list_1224 <- c(plot_list_1224, list(p))

}


}

pdf(paste(as.character(args[2]),"_199.pdf",sep=""), width= 12, height = 7)
do.call(grid.arrange, c(plot_list_199, list(nrow = 1) ))
dev.off()


pdf(paste(as.character(args[2]),"_1224.pdf",sep=""), width= 12, height = 7)
do.call(grid.arrange, c(plot_list_1224, list(nrow = 1) ))
dev.off()
