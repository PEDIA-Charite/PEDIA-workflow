args <- commandArgs(trailingOnly = TRUE)

print(args)

library(foreign)
library(ggplot2)
library(reshape)
library(reshape2)
library(dplyr)
library(PerfMeas)
library(precrec);
library(readr)


colours <- c("#1b9e77","#d95f02","#7570b3","#e6ab02","#e7298a","#a6761d","#000000", "#B6B6B6")


#plot style
size.text = 25;
size.legend = 25;
size.line = 1;
size.geom_line = 2;
standard_style <- theme_bw() + theme(axis.text = element_text(colour = "black",size=size.text), axis.title = element_text(colour = "black",size=size.text), plot.title = element_text(size = size.text, face="bold"), panel.grid.major = element_blank() , panel.grid.minor = element_blank(), panel.border = element_blank(), axis.ticks = element_line(colour = "black", size=1), axis.line = element_line(colour = "black", size=size.line), legend.key =  element_blank(), legend.text = element_text(size=size.legend), legend.position="top", legend.direction="horizontal", legend.key.size = unit(2, 'lines'), legend.title=element_blank())+ theme(axis.line.x = element_line(color="black", size = size.line), axis.line.y = element_line(color="black", size = size.line))
angle_style <- standard_style + theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size = 15))


#create balanced accuracy
createBalancedAccuracy <- function(t){
  return((0.5*t$`True Positives`)/(t$`True Positives`+t$`False Negatives`) + (0.5*t$`True Negatives`)/(t$`True Negatives`+t$`False Positives`))
}

# arff to performance
arffToPerformance <- function(evaluation){
  positives <- max(evaluation$`True Positives`)
  total <- positives + max(evaluation$`False Positives`)

  result <- evaluation %>% select(one_of(c("True Positives","False Positives"))) %>% dplyr::rename(TP = `True Positives`, FP =  `False Positives`) %>% group_by(TP) %>% summarise_each(funs(min)) %>% mutate(value=(TP+FP)/total, variable=TP/positives) %>% filter(variable!=0)
  return(result)
}

# start loading data
data <- data.frame()
percentile <- data.frame()
orderROC <- c()
orderPRC <- c()
folder <- args[length(args)]
names <- c()
for(i in seq(1, length(args)-1, by = 3)) {


	#negLabel <- 'negatives'
	#predictions for aurocs/auprcs
	predictions <- read.table(args[i], sep="\t", quote="");

	colnames(predictions) <- c("label","prediction");

	#predictions <- predictions %>% mutate(label= ifelse(label == negLabel,0,1));

	sscurves <- evalmod(scores = predictions$prediction, labels = predictions$label);
	m<-attr(sscurves,"auc",exact=FALSE);
	AUROC <-  format(round(m[1,"aucs"],3), nsmall=3);
	AUPRC <-  format(round(m[2,"aucs"],3), nsmall=3);


	# read arff data for plots
	tmpData <- read.arff(args[i+1])

	#balanced accuracy
	tmpData$`Balanced Accuracy` <- createBalancedAccuracy(tmpData)
	tmpData$Threshold <-  reshape::rescaler.default(tmpData$Threshold, type = "range")

	#classes for plotting
	tmpData$CLASS <- args[i+2]
	tmpData$AUROC <- paste0(args[i+2]," (",AUROC,")")
	tmpData$AUPRC <- paste0(args[i+2]," (",AUPRC,")")

	#percentile curves
	positives <- max(tmpData$`True Positives`)
  total <- positives + max(tmpData$`False Positives`)

  percentile <-  bind_rows(percentile, tmpData %>% dplyr::rename(TP = `True Positives`, FP =  `False Positives`) %>% group_by(TP) %>% summarise_each(funs(min)) %>% mutate(value=(TP+FP)/total, variable=TP/positives))

	data <- bind_rows(data,tmpData)

	names <- append(names,args[i+2])
	orderROC <- append(orderROC,tmpData$AUROC[1])
	orderPRC <- append(orderPRC,tmpData$AUPRC[1])
}

name <- paste(names, collapse=".")
name <- gsub(" ", "-", name, fixed = TRUE)

# ROC
data$AUROC1 <- factor(data$AUROC, orderROC);

p <- ggplot(data, aes( x= `False Positive Rate`, y = `True Positive Rate`, group = AUROC, colour =  AUROC1))
p <- p + geom_line(size=size.geom_line) + labs(x= "False positive rate", y="True positive rate", colour = "Classifiers") + standard_style + guides(col = guide_legend(ncol = 2))
p <- p + scale_colour_manual(values=colours,breaks = orderROC)

ggsave(p, file=paste0(folder,"/","ROC.",name,".pdf"),width = 10.0,height = 10)

# PRC
data$AUPRC1 <- factor(data$AUPRC, orderPRC);
dataPRC <- data[data$Recall != 0,]

p <- ggplot(dataPRC, aes( x= dataPRC$Recall, y = dataPRC$Precision, group = AUPRC, colour = AUPRC1))
p <- p + geom_line(size=size.geom_line) + labs(x= "Recall", y="Precision", colour = "Classifiers")+ standard_style  + guides(col = guide_legend(ncol = 2))
p <- p + scale_colour_manual(values=colours,breaks = orderPRC)

ggsave(p, file=paste0(folder,"/","PRC.",name,".pdf"),width = 10.0,height = 10)


# Percentile
percentile$CLASS1 <- factor(percentile$CLASS, names);
p <- ggplot(percentile%>% filter(variable!=0),  aes(x=value, y=variable, colour=CLASS1))
p <- p + geom_line(size=size.geom_line)  + labs(y= "Sensitivity", x="Quantile of top ranked variants") + standard_style +  scale_x_continuous(breaks=c(10^seq(0,-8,by=-2)), trans="log10",labels=c(expression(10^0), expression(10^-2), expression(10^-4), expression(10^-6), expression(10^-8))) + guides(col = guide_legend(ncol = 2))
p <- p + scale_colour_manual(values=colours,breaks = names)
# scale_x_continuous(breaks=c(10^seq(0,-8,by=-2)), trans="log10",label=scientific_10)
ggsave(p, file=paste0(folder,"/","percentile.",name,".pdf"),width = 10.0,height = 10.0)




# Precision
data$CLASS1 <- factor(data$CLASS, names);
p <- ggplot(data, aes( x= Threshold, y = Precision, group = CLASS, colour = CLASS1))
p <- p + geom_line(size=size.geom_line) + labs(x= "Threshold", y="Precision", colour = "Classifiers")+ standard_style  + guides(col = guide_legend(ncol = 3))
p <- p + scale_colour_manual(values=colours,breaks = names)
ggsave(p, file=paste0(folder,"/","precision.",name,".pdf"),width = 10.0,height = 10)

# Recall

p <- ggplot(data, aes( x= Threshold, y = Recall, group = CLASS, colour = CLASS1))
p <- p + geom_line(size=size.geom_line) + labs(x= "Threshold", y="Recall", colour = "Classifiers")+ standard_style  + guides(col = guide_legend(ncol = 3))
p <- p + scale_colour_manual(values=colours,breaks = names)
ggsave(p, file=paste0(folder,"/","recall.",name,".pdf"),width = 10.0,height = 10)

# F-Score

p <- ggplot(data, aes( x= Threshold, y = FMeasure, group = CLASS, colour = CLASS1))
p <- p + geom_line(size=size.geom_line) + labs(x= "Threshold", y="F-score", colour = "Classifiers")+ standard_style  + guides(col = guide_legend(ncol = 3))
p <- p + scale_colour_manual(values=colours,breaks = names)
ggsave(p, file=paste0(folder,"/","f-score.",name,".pdf"),width = 10.0,height = 10)

# Balanced accuracy

p <- ggplot(data, aes( x= Threshold, y = `Balanced Accuracy`, group = CLASS, colour = CLASS1))
p <- p + geom_line(size=size.geom_line) + labs(x= "Threshold", y="Balanced Accuracy", colour = "Classifiers")+ standard_style  + guides(col = guide_legend(ncol = 3))
p <- p + scale_colour_manual(values=colours,breaks = names)
ggsave(p, file=paste0(folder,"/","balanced-accuracy.",name,".pdf"),width = 10.0,height = 10)
