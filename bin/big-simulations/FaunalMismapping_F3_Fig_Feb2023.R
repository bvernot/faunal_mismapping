## This script creates the two heatmap images for the Faunal Mismapping 
## Faunalzation, the running of f3 in all combinations, and the creation of the F3 tables used as input for this script are outlined at: https://github.com/bvernot/faunal_mismapping/blob/main/bin/big-simulations/k_faunalize-f3.sh

rm(list=ls())

## Install pacakged
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}


## Set working directoty
work=paste("/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/13_F3Stats//")

## Read in colour files for data labels
colour_file=paste("/home/niall_cooke/Documents/SpeedRun_Notes/Heatmaps/Faunalizer_Labels_Colours.csv",sep="")
colour_data<-read.csv(colour_file,
                      col.names=c("Population.", "NotLenientLabel","NotLenientLabel_60","LenientLabel","Colour"),
                      header=FALSE)

## Set the palette for the heatmap (low, high)
my_palette <- colorRampPalette(c("red", "yellow"))(n = 30)

## The analysis was done in several different stages - the name of the one that is to be included in the paper is "4_SNP-subset-60-per-cent-modern-europe"
analysis_name=paste("4_SNP-subset-60-per-cent-modern-europe")



###### 1) Control image (i.e. pre-faunalization)

## Read in the data
control_file=paste(work,analysis_name,"/control/",analysis_name,"-control-HEATMAP_TABLE.csv",sep="")
control_data<-read.csv(control_file, sep = "", header = TRUE)

## Replace f3 vales of zero with "NA"
control_data[control_data == 0] <- NA

## Merge control data to the colour file
control_data_with_colour = merge(control_data, colour_data, by="Population.")

## Isolate the values for the heatmap
mat_data<-data.matrix(control_data[,2:(ncol(control_data))])

## Add the names of the individuals to the the matrix of values
rnames<-control_data_with_colour[,1]
rownames(mat_data) <- rnames

## Isolate the correspndingcolour values in the right order
colour_ordered<-as.vector(control_data_with_colour$Colour)


## Plot as a pdf (or a png)
#png(paste(work,analysis_name,"/control/",analysis_name,"-control-HEATMAP.png", sep=""), height=15, width=15, res=300, units="cm")

pdf(paste(work,analysis_name,"/control/",analysis_name,"-control-HEATMAP.pdf",sep=""))
heatmap.2(mat_data,
          trace="none",         # turns off trace lines inside the hat map
          scale="none",
          col=my_palette,
          cexRow=0.6,
          cexCol=0.6,
          srtRow = 0,
          srtCol = 90,
          keysize=1.25,
          key.title = NULL,
          key.xlab = "f3",
          key.ylab = NULL,
          density.info="none",
          Rowv = TRUE,
          Colv = TRUE,
          colRow = colour_ordered,
          colCol = colour_ordered,
          labRow = paste(as.vector(control_data_with_colour$LenientLabel),sep=""),
          labCol = paste(as.vector(control_data_with_colour$LenientLabel),sep=""),
          #Colv=if(symm)"Rowv" else TRUE,
          dendrogram="col",
          #notecol="green",
          margins =c(7,7))
dev.off()

###### 2) "Not lenient" image (i.e. post-faunalization)

## Read in the data
not_lenient_file=paste(work,analysis_name,"/not_lenient/",analysis_name,"-not_lenient-HEATMAP_TABLE.csv",sep="")
not_lenient_data<-read.csv(not_lenient_file, sep = "", header = TRUE)

## Replace f3 vales of zero with "NA"
not_lenient_data[not_lenient_data == 0] <- NA

## Merge control data to the colour file
not_lenient_data_with_colour = merge(not_lenient_data, colour_data, by="Population.")

## Isolate the values for the heatmap
mat_data<-data.matrix(not_lenient_data[,2:(ncol(not_lenient_data))])

## Add the names of the individuals to the the matrix of values
rnames<-not_lenient_data_with_colour[,1]
rownames(mat_data) <- rnames

## Isolate the correspnding colour values in the right order
colour_ordered<-as.vector(not_lenient_data_with_colour$Colour)

## Plot as a pdf (or a png)
#png(paste(work,analysis_name,"/not_lenient/",analysis_name,"-not_lenient-HEATMAP.png", sep=""), height=15, width=15, res=300, units="cm")
pdf(paste(work,analysis_name,"/not_lenient/",analysis_name,"-not_lenient-HEATMAP.pdf",sep=""))


#pdf("/home/niall_cooke/Documents/SpeedRun_Notes/test2.pdf")
heatmap.2(mat_data,
          trace="none",         # turns off trace lines inside the hat map
          scale="none",
          col=my_palette,
          cexRow=0.6,
          cexCol=0.6,
          srtRow = 0,
          srtCol = 90,
          keysize=1.25,
          key.title = NULL,
          key.xlab = "f3",
          key.ylab = NULL,
          density.info="none",
          Rowv = TRUE,
          Colv = TRUE,
          colRow = colour_ordered,
          colCol = colour_ordered,
          labRow = paste(as.vector(not_lenient_data_with_colour$NotLenientLabel_60),sep=""),
          labCol = paste(as.vector(not_lenient_data_with_colour$NotLenientLabel_60),sep=""),
          #Colv=if(symm)"Rowv" else TRUE,
          dendrogram="col",
          #notecol="green",
          margins =c(7,7))
dev.off()

