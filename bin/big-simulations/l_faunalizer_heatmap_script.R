## Install ggplots
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

rm(list=ls())

work=paste("/mnt/expressions/benjamin_vernot/faunal_mismapping/data/big-simulations/13_F3Stats//")
colour_file=paste(work,"Faunalizer_Labels_Colours.csv",sep="")

colour_data<-read.csv(colour_file,
                      col.names=c("Population.", "Label","Group","Colour"),
                      header=FALSE)

my_palette <- colorRampPalette(c("yellow", "red"))(n = 30)

############################## 1_SNP-subset-broad_analysis ############################## 
analysis_name=paste("1_SNP-subset-broad_analysis")

######################## Control ####################################
control_file=paste(work,analysis_name,"/control/",analysis_name,"-control-HEATMAP_TABLE.csv",sep="")

control_data<-read.csv(control_file, sep = "", header = TRUE)
control_data[control_data == 0] <- NA


control_data_with_colour = merge(control_data, colour_data, by="Population.")

rnames<-control_data_with_colour[,1]
mat_data<-data.matrix(control_data[,2:(ncol(control_data))])
rownames(mat_data) <- rnames
colour_ordered<-as.vector(control_data_with_colour$Colour)

png(paste(work,analysis_name,"/control/",analysis_name,"-control-HEATMAP.png", sep=""), height=15, width=15, res=300, units="cm")

heatmap.2(mat_data,
          trace="none",         # turns off trace lines inside the hat map
          scale="none",
          col=my_palette,
          cexRow=0.6,
          cexCol=0.6,
          srtRow = 0,
          srtCol = 90,
          keysize=1.25,
          #key.xlab = "Value",
          #key.ylab = NULL,
          #density.info="none",
          Rowv = TRUE,
          Colv = TRUE,
          colRow = colour_ordered,
          colCol = colour_ordered,
          labRow = paste(as.vector(control_data_with_colour$Label),sep=""),
          labCol = paste(as.vector(control_data_with_colour$Label),sep=""),
          #Colv=if(symm)"Rowv" else TRUE,
          #dendrogram="col",
          #notecol="green",
          margins =c(7,7))
dev.off()

######################## Not Lenient ######################## 

not_lenient_file=paste(work,analysis_name,"/not_lenient/",analysis_name,"-not_lenient-HEATMAP_TABLE.csv",sep="")
not_lenient_data<-read.csv(not_lenient_file, sep = "", header = TRUE)

not_lenient_data[not_lenient_data == 0] <- NA


not_lenient_data_with_colour = merge(not_lenient_data, colour_data, by="Population.")

rnames<-not_lenient_data_with_colour[,1]
mat_data<-data.matrix(not_lenient_data[,2:(ncol(not_lenient_data))])
rownames(mat_data) <- rnames
colour_ordered<-as.vector(not_lenient_data_with_colour$Colour)

png(paste(work,analysis_name,"/control/",analysis_name,"-not_lenient-HEATMAP.png", sep=""), height=15, width=15, res=300, units="cm")

heatmap.2(mat_data,
          trace="none",         # turns off trace lines inside the hat map
          scale="none",
          col=my_palette,
          cexRow=0.6,
          cexCol=0.6,
          srtRow = 0,
          srtCol = 90,
          keysize=1.25,
          #key.xlab = "Value",
          #key.ylab = NULL,
          #density.info="none",
          Rowv = TRUE,
          Colv = TRUE,
          colRow = colour_ordered,
          colCol = colour_ordered,
          labRow = paste(as.vector(not_lenient_data_with_colour$Label),sep=""),
          labCol = paste(as.vector(not_lenient_data_with_colour$Label),sep=""),
          #Colv=if(symm)"Rowv" else TRUE,
          #dendrogram="col",
          #notecol="green",
          margins =c(7,7))
dev.off()


######################## Lenient ######################## 
lenient_file=paste(work,analysis_name,"/lenient/",analysis_name,"-lenient-HEATMAP_TABLE.csv",sep="")
lenient_data<-read.csv(lenient_file, sep = "", header = TRUE)

lenient_data[lenient_data == 0] <- NA

lenient_data_with_colour = merge(lenient_data, colour_data, by="Population.")

rnames<-lenient_data_with_colour[,1]
mat_data<-data.matrix(lenient_data[,2:(ncol(lenient_data))])
rownames(mat_data) <- rnames
colour_ordered<-as.vector(lenient_data_with_colour$Colour)

png(paste(work,analysis_name,"/control/",analysis_name,"-lenient-HEATMAP.png", sep=""), height=15, width=15, res=300, units="cm")

heatmap.2(mat_data,
          trace="none",         # turns off trace lines inside the hat map
          scale="none",
          col=my_palette,
          cexRow=0.6,
          cexCol=0.6,
          srtRow = 0,
          srtCol = 90,
          keysize=1.25,
          #key.xlab = "Value",
          #key.ylab = NULL,
          #density.info="none",
          Rowv = TRUE,
          Colv = TRUE,
          colRow = colour_ordered,
          colCol = colour_ordered,
          labRow = paste(as.vector(lenient_data_with_colour$Label),sep=""),
          labCol = paste(as.vector(lenient_data_with_colour$Label),sep=""),
          #Colv=if(symm)"Rowv" else TRUE,
          #dendrogram="col",
          #notecol="green",
          margins =c(7,7))
dev.off()




############################## 2_SNP-subset-modern-europe ############################## 
analysis_name=paste("2_SNP-subset-modern-europe")

######################## Control ####################################
control_file=paste(work,analysis_name,"/control/",analysis_name,"-control-HEATMAP_TABLE.csv",sep="")

control_data<-read.csv(control_file, sep = "", header = TRUE)
control_data[control_data == 0] <- NA


control_data_with_colour = merge(control_data, colour_data, by="Population.")

rnames<-control_data_with_colour[,1]
mat_data<-data.matrix(control_data[,2:(ncol(control_data))])
rownames(mat_data) <- rnames
colour_ordered<-as.vector(control_data_with_colour$Colour)

png(paste(work,analysis_name,"/control/",analysis_name,"-control-HEATMAP.png", sep=""), height=15, width=15, res=300, units="cm")

heatmap.2(mat_data,
          trace="none",         # turns off trace lines inside the hat map
          scale="none",
          col=my_palette,
          cexRow=0.6,
          cexCol=0.6,
          srtRow = 0,
          srtCol = 90,
          keysize=1.25,
          #key.xlab = "Value",
          #key.ylab = NULL,
          #density.info="none",
          Rowv = TRUE,
          Colv = TRUE,
          colRow = colour_ordered,
          colCol = colour_ordered,
          labRow = paste(as.vector(control_data_with_colour$Label),sep=""),
          labCol = paste(as.vector(control_data_with_colour$Label),sep=""),
          #Colv=if(symm)"Rowv" else TRUE,
          #dendrogram="col",
          #notecol="green",
          margins =c(7,7))
dev.off()


######################## Not Lenient######################## 

not_lenient_file=paste(work,analysis_name,"/not_lenient/",analysis_name,"-not_lenient-HEATMAP_TABLE.csv",sep="")
not_lenient_data<-read.csv(not_lenient_file, sep = "", header = TRUE)

not_lenient_data[not_lenient_data == 0] <- NA


not_lenient_data_with_colour = merge(not_lenient_data, colour_data, by="Population.")

rnames<-not_lenient_data_with_colour[,1]
mat_data<-data.matrix(not_lenient_data[,2:(ncol(not_lenient_data))])
rownames(mat_data) <- rnames
colour_ordered<-as.vector(not_lenient_data_with_colour$Colour)

png(paste(work,analysis_name,"/control/",analysis_name,"-not_lenient-HEATMAP.png", sep=""), height=15, width=15, res=300, units="cm")

heatmap.2(mat_data,
          trace="none",         # turns off trace lines inside the hat map
          scale="none",
          col=my_palette,
          cexRow=0.6,
          cexCol=0.6,
          srtRow = 0,
          srtCol = 90,
          keysize=1.25,
          #key.xlab = "Value",
          #key.ylab = NULL,
          #density.info="none",
          Rowv = TRUE,
          Colv = TRUE,
          colRow = colour_ordered,
          colCol = colour_ordered,
          labRow = paste(as.vector(not_lenient_data_with_colour$Label),sep=""),
          labCol = paste(as.vector(not_lenient_data_with_colour$Label),sep=""),
          #Colv=if(symm)"Rowv" else TRUE,
          #dendrogram="col",
          #notecol="green",
          margins =c(7,7))
dev.off()


######################## Lenient ######################## 
lenient_file=paste(work,analysis_name,"/lenient/",analysis_name,"-lenient-HEATMAP_TABLE.csv",sep="")
lenient_data<-read.csv(lenient_file, sep = "", header = TRUE)

lenient_data[lenient_data == 0] <- NA

lenient_data_with_colour = merge(lenient_data, colour_data, by="Population.")

rnames<-lenient_data_with_colour[,1]
mat_data<-data.matrix(lenient_data[,2:(ncol(lenient_data))])
rownames(mat_data) <- rnames
colour_ordered<-as.vector(lenient_data_with_colour$Colour)

png(paste(work,analysis_name,"/control/",analysis_name,"-lenient-HEATMAP.png", sep=""), height=15, width=15, res=300, units="cm")

heatmap.2(mat_data,
          trace="none",         # turns off trace lines inside the hat map
          scale="none",
          col=my_palette,
          cexRow=0.6,
          cexCol=0.6,
          srtRow = 0,
          srtCol = 90,
          keysize=1.25,
          #key.xlab = "Value",
          #key.ylab = NULL,
          #density.info="none",
          Rowv = TRUE,
          Colv = TRUE,
          colRow = colour_ordered,
          colCol = colour_ordered,
          labRow = paste(as.vector(lenient_data_with_colour$Label),sep=""),
          labCol = paste(as.vector(lenient_data_with_colour$Label),sep=""),
          #Colv=if(symm)"Rowv" else TRUE,
          #dendrogram="col",
          #notecol="green",
          margins =c(7,7))
dev.off()

### Lenient
