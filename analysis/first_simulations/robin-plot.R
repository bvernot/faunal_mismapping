library(ggplot2)
library(stringr)

setwd("/mnt/expressions/benjamin_vernot/faunal_mismapping/data/simulated-reads/parsed-bams")

files <- list.files()
input_df <- data.frame()
summary_df <- data.frame()

for (file in files) {
    taxa <- str_split(file, "_", simplify = TRUE)[1]
    input <- read.delim(file)
    input$taxa <- taxa
    input_df <- rbind(input_df, input)
    summary_df <- rbind(summary_df, data.frame(taxa = taxa,
                                               read_count = nrow(input)))
}
