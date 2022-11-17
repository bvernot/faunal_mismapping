library(ggplot2)
library(stringr)
library(RColorBrewer)
library(cowplot)

setwd("/mnt/expressions/benjamin_vernot/faunal_mismapping/data/once-again/6_parsed-bams/0_ref")

summary_df <- data.frame()
directories <- list.files()

index <- 0
for (directory in directories) {
  index <- index + 1
  inner_directories <- list.files(directory)

  input_df <- data.frame()

  for (inner_file in inner_directories) {
    combined_file <- paste(directory, "/", inner_file, sep="")

    input <- read.delim(combined_file)
    input$taxa <- str_split(inner_file, "-", simplify = T)[1]
    input_df <- rbind(input_df, input)

    summary_df <- rbind(summary_df, data.frame(taxa=str_split(inner_file, "-", simplify = T)[1],
                                                       read_count=nrow(input)/1E7,
                                                       prop_MQ25=sum(input$score>=25)/1E7,
                                                       mean_NM=mean(input$mismatch[input$score>=25]),
                                                       ref_db=gsub("._", "", directory)))
  }
}

summary_df$taxa <- factor(summary_df$taxa, levels = unique(summary_df$taxa[order(summary_df$read_count, decreasing = T)]))
input_df$taxa <- factor(input_df$taxa, levels =  levels(summary_df$taxa))

ggplot(summary_df)+
  geom_col(aes(taxa, read_count), fill="#1E90FF", col="black",alpha=0.3)+
  geom_col(aes(taxa, prop_MQ25),fill="#1E90FF", col="black")+
  geom_text(aes(taxa, read_count, label=round(read_count, digits = 5)),nudge_y = 0.01)+
  scale_y_continuous(expand = c(0,0), limits = c(0,1.05))+
  scale_fill_brewer(palette = "Paired")+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.margin = unit(c(0,0.25,0.25,0.25), "cm"),
        panel.background = element_rect(fill = "white"))+
  xlab("Taxa")+
  ylab("Proportion of mapping reads") # +
#   facet_wrap(~ref_db, ncol=1)

summary_df

# NM_plot <- ggplot(input_df)+
#   geom_boxplot(aes(taxa, mismatch))+
#   geom_point(aes(taxa, mean_NM), fill="#1E90FF")+
#   scale_y_continuous(expand = c(0,0), limits = c(0,10), breaks = c(0,2,4,6,8,10))+
#   scale_fill_brewer(palette = "Paired")+
#   theme(panel.grid = element_blank(),
#         panel.border = element_rect(color = "black", fill=NA),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         plot.margin = unit(c(0.25,0.25,0,0.25), "cm"), 
#         panel.background = element_rect(fill = "white"))+
#   xlab("") +
#   ylab("Mean number of missmatches")

# plot_grid(NM_plot, bar_plot, ncol = 1, rel_heights = c(1,2), align = "v", axis = "rl")