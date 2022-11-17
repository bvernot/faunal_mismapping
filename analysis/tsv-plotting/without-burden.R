rm(list = ls())

library(data.table)
library(ggfortify)
library(ggplot2)
library(dplyr)
library(httpgd)
library(RColorBrewer)

file <- "/mnt/expressions/benjamin_vernot/faunal_mismapping/analysis/derived_allele_stats/derived_allele_stats.reads.tsv.gz" # nolint
file2 <- "/mnt/expressions/benjamin_vernot/faunal_mismapping/data/modified_refs/hominin_derived_sites/modified_sites.third.txt"

tsv_dt <- fread(file)
modified_dt <- fread(file2, col.names = c("chrom","pos","None","a1","a2","ref_fasta","ref_replace"))

colnames(tsv_dt)[15] <- "ref_type"

tsv_dt[b_3 >= 15, b_3 := 15]
tsv_dt <- na.omit(tsv_dt)

tsv_dt.merge <- merge(tsv_dt, modified_dt[,c(1,2,7)], by = c("chrom","pos"), sort = F)

tsv_dt[, .N, b_3]       # print amount in the b_3 group

summary_dt <- tsv_dt.merge[, .(b_neq_a1a2.sum = sum(base != a1 & base != a2),
                         b_neq_a1a2.total = .N,
                         base_eq_ref.sum = sum(base == ref),
                         base_eq_a4.sum = sum(base == ref_replace),
                         base_eq_a3.sum = sum(base != a1 & base != a2 & base != ref_replace)),
                         by = .(species, ref_type, b_3)]

# summary_dt[, b_neq_a1a2 := b_neq_a1a2.sum / b_neq_a1a2.total]
summary_dt[, base_eq_derived := base_eq_ref.sum / b_neq_a1a2.total]                                
summary_dt[, base_eq_a4 := base_eq_a4.sum / b_neq_a1a2.total]
summary_dt[, base_eq_a3 := base_eq_a3.sum / b_neq_a1a2.total]

spc <- "1b_tarsier"
summary_dt_melt <- melt(summary_dt, id.vars = c("species",
                                                "ref_type",
                                                "b_3",
                                                "b_neq_a1a2.sum",
                                                "b_neq_a1a2.total",
                                                "base_eq_ref.sum",
                                                "base_eq_a4.sum",
                                                "base_eq_a3.sum"))


ggplot(summary_dt_melt[species == spc], aes(b_3, value, col = variable)) +
    geom_line() +
    facet_wrap(~ref_type) +
    theme(text = element_text(size = 26))


second_melt <- melt(summary_dt, id.vars = c("species",
                                            "ref_type",
                                            "b_neq_a1a2.sum",
                                            "b_neq_a1a2.total",
                                            "base_eq_ref.sum",
                                            "base_eq_a4.sum",
                                            "base_eq_a3.sum"))


ggplot(summary_dt_melt[species == spc]) +
    geom_bar(aes(x = variable, y = value)) # +
#    facet_wrap(~ref_type) +
#    theme(text = element_text(size = 26))


summary_dt_melt
