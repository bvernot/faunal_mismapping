rm(list = ls())

library(data.table)
library(ggfortify)
library(ggplot2)
library(dplyr)
library(httpgd)
library(RColorBrewer)

file <- "/mnt/expressions/benjamin_vernot/faunal_mismapping/analysis/derived_allele_stats/derived_allele_stats.reads.tsv.gz" # nolint
file2 <- "/mnt/expressions/benjamin_vernot/faunal_mismapping/data/modified_refs/hominin_derived_sites/modified_sites.third.txt"

file <- "/Volumes/expressions/benjamin_vernot/faunal_mismapping/analysis/derived_allele_stats/derived_allele_stats.reads.tsv.gz" # nolint
file2 <- "/Volumes/expressions/benjamin_vernot/faunal_mismapping/data/modified_refs/hominin_derived_sites/modified_sites.third.txt"


tsv_dt <- fread(file)
modified_dt <- fread(file2, col.names = c("chrom","pos","None","a1","a2","ref_fasta","ref_replace"))

colnames(tsv_dt)[15] <- "ref_type"

tsv_dt[b_3 >= 15, b_3 := 15]
tsv_dt <- na.omit(tsv_dt)

tsv_dt.merge <- merge(tsv_dt, modified_dt[,.(chrom, pos, ref_replace)], by = c("chrom","pos"), sort = F)

tsv_dt[, .N, b_3]       # print amount in the b_3 group

summary_dt <- tsv_dt.merge[, .(b_neq_a1a2.sum = sum(base != a1 & base != a2),
                         N = .N,
                         N.a1a2 = sum(base == a1 | base == a2),
                         ref = sum(base == ref),
                         anc = sum(base != ref & (base == a1 | base == a2)),
                         a3 = sum(base == ref_replace),
                         a4 = sum(base != a1 & base != a2 & base != ref_replace)),
                         by = .(species, ref_type, b_3)]

# summary_dt[, b_neq_a1a2 := b_neq_a1a2.sum / N]
summary_dt[, base_eq_derived := ref / N]                                
summary_dt[, base_eq_a4 := a4 / N]
summary_dt[, base_eq_a3 := a3 / N]


summary_dt_melt <- melt(summary_dt, 
                        id.vars = c("species",
                                    "ref_type",
                                    "b_3"), 
                        measure.vars = c(
                          'ref', 'a4', 'a3'
                        ))


ggplot(summary_dt_melt[species == spc], aes(b_3, value, col = variable)) +
    geom_line() +
    facet_wrap(~ref_type) +
    theme(text = element_text(size = 26))


summary_dt_melt.sum <- melt(summary_dt, 
                        id.vars = c("species",
                                    "ref_type",
                                    "b_3", 'N', 'N.a1a2'), 
                        measure.vars = c(
                          'ref', 'a4', 'a3', 'anc'
                        ))

summary_dt_melt.stats1 <-
  summary_dt_melt.sum[b_3 > 0,
                    .(value = sum(value) / sum(N),
                      value.a1a2 = sum(value) / sum(N.a1a2)), 
                    .(species, ref_type, variable)]

ggplot(summary_dt_melt.stats1,
       aes(x = species, y = value, fill = variable)) +
  geom_col() +
  facet_wrap(~ref_type) +
  theme(text = element_text(size = 26),
        axis.text.x = element_text(angle=-70))


ggplot(summary_dt_melt.stats1[species %like% 'dog' & variable %like% 'ref'],
       aes(x=ref_type, y=value)) + geom_col()


plot.spc <- c("1a_tamarin","1b_tarsier","2_dog","3_bear","4_bison","5_mouse")
ggplot(summary_dt_melt.stats1[species %in% plot.spc],
       aes(x = species, y = value, fill = variable)) +
  geom_col() +
  facet_wrap(~ref_type) +
  theme(text = element_text(size = 26),
        axis.text.x = element_text(angle=-70))

ggplot(summary_dt_melt.stats1[species %in% plot.spc & 
                                variable %in% c('ref', 'anc') &
                                ref_type != '2_N'],
       aes(x = species, y = value.a1a2, fill = variable)) +
  geom_col() +
  facet_wrap(~ref_type) +
  theme(text = element_text(size = 26),
        axis.text.x = element_text(angle=-70))



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
