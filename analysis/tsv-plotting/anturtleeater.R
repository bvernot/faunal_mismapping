library(data.table)
library(ggplot2)
library(dplyr)
library(httpgd)
library(RColorBrewer)

file <- "/mnt/expressions/benjamin_vernot/faunal_mismapping/analysis/derived_allele_stats/derived_allele_stats.reads.tsv.gz" # nolint
tsv_dt <- fread(file)
colnames(tsv_dt)[15] <- "ref_type"

tsv_dt[b_3 >= 15, b_3 := 15]
tsv_dt <- na.omit(tsv_dt)

tsv_dt[, .N, b_3]       # print amount in the b_3 group

summary_dt <- tsv_dt[, .(b_neq_a1a2.sum = sum(base != a1 & base != a2),
                         b_neq_a1a2.total = .N,
                         base_eq_ref.sum = sum(base == ref)),
                         by = .(species, ref_type, b_3)]

summary_dt[, b_neq_a1a2 := b_neq_a1a2.sum / b_neq_a1a2.total]

summary_dt <- summary_dt[, base_eq_ref := summary_dt[, base_eq_ref.sum]
                                        / summary_dt[, b_neq_a1a2.total]]                                    

spc <- "5_mouse"
ggplot(summary_dt, aes(x = b_3, y = b_neq_a1a2)) +
    geom_point() +
    geom_line(col = "red") +
    facet_grid(species~ref_type) +
    geom_line(aes(x = b_3, y = base_eq_ref), col = "#4a4a4a")


ggplot(summary_dt, aes(x = b_3, y = b_neq_a1a2.total, color = species)) +
    geom_point() +
    scale_y_continuous(trans='log10') +
    geom_line() +
    scale_color_brewer(palette = "Paired") +
    facet_wrap(~ref_type)#

spc <- "4_bison"
ggplot(summary_dt[, "ref_type"] == "1_third", aes(x = b_3, y = b_neq_a1a2, spec = species)) +
    geom_point() +
    geom_line() +
    geom_line(aes(x = b_3, y = base_eq_ref))

unique(tsv_dt[, "species"])


summary_dt[, ref_type, by = "1_third"]

summary_dt
