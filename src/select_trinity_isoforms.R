#!/usr/bin/env Rscript

#######
# LOG #
#######

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

#############
# LIBRARIES #
#############

library(tidyverse)

###########
# GLOBALS #
###########

abundance_file <- snakemake@input[["abundance"]]

########
# MAIN #
########

isoform_list <- read_tsv(abundance_file)

isoform_by_expression <- isoform_list %>% 
    group_by(gene_id) %>% 
    slice(which.max(IsoPct)) %>% 
    ungroup() %>% 
    select(transcript_id) %>% 
    distinct(transcript_id)


isoform_by_length <- isoform_list %>% 
    group_by(gene_id) %>% 
    slice(which.max(length)) %>% 
    ungroup() %>% 
    select(transcript_id) %>% 
    distinct(transcript_id)

# write output
write_delim(isoform_by_length,
            snakemake@output[["length"]],
            col_names = FALSE)

write_delim(isoform_by_expression,
            snakemake@output[["expression"]],
            col_names = FALSE)

# write log
sessionInfo()
