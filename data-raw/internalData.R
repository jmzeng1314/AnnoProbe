## code to prepare `DATASET` dataset goes here
load("data-raw/exists_anno_list.rda")
load("data-raw/gpl_list.rda")
load("data-raw/series.accession.rda")
usethis::use_data(exists_anno_list, gpl_list, series.accession, internal = TRUE, overwrite = TRUE)
