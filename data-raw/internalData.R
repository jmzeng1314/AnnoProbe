## code to prepare `DATASET` dataset goes here
load("data-raw/exists_anno_list.rda")
load("data-raw/gpl_list.rda")
load("data-raw/series.accession.rda")
load("data-raw/humanGTF.rda")
load("data-raw/mouseGTF.rda")
load("data-raw/ratGTF.rda")
usethis::use_data(exists_anno_list, gpl_list, series.accession,
                  humanGTF, mouseGTF, ratGTF,
                  internal = TRUE, overwrite = TRUE)
