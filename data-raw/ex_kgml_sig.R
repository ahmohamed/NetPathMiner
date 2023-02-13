library(NetPathMiner)
ex_kgml_sig <- KGML2igraph(c('data-raw/hsa04014.xml', 'data-raw/hsa04062.xml'), parse.as = 'signaling')
usethis::use_data(ex_kgml_sig, overwrite = TRUE)
