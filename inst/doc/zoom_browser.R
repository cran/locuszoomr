## ----eval=FALSE---------------------------------------------------------------
# # FTP download full summary statistics (308 MB) from SLE GWAS from
# # https://www.ebi.ac.uk/gwas/studies/GCST003156
# library(data.table)
# SLE_gwas <- fread('../bentham_2015_26502338_sle_efo0002690_1_gwas.sumstats.tsv')

## ----eval=FALSE---------------------------------------------------------------
# library(locuszoomr)
# library(EnsDb.Hsapiens.v75)
# 
# zoom(SLE_gwas, ens_db = "EnsDb.Hsapiens.v75")

## ----out.width='100%', echo=FALSE---------------------------------------------
knitr::include_graphics("zoom1.png")

## ----eval=FALSE---------------------------------------------------------------
# library(rtracklayer)
# recomb.hg19 <- import.bw("/../hapMapRelease24CombinedRecombMap.bw")
# 
# zoom(SLE_gwas, ens_db = "EnsDb.Hsapiens.v75", recomb = recomb.hg19)

## ----eval=FALSE---------------------------------------------------------------
# zoom(SLE_gwas, ens_db = "EnsDb.Hsapiens.v75", recomb = recomb.hg19, ld_token = "my_token")

## ----out.width='100%', echo=FALSE---------------------------------------------
knitr::include_graphics("zoom2.png")

## ----eval=FALSE---------------------------------------------------------------
# # load ensembl Homo sapiens database v110
# library(AnnotationHub)
# ah <- AnnotationHub()
# query(ah, c("EnsDb", "Homo sapiens"))  # shows latest versions
# ensDb_v110 <- ah[["AH113665"]]
# 
# # load hg38 recombination rate data downloaded from UCSC
# recomb1000G <- import.bw("/../recomb1000GAvg.bw")
# 
# # double GWAS, eQTL first
# zoom(eqtlres, data2 = RA_gwas, ens_db = ensDb_v110,
#      eqtl_gene = "gene", beta = c("beta", NA),
#      recomb = recomb1000G, ld_token = "my_token")

## ----out.width='100%', echo=FALSE---------------------------------------------
knitr::include_graphics("dbl_gwas_manh.png")
knitr::include_graphics("dbl_loc_eqtl.png")

