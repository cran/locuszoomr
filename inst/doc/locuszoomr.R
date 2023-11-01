## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE
)
library(locuszoomr)

## ----eval = FALSE-------------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#  BiocManager::install("ensembldb")
#  BiocManager::install("EnsDb.Hsapiens.v75")

## ----eval = FALSE-------------------------------------------------------------
#  install.packages("locuszoomr")

## ----eval = FALSE-------------------------------------------------------------
#  devtools::install_github("myles-lewis/locuszoomr")

## -----------------------------------------------------------------------------
library(locuszoomr)
data(SLE_gwas_sub)  ## limited subset of data from SLE GWAS
head(SLE_gwas_sub)

## ----warning=FALSE, message=FALSE, fig.dim=c(8, 7)----------------------------
library(EnsDb.Hsapiens.v75)
loc <- locus(SLE_gwas_sub, gene = 'UBE2L3', flank = 1e5,
             ens_db = "EnsDb.Hsapiens.v75")
summary(loc)
locus_plot(loc)

## ----eval=FALSE---------------------------------------------------------------
#  library(AnnotationHub)
#  ah <- AnnotationHub()
#  query(ah, c("EnsDb", "Homo sapiens"))

## ----eval=FALSE---------------------------------------------------------------
#  ## AnnotationHub with 25 records
#  ## # snapshotDate(): 2023-04-25
#  ## # $dataprovider: Ensembl
#  ## # $species: Homo sapiens
#  ## # $rdataclass: EnsDb
#  ## # additional mcols(): taxonomyid, genome, description, coordinate_1_based, maintainer,
#  ## #   rdatadateadded, preparerclass, tags, rdatapath, sourceurl, sourcetype
#  ## # retrieve records with, e.g., 'object[["AH53211"]]'
#  ##
#  ##              title
#  ##   AH53211  | Ensembl 87 EnsDb for Homo Sapiens
#  ##   ...        ...
#  ##   AH100643 | Ensembl 106 EnsDb for Homo sapiens
#  ##   AH104864 | Ensembl 107 EnsDb for Homo sapiens
#  ##   AH109336 | Ensembl 108 EnsDb for Homo sapiens
#  ##   AH109606 | Ensembl 109 EnsDb for Homo sapiens
#  ##   AH113665 | Ensembl 110 EnsDb for Homo sapiens

## ----eval=FALSE---------------------------------------------------------------
#  ensDb_v106 <- ah[["AH100643"]]
#  
#  # built-in mini dataset
#  data("SLE_gwas_sub")
#  loc <- locus(SLE_gwas_sub, gene = 'UBE2L3', fix_window = 1e6,
#               ens_db = ensDb_v106)
#  locus_plot(loc)

## ----eval = FALSE-------------------------------------------------------------
#  # Locus plot using SLE GWAS data from Bentham et al 2015
#  # FTP download full summary statistics from
#  # https://www.ebi.ac.uk/gwas/studies/GCST003156
#  library(data.table)
#  SLE_gwas <- fread('../bentham_2015_26502338_sle_efo0002690_1_gwas.sumstats.tsv')
#  loc <- locus(SLE_gwas, gene = 'UBE2L3', flank = 1e5,
#               ens_db = "EnsDb.Hsapiens.v75")
#  loc <- link_LD(loc, LDtoken = "your_token")
#  locus_plot(loc)

## ----fig.dim=c(8, 7)----------------------------------------------------------
loc <- locus(SLE_gwas_sub, gene = 'UBE2L3', flank = 1e5, LD = "r2",
             ens_db = "EnsDb.Hsapiens.v75")
locus_plot(loc, labels = c("index", "rs5754467"))

## ----eval = FALSE-------------------------------------------------------------
#  # Filter by gene biotype
#  locus_plot(loc, filter_gene_biotype = "protein_coding")
#  
#  # Custom selection of genes using gene names
#  locus_plot(loc, filter_gene_name = c('UBE2L3', 'RIMBP3C', 'YDJC', 'PPIL2',
#                                       'PI4KAP2', 'MIR301B'))

## ----fig.dim=c(7, 3.5)--------------------------------------------------------
genetracks(loc)

## ----fig.dim=c(7, 2.5)--------------------------------------------------------
# Limit the number of tracks
# Filter by gene biotype
# Customise colours
genetracks(loc, maxrows = 3, filter_gene_biotype = 'protein_coding',
           gene_col = 'grey', exon_col = 'orange', exon_border = 'darkgrey')

## ----eval=FALSE---------------------------------------------------------------
#  # obtain GTEx eQTL data from LDlinkR
#  # needs personal access token for LDlink
#  loc2 <- locus(SLE_gwas_sub, gene = 'IRF5', flank = c(7e4, 2e5), LD = "r2",
#                ens_db = "EnsDb.Hsapiens.v75")
#  loc2 <- link_eqtl(loc2, LDtoken = "..")
#  head(loc2$LDexp)  # show eQTL data

## ----eval=FALSE---------------------------------------------------------------
#  table(loc2$LDexp$Gene_Symbol)
#  ##      AHCYL2           IRF5            KCP        METTL2B           NRF1   RP11-212P7.2
#  ##           5           1618              1              2              9              3
#  ## ...
#  
#  table(loc2$LDexp$Tissue)
#  ##             Adipose - Subcutaneous          Adipose - Visceral (Omentum)
#  ##                                 78                                    69
#  ##                      Adrenal Gland                        Artery - Aorta
#  ##                                 38                                   109
#  ## ...

## ----eval=FALSE---------------------------------------------------------------
#  # just plot eQTL results
#  eqtl_plot(loc2, tissue = "Whole Blood", eqtl_gene = "IRF5")
#  
#  # overlay eqtl plot
#  overlay_plot(loc2, eqtl_gene = "IRF5")

## ----out.width='60%', fig.align="center", echo=FALSE--------------------------
knitr::include_graphics("eqtl_overlay.png")

## ----fig.dim=c(8, 7)----------------------------------------------------------
locb <- locus(SLE_gwas_sub, gene = 'UBE2L3', flank = 1e5, yvar = "beta",
              ens_db = "EnsDb.Hsapiens.v75")
locus_plot(locb)

## ----eval = FALSE-------------------------------------------------------------
#  genes <- c("STAT4", "IRF5", "UBE2L3")
#  
#  # generate list of 'locus' class objects, one for each gene
#  loclist <- lapply(genes, locus,
#                    data = SLE_gwas_sub,
#                    ens_db = "EnsDb.Hsapiens.v75",
#                    LD = "r2")
#  
#  ## produce 3 locus plots, one on each page
#  pdf("myplot.pdf")
#  multi_layout(loclist)
#  dev.off()
#  
#  ## place 3 locus plots in a row on a single page
#  pdf("myplot.pdf")
#  multi_layout(loclist, ncol = 3, labels = "index")
#  dev.off()
#  
#  ## full control
#  loc2 <- locus(SLE_gwas_sub, gene = 'IRF5', flank = c(7e4, 2e5), LD = "r2",
#                ens_db = "EnsDb.Hsapiens.v75")
#  loc3 <- locus(SLE_gwas_sub, gene = 'STAT4', flank = 1e5, LD = "r2",
#                ens_db = "EnsDb.Hsapiens.v75")
#  
#  pdf("myplot.pdf", width = 9, height = 6)
#  multi_layout(ncol = 3,
#               plots = {
#                 locus_plot(loc, use_layout = FALSE, legend_pos = 'topleft')
#                 locus_plot(loc2, use_layout = FALSE, legend_pos = NULL)
#                 locus_plot(loc3, use_layout = FALSE, legend_pos = NULL,
#                            labels = "index")
#               })
#  dev.off()

## ----echo=FALSE, message=FALSE, fig.dim = c(9,5)------------------------------
loc2 <- locus(SLE_gwas_sub, gene = 'IRF5', flank = c(7e4, 2e5), LD = "r2",
              ens_db = "EnsDb.Hsapiens.v75")
loc3 <- locus(SLE_gwas_sub, gene = 'STAT4', flank = 1e5, LD = "r2",
              ens_db = "EnsDb.Hsapiens.v75")

multi_layout(ncol = 3,
             plots = {
               locus_plot(loc, use_layout = FALSE, legend_pos = 'topleft')
               locus_plot(loc2, use_layout = FALSE, legend_pos = NULL)
               locus_plot(loc3, use_layout = FALSE, legend_pos = NULL,
                          labels = "index")
             })

## ----eval=FALSE---------------------------------------------------------------
#  pdf("myplot2.pdf", width = 6, height = 8)
#  # set up layered plot with 2 plots & a gene track; store old par() settings
#  oldpar <- set_layers(2)
#  scatter_plot(loc, xticks = FALSE)
#  line_plot(loc, col = "orange", xticks = FALSE)
#  genetracks(loc)
#  par(oldpar)  # revert par() settings
#  dev.off()

## ----echo=FALSE, fig.dim = c(6,8)---------------------------------------------
oldpar <- set_layers(2)
scatter_plot(loc, xticks = FALSE)
line_plot(loc, col = "orange", xticks = FALSE)
genetracks(loc)
par(oldpar)

## ----eval=FALSE---------------------------------------------------------------
#  dat <- SLE_gwas_sub
#  dat$p2 <- -log10(dat$p * 0.1)
#  locp <- locus(dat, gene = 'UBE2L3', flank = 1e5)
#  locp2 <- locus(dat, gene = 'UBE2L3', flank = 1e5, yvar = "p2")
#  
#  # set up overlaid plot with 1 plot & a gene track; store old par() settings
#  oldpar <- set_layers(1)
#  scatter_plot(locp, xticks = FALSE, pcutoff = NULL, ylim = c(0, 16))
#  scatter_plot(locp2, xticks = FALSE, pcutoff = NULL, chromCol = "orange",
#               pch = 22, add = TRUE)
#  genetracks(loc)
#  par(oldpar)  # revert par() settings

## ----echo=FALSE, message=FALSE, fig.dim = c(7,7)------------------------------
dat <- SLE_gwas_sub
dat$p2 <- -log10(dat$p * 0.1)
locp <- locus(dat, gene = 'UBE2L3', flank = 1e5, ens_db = "EnsDb.Hsapiens.v75")
locp2 <- locus(dat, gene = 'UBE2L3', flank = 1e5, yvar = "p2",
               ens_db = "EnsDb.Hsapiens.v75")

# set up overlaid plot with 1 plot & a gene track; store old par() settings
oldpar <- set_layers(1)
scatter_plot(locp, xticks = FALSE, pcutoff = NULL, ylim = c(0, 16))
scatter_plot(locp2, xticks = FALSE, pcutoff = NULL, chromCol = "orange",
             pch = 22, add = TRUE)
genetracks(loc)
par(oldpar)  # revert par() settings

## ----message=FALSE, fig.dim = c(7,7)------------------------------------------
# add vertical lines for gene of interest under the main plot
pf <- quote({
  v <- locp$TX[locp$TX$gene_name == "UBE2L3", c("start", "end")]
  abline(v = v, col = "orange")
})

pl <- quote({
  # add custom text label for index SNP
  lx <- locp$data$pos[locp$data$rsid == locp$index_snp]
  ly <- locp$data$logP[locp$data$rsid == locp$index_snp]
  text(lx, ly, locp$index_snp, pos = 4, cex = 0.8)
  # add extra points
  px <- rep(22.05e6, 3)
  py <- 10:12
  points(px, py, pch = 21, bg = "green")
  # add custom legend
  legend("topleft", legend = c("group A", "group B"),
         pch = 21, pt.bg = c("blue", "green"), bty = "n")
})

locus_plot(locp, pcutoff = NULL, panel.first = pf, panel.last = pl)

## ----eval=FALSE---------------------------------------------------------------
#  locus_ggplot(loc)

## ----eval=FALSE---------------------------------------------------------------
#  grid::grid.newpage()
#  gg_genetracks(loc)

## ----eval=FALSE---------------------------------------------------------------
#  p <- gg_scatter(loc)
#  p

## ----eval=FALSE---------------------------------------------------------------
#  gg_addgenes(p, loc)

## ----fig.dim = c(12, 7)-------------------------------------------------------
library(cowplot)
p1 <- locus_ggplot(loc, draw = FALSE)
p2 <- locus_ggplot(loc2, legend_pos = NULL, draw = FALSE)
plot_grid(p1, p2, ncol = 2)

## ----eval=FALSE---------------------------------------------------------------
#  library(ggpubr)
#  pdf("my_ggplot.pdf", width = 10)
#  ggarrange(p1, p2, ncol = 2)
#  dev.off()
#  
#  library(gridExtra)
#  pdf("my_ggplot.pdf", width = 10)
#  grid.arrange(p1, p2, ncol = 2)
#  dev.off()

## ----eval=FALSE---------------------------------------------------------------
#  locus_plotly(loc2)

## ----out.width='70%', fig.align="center", echo=FALSE--------------------------
knitr::include_graphics("plotly.png")

