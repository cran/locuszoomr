% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/locus_ggplot.R
\name{locus_ggplot}
\alias{locus_ggplot}
\title{Locus plot using ggplot2}
\usage{
locus_ggplot(
  loc,
  height = unit(5, "cm"),
  index_snp = loc$index_snp,
  pcutoff = 5e-08,
  scheme = c("royalblue", "red", "purple"),
  size = 2,
  cex.axis = 1,
  cex.lab = 1,
  xlab = NULL,
  ylab = NULL,
  xticks = "top",
  border = FALSE,
  showLD = TRUE,
  LD_scheme = c("grey", "royalblue", "cyan2", "green3", "orange", "red", "purple"),
  legend_pos = "topleft",
  draw = TRUE,
  ...
)
}
\arguments{
\item{loc}{Object of class 'locus' to use for plot. See \code{\link[=locus]{locus()}}.}

\item{height}{A unit object specifying height of the lower gene track.}

\item{index_snp}{Specifies index SNP to be shown in a different colour and
symbol. Defaults to the SNP with the lowest p-value. Set to \code{NULL} to not
show this.}

\item{pcutoff}{Cut-off for p value significance. Defaults to p = 5e-08. Set
to \code{NULL} to disable.}

\item{scheme}{Vector of 3 colors if LD is not shown: 1st = normal points, 2nd
= colour for significant points, 3rd = index SNP.}

\item{size}{Specifies size for points.}

\item{cex.axis}{Specifies font size for axis numbering.}

\item{cex.lab}{Specifies font size for axis titles.}

\item{xlab}{x axis title.}

\item{ylab}{y axis title.}

\item{xticks}{Logical whether x axis numbers and axis title are plotted.}

\item{border}{Logical whether a bounding box is plotted around top and bottom
plots.}

\item{showLD}{Logical whether to show LD with colours}

\item{LD_scheme}{Vector of colours for plotting LD. The first colour is for
SNPs which lack LD information. The next 5 colours are for r2 or D' LD
results ranging from 0 to 1 in intervals of 0.2. The final colour is for
the index SNP.}

\item{legend_pos}{Position of legend e.g. "topleft", "topright" or ggplot2
settings. Set to \code{NULL} to hide legend.}

\item{draw}{Logical whether to call \code{\link[=grid.draw]{grid.draw()}} to draw the plot.}

\item{...}{Additional arguments passed to \code{\link[=gg_genetracks]{gg_genetracks()}} to control
colours of gene tracks etc.}
}
\value{
Returns a ggplot2 plot containing a scatter plot with genetracks
underneath.
}
\description{
Genomic locus plot similar to locuszoom.
}
\examples{
if(require(EnsDb.Hsapiens.v75)) {
data(SLE_gwas_sub)
loc <- locus(SLE_gwas_sub, gene = 'IRF5', flank = c(7e4, 2e5), LD = "r2",
             ens_db = "EnsDb.Hsapiens.v75")
locus_ggplot(loc)
}
}
\seealso{
\code{\link[=gg_scatter]{gg_scatter()}} \code{\link[=gg_genetracks]{gg_genetracks()}}
}
