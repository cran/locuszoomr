% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multi_layout.R
\name{multi_layout}
\alias{multi_layout}
\title{Layout multiple locus plots}
\usage{
multi_layout(
  plots,
  nrow = 1,
  ncol = 1,
  heights = c(3, 2),
  legend_pos = "topleft",
  ...
)
}
\arguments{
\item{plots}{Either an 'expression' to be evaluated which is a series of
calls to \code{\link[=locus_plot]{locus_plot()}} or similar plotting functions, or a list of 'locus'
class objects which are plotted in sequence.}

\item{nrow}{Number of rows of plots}

\item{ncol}{Number of columns of plots}

\item{heights}{Vector of length 2 specifying height for plot and gene tracks}

\item{legend_pos}{A keyword either "topleft" or "topright" or \code{NULL} to hide
the legend. Not invoked if \code{plots} is an expression. The legend is only
shown on one plot on each page.}

\item{...}{Optional arguments passed to \code{\link[=locus_plot]{locus_plot()}} if \code{plots} contains a
list}
}
\value{
No return value.
}
\description{
Produces pages with multiple locus plots on.
}
\examples{
if(require(EnsDb.Hsapiens.v75)) {

data(SLE_gwas_sub)
genes <- c("STAT4", "UBE2L3", "IRF5")
loclist <- lapply(genes, locus,
                  data = SLE_gwas_sub,
                  ens_db = "EnsDb.Hsapiens.v75",
                  LD = "r2")
## produce 3 locus plots, one on each page
multi_layout(loclist)

## place 3 locus plots in a row on a single page
multi_layout(loclist, ncol = 3)

## full control
loc <- locus(SLE_gwas_sub, gene = 'STAT4', flank = 1e5, LD = "r2",
             ens_db = "EnsDb.Hsapiens.v75")
loc2 <- locus(SLE_gwas_sub, gene = 'IRF5', flank = c(7e4, 2e5), LD = "r2",
              ens_db = "EnsDb.Hsapiens.v75")
loc3 <- locus(SLE_gwas_sub, gene = 'UBE2L3', LD = "r2",
              ens_db = "EnsDb.Hsapiens.v75")
multi_layout(ncol = 3,
             plots = {
               locus_plot(loc, use_layout = FALSE, legend_pos = 'topleft')
               locus_plot(loc2, use_layout = FALSE, legend_pos = NULL)
               locus_plot(loc3, use_layout = FALSE, legend_pos = NULL)
             })

}
}
\seealso{
\code{\link[=locus_plot]{locus_plot()}}
}
