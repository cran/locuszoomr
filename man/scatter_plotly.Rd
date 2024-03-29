% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scatter_plotly.R
\name{scatter_plotly}
\alias{scatter_plotly}
\title{Locus scatter plotly}
\usage{
scatter_plotly(
  loc,
  index_snp = loc$index_snp,
  pcutoff = 5e-08,
  chromCol = "royalblue",
  sigCol = "red",
  xlab = NULL,
  ylab = NULL,
  yzero = (loc$yvar == "logP"),
  showLD = TRUE,
  LD_scheme = c("grey", "royalblue", "cyan2", "green3", "orange", "red", "purple"),
  marker_outline = "black",
  marker_size = 7,
  recomb_col = "blue"
)
}
\arguments{
\item{loc}{Object of class 'locus' to use for plot. See \link{locus}.}

\item{index_snp}{Specifies index SNP to be shown in a different colour and
symbol. Defaults to the SNP with the lowest p-value. Set to \code{NULL} to not
show this.}

\item{pcutoff}{Cut-off for p value significance. Defaults to p = 5e-08. Set
to \code{NULL} to disable.}

\item{chromCol}{Colour for normal points if \code{LD} is \code{FALSE} when the locus
object is made.}

\item{sigCol}{Colour for significant points if \code{LD} is \code{FALSE}.}

\item{xlab}{x axis title.}

\item{ylab}{y axis title.}

\item{yzero}{Logical whether to force y axis limit to include y=0.}

\item{showLD}{Logical whether to show LD with colours}

\item{LD_scheme}{Vector of colours for plotting LD. The first colour is for
SNPs which lack LD information. The next 5 colours are for r^2 or D' LD
results ranging from 0 to 1 in intervals of 0.2. The final colour is for
the index SNP.}

\item{marker_outline}{Specifies colour for outlining points.}

\item{marker_size}{Value for size of markers in plotly units.}

\item{recomb_col}{Colour for recombination rate line if recombination rate
data is present. Set to \code{NA} to hide the line. See \code{\link[=link_recomb]{link_recomb()}} to add
recombination rate data.}
}
\value{
A \code{plotly} scatter plot.
}
\description{
Produces a scatter plot from a 'locus' class object using plotly.
}
\seealso{
\code{\link[=locus]{locus()}} \code{\link[=locus_plotly]{locus_plotly()}}
}
