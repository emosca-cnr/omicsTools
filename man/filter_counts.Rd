% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_counts.R
\name{filter_counts}
\alias{filter_counts}
\title{Filter and normalize count data}
\usage{
filter_counts(
  X = NULL,
  out_dir = NULL,
  col = NULL,
  width = 180,
  height = 180,
  res = 300,
  units = "mm",
  prior.count = 3,
  filter.by.cpm = TRUE,
  min.cpm = 3,
  min.counts = 10,
  subset.remove = NULL,
  min.samples = NULL
)
}
\arguments{
\item{X}{genes-by-samples count matrix}

\item{out_dir}{output directory}

\item{col}{vector of colors; this must match the number of columns of X}

\item{width}{image width}

\item{height}{image height}

\item{res}{resolution}

\item{units}{image units}

\item{prior.count}{prior.count for cpm() function}

\item{filter.by.cpm}{if TRUE filter by cmp rather then raw counts}

\item{min.cpm}{mininum count per million a gene must have in at least min_sample_perc percentage of samples}

\item{min.counts}{minimum number of counts}

\item{subset.remove}{these samples will be removed}

\item{min.samples}{minimum number of samples}
}
\value{
filtered count matrix
}
\description{
Filter and normalize count data
}
