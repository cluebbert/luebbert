% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/memesave.R
\name{memesave}
\alias{memesave}
\title{GGsave for me! ggsave with some defaults}
\usage{
memesave(
  myfilename,
  device = "png",
  units = "in",
  height = 10,
  width = 15,
  bg = "white",
  ...
)
}
\arguments{
\item{myfilename}{filename to output}

\item{device}{device (can change this)}

\item{units}{defaults to inches}

\item{height}{some sensible and large default}

\item{width}{see height}

\item{bg}{no transparent background pngs}

\item{...}{extra arguments passed to ggplot2::ggsave()}
}
\value{
write plot
}
\description{
GGsave for me! ggsave with some defaults
}
