% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/moveVIZ.R
\name{moveVIZ}
\alias{moveVIZ}
\title{Creates a plot/table of moveSIM() results}
\usage{
moveVIZ(
  data,
  type = "plot",
  title = "moveSIM results",
  aspect_ratio = 1,
  xlim = NULL,
  ylim = NULL
)
}
\arguments{
\item{data}{Data to be plotted, this object should be the output from
moveSIM().}

\item{type}{"plot" or "summary_table", default "plot".}

\item{title}{Title for the plot that is output.}

\item{aspect_ratio}{Aspect ratio, defaults to 1.}

\item{xlim}{Optionally specify desired x limits as a numeric vector: c(low,hi)}

\item{ylim}{Optionally specify desired y limits as a numeric vector: c(low,hi)}
}
\value{
Plot or table displaying moveSIM() results.
}
\description{
When type="plot", function plots the movement tracks versus the the straight
line track between the origin and destination (unless the destination was
unspecified in the call to moveSIM(), then straight line track is omitted).
When type="summary_table", a summary table is output.
}
\examples{

# 1. Define Population and Run moveSIM()

pop1 <- as.species(x=-100, y=55)

EX2=moveSIM(replicates=2,days=5,env_rast=ex_raster, search_radius=550,
sigma=.1, dest_x=-108.6, dest_y=26.2, mot_x=.8, mot_y=.8,
modeled_species=pop1,optimum=.6, n_failures=5, fail_thresh=.40,
direction="R",write_results=FALSE,single_rast=TRUE,mortality = TRUE)

# 2. Run moveVIZ() on your result
moveVIZ(EX2,title="Visualizing MoveSIM results",type="plot")

moveVIZ(EX2, type="summary_table")

}
