\name{BCEF_PTC}

\docType{data}

\alias{BCEF_PTC}

\title{Percent Tree Cover over Bonanza Creek Experimental Forest}

\description{
Percent tree cover (Hansen et al. 2013) over a subset of Bonanza Creek
Experimental Forest, AK. Can be used as the set of prediction locations
and covariate for models fit using the \code{BCEF} data.

The \code{BCEF_PTC} matrix columns are longitude (x), latitude (y), and Landsat derived percent tree cover
(PTC) for 237,617 observations.  Longitude and latitude are
in Albers Equal Area (proj4string
"+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs"). 
}

\usage{
data(BCEF_PTC)
}

\references{
Hansen, M.C., Potapov, P.V., Moore, R., Hancher, M., Turubanova,
S.A., Tyukavina, A.,Thau, D., Stehman, S.V., Goetz, S.J., Loveland,
T.R., Kommareddy, A., Egorov, A., Chini, L., Justice, C.O., and
Townshend, J.R.G. (2013), High-Resolution Global Mapsof 21st-Century
Forest Cover Change, Science, 342, 850-853.
}

\format{
  A matrix containing 237,617 rows and 3 columns named longitude,
  latitude, and PTC. 
}

\seealso{
  \code{\link{BCEF}}.
}

\keyword{datasets}
