\name{BCEF}

\docType{data}

\alias{BCEF}

\title{Forest Canopy Height from NASA Goddard's LiDAR Hyperspectral and
  Thermal (G-LiHT) over Bonanza Creek Experimental Forest}

\description{
Forest canopy height (FCH) estimates from NASA Goddard's LiDAR Hyperspectral and
Thermal (G-LiHT; Cook et al. 2013) Airborne Imager and percent tree
cover (Hansen et al. 2013) over a subset of
Bonanza Creek Experimental Forest, AK, collected in Summer 2014.

The \code{BCEF} matrix columns are longitude (x), latitude (y), forest
canopy height (FCH) in meters from ground, and Landsat derived percent tree cover
(PTC) for 188,717 observations.  Longitude and latitude are
in Albers Equal Area (proj4string
"+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs"). The last column (holdout) identifies a subset of data useful for assessing wall-to-wall predictive performance.
}

\usage{
data(BCEF)
}

\source{G-LiHT data were downloaded from \url{https://gliht.gsfc.nasa.gov}.}

\references{
Cook, B.D., L.W. Corp, R.F. Nelson, E.M. Middleton, D.C. Morton,
J.T. McCorkel, J.G. Masek, K.J. Ranson, and V. Ly. (2013) NASA Goddard's
Lidar, Hyperspectral and Thermal (G-LiHT) airborne imager. Remote
Sensing 5:4045-4066.

Hansen, M.C., Potapov, P.V., Moore, R., Hancher, M., Turubanova,
S.A., Tyukavina, A.,Thau, D., Stehman, S.V., Goetz, S.J., Loveland,
T.R., Kommareddy, A., Egorov, A., Chini, L., Justice, C.O., and
Townshend, J.R.G. (2013), High-Resolution Global Mapsof 21st-Century
Forest Cover Change, Science, 342, 850-853.
}

\format{
  A matrix containing 188,717 rows and 6 columns named longitude,
  latitude, FCH, PTC, and holdout. 
}

\keyword{datasets}
