\name{CHM}

\docType{data}

\alias{CHM}

\title{Canopy Height Model from NASA Goddard's LiDAR Hyperspectral and
  Thermal (G-LiHT)}

\description{
Canopy Height Model (CHM) from NASA Goddard's LiDAR Hyperspectral and
Thermal (G-LiHT; Cook et al. 2013) Airborne Imager over a subset of
Harvard Forest Simes Tract, MA, collected in Summer 2012.

The \code{CHM} matrix columns are longitude, latitude, and canopy height
(m) from ground for 17,23,137 observations. Longitude and latitude are
in UTM Zone 18 (proj4string "+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0").
}

\usage{
data(CHM)
}

\source{Data were downloaded from \url{https://gliht.gsfc.nasa.gov} with metadata available at \url{ftp://fusionftp.gsfc.nasa.gov/G-LiHT/Simes_Jun2012/metadata/Simes_Jun2012_metadata.pdf}.}

\references{
Cook, B.D., L.W. Corp, R.F. Nelson, E.M. Middleton, D.C. Morton,
J.T. McCorkel, J.G. Masek, K.J. Ranson, and V. Ly. (2013) NASA Goddard's
Lidar, Hyperspectral and Thermal (G-LiHT) airborne imager. Remote Sensing 5:4045-4066.
}

\format{
  A matrix containing 17,23,137 rows and 3 columns named longitude, latitude, and CHM. 
}

\keyword{datasets}
