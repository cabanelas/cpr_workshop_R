## -----------------------------------------------------------------------------
## f_CPRBeam_SolarAzEl()
## R port of f_CPRBeam_SolarAzEl.m
##   original MATLAB: Darin C. Koblick 2009/2013;
##    algorithm from http://stjarnhimlen.se/comp/tutorial.html#5
##
## Given a UTC time and a site location, returns the solar azimuth and elevation
## (altitude) angles at that site.
##
## Inputs
##   utc : POSIXct vector, assumed to be in UTC (build with tz = "UTC")
##   lat : site latitude  in degrees  (-90..90,   S(-) N(+))
##   lon : site longitude in degrees  (-180..180, W(-) E(+))
##   alt : site altitude above sea level in km (usually 0)
##   lat/lon/alt may be scalars or vectors the same length as utc.
##
## Output
##   data.frame with columns:
##     Az : solar azimuth   (degrees)
##     El : solar elevation (degrees above the horizon; negative = below)
##
## Verified against the original function's documented example:
##   f_CPRBeam_SolarAzEl(as.POSIXct("1991-05-19 13:00:00", tz="UTC"), 50, 10, 0)
##   -> Az ~ 223.60, El ~ 53.41
##
## Alternative: oce::sunAngle(utc, longitude = lon, latitude = lat) returns the
## same quantities ($azimuth, $altitude) from a maintained package. Results agree
## to a fraction of a degree; this port is kept so workshop output is identical.
## -----------------------------------------------------------------------------

f_CPRBeam_SolarAzEl <- function(utc, lat, lon, alt = 0) {

  d2r <- pi / 180

  # --- calendar components from the UTC timestamp ---
  lt  <- as.POSIXlt(utc, tz = "UTC")
  yr  <- lt$year + 1900
  mo  <- lt$mon + 1
  day <- lt$mday
  hr  <- lt$hour
  mn  <- lt$min
  sc  <- lt$sec

  # --- Julian date (custom formula from the original juliandate subfunction) ---
  yy <- yr
  mm <- mo
  idx <- mm <= 2
  yy[idx] <- yy[idx] - 1
  mm[idx] <- mm[idx] + 12
  jd <- floor(365.25 * (yy + 4716)) + floor(30.6001 * (mm + 1)) + 2 -
        floor(yy / 100) + floor(floor(yy / 100) / 4) +
        day - 1524.5 + (hr + mn / 60 + sc / 3600) / 24

  d <- jd - 2451543.5

  # --- Keplerian elements for the Sun (geocentric) ---
  w      <- 282.9404 + 4.70935e-5 * d          # longitude of perihelion
  e      <- 0.016709 - 1.151e-9  * d           # eccentricity
  M      <- (356.0470 + 0.9856002585 * d) %% 360  # mean anomaly
  L      <- w + M                               # Sun's mean longitude
  oblecl <- 23.4393 - 3.563e-7 * d              # obliquity of the ecliptic

  # --- eccentric anomaly (auxiliary angle) ---
  E <- M + (180 / pi) * e * sin(M * d2r) * (1 + e * cos(M * d2r))

  # --- rectangular coords in the ecliptic plane ---
  x <- cos(E * d2r) - e
  y <- sin(E * d2r) * sqrt(1 - e^2)

  r <- sqrt(x^2 + y^2)
  v <- atan2(y, x) * (180 / pi)      # true anomaly
  lon_sun <- v + w                   # longitude of the Sun

  xeclip <- r * cos(lon_sun * d2r)
  yeclip <- r * sin(lon_sun * d2r)
  zeclip <- 0.0

  # --- rotate to equatorial rectangular coords ---
  xequat <- xeclip
  yequat <- yeclip * cos(oblecl * d2r) + zeclip * sin(oblecl * d2r)
  # NOTE: the original hardcodes 23.4406 on this line (not oblecl); kept as-is
  # so results match the MATLAB exactly.
  zequat <- yeclip * sin(23.4406 * d2r) + zeclip * cos(oblecl * d2r)

  # --- RA and declination (with altitude correction rolled into r) ---
  r     <- sqrt(xequat^2 + yequat^2 + zequat^2) - (alt / 149598000)
  RA    <- atan2(yequat, xequat) * (180 / pi)
  delta <- asin(zequat / r) * (180 / pi)

  # --- local sidereal time and hour angle ---
  UTH     <- hr + mn / 60 + sc / 3600
  GMST0   <- ((L + 180) %% 360) / 15
  SIDTIME <- GMST0 + UTH + lon / 15
  HA      <- SIDTIME * 15 - RA

  # --- to horizontal (topocentric) coords ---
  x <- cos(HA * d2r) * cos(delta * d2r)
  y <- sin(HA * d2r) * cos(delta * d2r)
  z <- sin(delta * d2r)

  xhor <- x * cos((90 - lat) * d2r) - z * sin((90 - lat) * d2r)
  yhor <- y
  zhor <- x * sin((90 - lat) * d2r) + z * cos((90 - lat) * d2r)

  Az <- atan2(yhor, xhor) * (180 / pi) + 180
  El <- asin(pmin(pmax(zhor, -1), 1)) * (180 / pi)  # clamp guards asin domain

  data.frame(Az = Az, El = El)
}
