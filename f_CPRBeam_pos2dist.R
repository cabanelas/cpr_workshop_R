## -----------------------------------------------------------------------------
## f_CPRBeam_pos2dist()
## R port of f_CPRBeam_pos2dist.m (Flora Sun, Univ. of Toronto, 2004)
## by: Alexandra Cabanelas -- NOV 2025
##
## Distance (in km) between two points on Earth given lat/lon in decimal degrees.
##   method 1 = planar approximation (good for points within ~tens of km)
##   method 2 = spherical geodesic   (farther apart; ignores Earth flattening)
##
## Vectorized: any argument may be a scalar or an equal-length vector, so you can
## pass one grid node (lat1, lon1 scalars) against all samples (lat2, lon2 vectors)
## in a single call.
## -----------------------------------------------------------------------------

f_CPRBeam_pos2dist <- function(lat1, lon1, lat2, lon2, method = 1) {

  # Normalise longitudes to 0..360, as in the original
  lon1 <- ifelse(lon1 < 0, lon1 + 360, lon1)
  lon2 <- ifelse(lon2 < 0, lon2 + 360, lon2)

  if (method == 1) {
    # --- Planar approximation ---
    km_per_deg_la <- 111.3237
    km_per_deg_lo <- 111.1350
    km_la <- km_per_deg_la * (lat1 - lat2)
    # always take the shorter arc in longitude
    dlon  <- abs(lon1 - lon2)
    dlon  <- ifelse(dlon > 180, dlon - 180, dlon)
    km_lo <- km_per_deg_lo * dlon * cos((lat1 + lat2) * pi / 360)
    dist  <- sqrt(km_la^2 + km_lo^2)

  } else {
    # --- Spherical geodesic ---
    R_aver <- 6374
    d2r    <- pi / 180
    a <- cos(lat1 * d2r) * cos(lat2 * d2r) * cos((lon1 - lon2) * d2r) +
         sin(lat1 * d2r) * sin(lat2 * d2r)
    a <- pmin(pmax(a, -1), 1)   # clamp to [-1, 1] so acos() never returns NaN
    dist <- R_aver * acos(a)    # (added over the MATLAB, which can NaN when node == sample)
  }

  dist
}
