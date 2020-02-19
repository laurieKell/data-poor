trjs=dlply(myers,.(assessid), with, TrajFromCoords(data.frame(x=bbmsy,y=sp)))

angular<-function(x,y){
   
    angle =as.numeric(coord2rad(x,y))
    angvec=angle[-1]-angle[-length(angle)]
    
    data.frame(angle=angle[-length(angle)],angvec=angvec)}

trj2=dlply(myers,.(assessid), with, angular(bbmsy,sp))
trj3=mlply(seq(length(trjs)),function(x) cbind(trjs[[x]],trj2[[x]]))

characteriseTrajectory <- function(trj) {
  # Measures of speed
  derivs <- TrajDerivatives(trj)
  mean_speed <- mean(derivs$speed)
  sd_speed <- sd(derivs$speed)
  
  # Measures of straightness
  sinuosity <- TrajSinuosity(trj)
  resampled <- TrajRediscretize(trj, .001)
  Emax <- TrajEmax(resampled)
  
  # Periodicity
  corr <- TrajDirectionAutocorrelations(resampled, 60)
  first_min <- TrajDAFindFirstMinimum(corr)
  
  angle = mean(trj$angvec,na.rm=T)
  
  # Return a list with all of the statistics for this trajectory
  list(mean_speed = mean_speed,
       sd_speed = sd_speed,
       sinuosity = sinuosity,
       Emax = Emax,
       angvec=angle,
       min_deltaS = first_min[1],
       min_C = first_min[2]
  )}

stats<-TrajsMergeStats(trj3, characteriseTrajectory)
stats$assessid=names(trjs)

customPcaPlot <- function(x, xlabs, xcols, choices = 1L:2L, ycol = "#ff2222aa", ...) {
  # Draw points
  pts <- t(t(x$x[, choices]))
  plot(pts, type = "p", 
       xlim = extendrange(pts[, 1L]), ylim = extendrange(pts[, 2L]), 
       asp = 1,
       xlab = "PC1", ylab = "PC2", pch = 16, col = xcols, ...)
  text(pts, labels = xlabs, pos = 1, ...)
  
  # Draw arrows
  axs <- t(t(x$rotation[, choices])) * 3.5
  text(axs, labels = dimnames(axs)[[1L]], col = ycol, ...)
  arrows(0, 0, axs[, 1L] * .8, axs[, 2L] * .8, length = .1, col = ycol)
}
# ---

# Here we are operating on the statistics from the previous example

# Get rid of NAs in stats because prcomp can't handle them.
# First fix min_deltaS, and add an extra column which flags non-periodic trajectories
pcaStats <- TrajsStatsReplaceNAs(stats, "min_deltaS", 
                                 replacementValue = 2 * max(stats$min_deltaS, na.rm = TRUE),
                                 flagColumn = "no_first_min")
# Also get rid of NAs in min_C - no need to add another column since it would duplicate no_first_min
pcaStats <- TrajsStatsReplaceNAs(pcaStats, "min_C",
                                 replacementValue = 2 * max(stats$min_C, na.rm = TRUE))
# Perform the PCA
PCA <- prcomp(pcaStats, scale. = TRUE)
# Plot it using custom plotting function. Could just call biplot instead
customPcaPlot(PCA, tracks$category, tracks$col, cex = .8)
legend("bottomleft", c("Spider", "Mimic", "Ant"), pch = 16, 
       col = c('red', 'blue', 'black'), inset = c(0.01, .02))





