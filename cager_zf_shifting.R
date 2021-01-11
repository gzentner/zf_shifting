library(ZebrafishDevelopmentalCAGE)
library(CAGEr)

# Load data
data(ZebrafishSamples)
data(ZebrafishCAGE)

# Create character vector of sample names for simultaneous import
zf_samples <- as.character(ZebrafishSamples$sample)

# Import samples 
zf_cageset <- importPublicData(source = "ZebrafishDevelopment", 
                               dataset = "ZebrafishCAGE", 
                               group = "development", 
                               sample = zf_samples)

corr.m <- plotCorrelation2(zf_cageset, samples = c("zf_unfertilized_egg", "zf_fertilized_egg"),
                           tagCountThreshold = 1, applyThresholdBoth = FALSE, method = "pearson")

librarySizes(zf_cageset)

# Normalize
plotReverseCumulatives(zf_cageset, fitInRange = c(5,1000), onePlot = TRUE)

normalizeTagCount(zf_cageset, method = "powerLaw", fitInRange = c(5,1000), alpha = 1.21, T = 1e6)

# Cluster TSSs into TSRs
clusterCTSS(object = zf_cageset, threshold = 1, thresholdIsTpm = TRUE, nrPassThreshold = 1,
            method = "distclu", maxDist = 25, removeSingletons = FALSE, keepSingletonsAbove = 5)

