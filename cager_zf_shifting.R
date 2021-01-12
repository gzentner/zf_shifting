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

librarySizes(zf_cageset)

# Normalize
plotReverseCumulatives(zf_cageset, fitInRange = c(5,1000), onePlot = TRUE)

normalizeTagCount(zf_cageset, method = "powerLaw", fitInRange = c(5,1000), alpha = 1.21, T = 1e6)

# Cluster TSSs into TSRs
clusterCTSS(object = zf_cageset, threshold = 1, thresholdIsTpm = TRUE, nrPassThreshold = 1,
            method = "distclu", maxDist = 25, removeSingletons = FALSE, keepSingletonsAbove = 5,
            useMulticore = TRUE)

cumulativeCTSSdistribution(zf_cageset, clusters = "tagClusters", useMulticore = TRUE)
quantilePositions(zf_cageset, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)

aggregateTagClusters(zf_cageset, tpmThreshold = 5, maxDist = 100, qLow = 0.1, qUp = 0.9)

cumulativeCTSSdistribution(zf_cageset, clusters = "consensusClusters", useMulticore = TRUE)

scoreShift(zf_cageset, groupX = "zf_unfertilized_egg", groupY = "zf_prim20", testKS = TRUE, useTpmKS = FALSE)

shifting.promoters <- getShiftingPromoters(zf_cageset,
      tpmThreshold = 10,
      fdrThreshold = 0.01)
head(shifting.promoters)