##
## aim:
## calculate coherence between NGT and Arctic2k and between NGT and 20CR
## relation:
## NGT paper; https://github.com/EarthSystemDiagnostics/exNGT
##

path <- "C:/Users/mhoerhol/Desktop/NGTRCode_GIT/exNGT" #Maria
path <- "~/programming/R/exNGT" #Thomas

setwd(path)
source("init.R")

# ------------------------------------------------------------------------------
# Load data

NGT <- processNGT() %>%
  stackNGT()

Arctic2k <- readArctic2k() %>%
  extendWithHadCrut()

TCR <- readTwenCR()

# ------------------------------------------------------------------------------
# Cut out range of overlapping years to calculate coherence

t1 <- 2011 : 1000

ngt1 <- subsetData(NGT, t = t1, var = "stack")
a2k1 <- subsetData(Arctic2k, t = t1, var = "TempAnomaly")

t2 <- 2011 : 1836

ngt2 <- subsetData(NGT, t = t2, var = "stack")
tcr2 <- subsetData(TCR, t = t2, var = "t2m")

# ------------------------------------------------------------------------------
# Estimate mangitude-squared coherence with significance

nmc <- 1000

ngt2a2k <- coherence(ngt1, a2k1, spans = 61, nmc = nmc)
ngt2tcr <- coherence(ngt2, tcr2, spans = 11, nmc = nmc)

# ------------------------------------------------------------------------------
# Save for use in main plot

list(ngt2a2k = ngt2a2k, ngt2tcr = ngt2tcr) %>%
  saveRDS(file = "out/coherence.rds")
