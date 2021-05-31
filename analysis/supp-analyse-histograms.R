##
## aim:
## analyse the histograms for different stacking possibilities and options.
## relation:
## NGT paper supplementary; https://github.com/EarthSystemDiagnostics/exNGT
##

path <- "C:/Users/mhoerhol/Desktop/NGTRCode_GIT/exNGT" #Maria
path <- "~/programming/R/exNGT" #Thomas

setwd(path)
source("init.R")

# ------------------------------------------------------------------------------
# Compare main method to most relevant other options in a plot

filter.window <- 11

xmain = c("a. NGT main stack",
          "b. NGT main stack (forward diffused)",
          "c. NGT \"frozen\" stack (fixed N)",
          "d. NGT simple stack",
          "e. NGT simple stack (forward diffused)")

grfxtools::Quartz(file = "./fig/supplement-histograms.pdf",
                  height = 8, width = 12,
                  mfrow = c(2, 3), mar = c(5, 5, 4, 0.5))

# 1. Main method
plotHistogram(stack.method = "main", filter.window = filter.window,
              plot.legend = FALSE, xmain = xmain[1])

# 2. Main stack after full forward diffusion
plotHistogram(stack.method = "main", filter.window = filter.window,
              diffuse = TRUE, plot.legend = FALSE, xmain = xmain[2])

# 3. Main stack for constant number of records (n = 5)
plotHistogram(stack.method = "fix_N", filter.window = filter.window,
              nfix = 5, plot.legend = FALSE, xmain = xmain[3])

# 4. Simply stack all records
plotHistogram(stack.method = "stack_all", filter.window = filter.window,
              plot.legend = FALSE, xmain = xmain[4])

# 5. Simply stack with full forward diffusion
plotHistogram(stack.method = "stack_all", filter.window = filter.window,
              diffuse = TRUE, plot.legend = FALSE, xmain = xmain[5])

op <- par(cex = 1, mar = c(0, 0, 0, 0))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")

lg <- c("Pre-industrial distribution\n(1000-1800 CE)", "Gaussian fit",
        "2.5, 50, 97.5 %\nquantiles", "Recent value\n(2001-2011 CE)")

legend("topleft", legend = lg, lty = c(1, 1, 5, 1), lwd = c(10, 2, 1, 2.5),
       col = c(adjustcolor("black", 0.2), "black", "black", "firebrick4"),
       bty = "n", adj = c(0, 0.75), y.intersp = 1.5, seg.len = 2.5,
       inset = c(0.1, -0.05))

par(op)
dev.off()
