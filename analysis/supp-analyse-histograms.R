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

asp <- 2 / 3
w <- 18
natfig(file = "./fig/supplement-histograms.pdf",
       height = asp * w, width = w,
       mfrow = c(2, 3), mar = c(5, 5, 4, 0.5))
par(cex = 1)

# 1. Main method
plotHistogram(stack.method = "main", filter.window = filter.window,
              plot.legend = FALSE, xmain = xmain[1], q.lty = 2)

# 2. Main stack after full forward diffusion
plotHistogram(stack.method = "main", filter.window = filter.window,
              diffuse = TRUE, plot.legend = FALSE, xmain = xmain[2], q.lty = 2)

# 3. Main stack for constant number of records (n = 5)
plotHistogram(stack.method = "fix_N", filter.window = filter.window,
              nfix = 5, plot.legend = FALSE, xmain = xmain[3], q.lty = 2)

# 4. Simply stack all records
plotHistogram(stack.method = "stack_all", filter.window = filter.window,
              plot.legend = FALSE, xmain = xmain[4], q.lty = 2)

# 5. Simply stack with full forward diffusion
plotHistogram(stack.method = "stack_all", filter.window = filter.window,
              diffuse = TRUE, plot.legend = FALSE, xmain = xmain[5], q.lty = 2)

op <- par(cex = 1, mar = c(0, 0, 0, 0))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")

lg <- c("Pre-industrial distribution\n(1000-1800 CE)", "Gaussian fit\n",
        "p = 0.95 and p = 0.99\n", "2001-2011 CE\naverage")

legend("topleft", legend = lg, lty = c(1, 1, 2, 1), lwd = Lwd(c(10, 2, 1, 3)),
       col = c(adjustcolor("black", 0.2), "black", "black", "orange2"),
       bty = "n", adj = c(0, 0.8), y.intersp = 1.5, seg.len = 2.5,
       inset = c(0.2, -0.05), cex = 7 / 6)

par(op)
dev.off()
