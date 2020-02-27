path <- "C:/Users/mhoerhol/Desktop/NGTRCode_GIT/exNGT" #Maria
#path <- "~/programming/R/exNGT" #Thomas

#setwd(path)

# Load all exNGT library functions
source("init.R")



NGT <- readNGT()

Arctic2k <- readArctic2k()

# Get anomaly
NGT <- processNGT()

#Get stack
stackID_12 <- c("Year","B18_12", "B21_12", "B23_12", "B26_12", "NGRIP_12")
stackdata_12 <-NGT[,stackID_12]
stack_12 <- rowMeans(stackdata_12[,2:6], na.rm=TRUE)
stackID_94 <- c("Year","B18", "B21", "B23", "B26", "NGRIP")
stackdata_94 <-NGT[,stackID_94]
stack <- rowMeans(stackdata_94[,2:6], na.rm=TRUE)

#Add stacks to NGT
NGT_ext <- cbind(NGT, stack, stack_12)
write.csv(NGT_ext, file="NGT_ext.csv")
#read into Igor as delimited text - will skip the first column
#Opening in Excel, replace all NA as NaN, then Igor can read all columns

#Apply Filter
NGT_3filt <- filterData(NGT_ext, window=3)
write.csv(NGT_3filt, file="NGT_3filt.csv")
NGT_5filt <- filterData(NGT_ext, window=5)
write.csv(NGT_5filt, file="NGT_5filt.csv")
NGT_7filt <- filterData(NGT_ext, window=7)
write.csv(NGT_7filt, file="NGT_7filt.csv")
NGT_11filt <- filterData(NGT_ext, window=11)
write.csv(NGT_11filt, file="NGT_11filt.csv")
NGT_21filt <- filterData(NGT_ext, window=21)
write.csv(NGT_21filt, file="NGT_21filt.csv")

#Calculate Overlap statistics

NGT_overlap1_B18 <- calculateOverlapStatistics(NGT_ext, site = "B18")
NGT_overlap3_B18 <- calculateOverlapStatistics(NGT_3filt, site = "B18")
NGT_overlap5_B18 <- calculateOverlapStatistics(NGT_5filt, site = "B18")
NGT_overlap7_B18 <- calculateOverlapStatistics(NGT_7filt, site = "B18")
NGT_overlap11_B18 <- calculateOverlapStatistics(NGT_11filt, site = "B18")
NGT_overlap21_B18 <- calculateOverlapStatistics(NGT_21filt, site = "B18")

NGT_overlap1_B21 <- calculateOverlapStatistics(NGT_ext, site = "B21")
NGT_overlap3_B21 <- calculateOverlapStatistics(NGT_3filt, site = "B21")
NGT_overlap5_B21 <- calculateOverlapStatistics(NGT_5filt, site = "B21")
NGT_overlap7_B21 <- calculateOverlapStatistics(NGT_7filt, site = "B21")
NGT_overlap11_B21 <- calculateOverlapStatistics(NGT_11filt, site = "B21")
NGT_overlap21_B21 <- calculateOverlapStatistics(NGT_21filt, site = "B21")

NGT_overlap1_B23 <- calculateOverlapStatistics(NGT_ext, site = "B23")
NGT_overlap3_B23 <- calculateOverlapStatistics(NGT_3filt, site = "B23")
NGT_overlap5_B23 <- calculateOverlapStatistics(NGT_5filt, site = "B23")
NGT_overlap7_B23 <- calculateOverlapStatistics(NGT_7filt, site = "B23")
NGT_overlap11_B23 <- calculateOverlapStatistics(NGT_11filt, site = "B23")
NGT_overlap21_B23 <- calculateOverlapStatistics(NGT_21filt, site = "B23")

NGT_overlap1_B26 <- calculateOverlapStatistics(NGT_ext, site = "B26")
NGT_overlap3_B26 <- calculateOverlapStatistics(NGT_3filt, site = "B26")
NGT_overlap5_B26 <- calculateOverlapStatistics(NGT_5filt, site = "B26")
NGT_overlap7_B26 <- calculateOverlapStatistics(NGT_7filt, site = "B26")
NGT_overlap11_B26 <- calculateOverlapStatistics(NGT_11filt, site = "B26")
NGT_overlap21_B26 <- calculateOverlapStatistics(NGT_21filt, site = "B26")

NGT_overlap1_NGRIP <- calculateOverlapStatistics(NGT_ext, site = "NGRIP")
NGT_overlap3_NGRIP <- calculateOverlapStatistics(NGT_3filt, site = "NGRIP")
NGT_overlap5_NGRIP <- calculateOverlapStatistics(NGT_5filt, site = "NGRIP")
NGT_overlap7_NGRIP <- calculateOverlapStatistics(NGT_7filt, site = "NGRIP")
NGT_overlap11_NGRIP <- calculateOverlapStatistics(NGT_11filt, site = "NGRIP")
NGT_overlap21_NGRIP <- calculateOverlapStatistics(NGT_21filt, site = "NGRIP")

NGT_overlap1_stack <- calculateOverlapStatistics(NGT_ext, site = "stack")
NGT_overlap3_stack <- calculateOverlapStatistics(NGT_3filt, site = "stack")
NGT_overlap5_stack <- calculateOverlapStatistics(NGT_5filt, site = "stack")
NGT_overlap7_stack <- calculateOverlapStatistics(NGT_7filt, site = "stack")
NGT_overlap11_stack <- calculateOverlapStatistics(NGT_11filt, site = "stack")
NGT_overlap21_stack <- calculateOverlapStatistics(NGT_21filt, site = "stack")

#-------------------------------------
#Continue with filter length 11
#extract signle drill sites
stackID_single <- c("Year","B16","B17","B20", "B22","B27.28","B29", "B30","GISP2", "GRIP", "NEGIS","NEEM")
stackdata_single <-NGT_11filt[,stackID_single]

#Merge re-drilled sites with different approaches
NGT11_red_AdjMean_Start <- mergeCores(NGT_11filt, method=1, adjustMean=TRUE)
NGT11_red_AdjMean_End <- mergeCores(NGT_11filt, method=2, adjustMean=TRUE)
NGT11_red_NoAdjMean_Start <- mergeCores(NGT_11filt, method=1, adjustMean=FALSE)
NGT11_red_NoAdjMean_End <- mergeCores(NGT_11filt, method=2, adjustMean=FALSE)

#constitute the different stacks, note that last column is stack from beginning
NGT11_NoAdjMean_start <-cbind(stackdata_single, NGT11_red_NoAdjMean_Start) 
NGT11_NoAdjMean_end <-cbind(stackdata_single, NGT11_red_NoAdjMean_End) 
NGT11_AdjMean_start <-cbind(stackdata_single, NGT11_red_AdjMean_Start) 
NGT11_AdjMean_end <-cbind(stackdata_single, NGT11_red_AdjMean_End)

#rowlmean without last column, last column = stack

stack_NoAdjMean_start <- rowMeans(NGT11_NoAdjMean_start[,2:17], na.rm=TRUE) #stack 1
stack_NoAdjMean_end <- rowMeans(NGT11_NoAdjMean_end[,2:17], na.rm=TRUE) #stack 2
stack_AdjMean_start <- rowMeans(NGT11_AdjMean_start[,2:17], na.rm=TRUE) #stack 3
stack_AdjMean_end <- rowMeans(NGT11_AdjMean_end[,2:17], na.rm=TRUE) #stack 4

#last column with different merging approaches
stackS1_NoAdjMean_start<-NGT11_NoAdjMean_start[,18] #stack 5
stackS1_NoAdjMean_end<-NGT11_NoAdjMean_end[,18] #stack 6
stackS1_AdjMean_start<-NGT11_AdjMean_start[,18] #stack 7
stackS1_AdjMean_end<-NGT11_AdjMean_end[,18] #stack 8

# stacking over all single cores
stackAllSingle <-rowMeans(NGT[,2:22], na.rm=TRUE)  #stack 9
#stackAllSingle11 <-filterData(stackAllSingle, window=11)

yearID <-("Year")
yeardata <- NGT[yearID]
stackTable <- cbind.data.frame(yeardata, stack_NoAdjMean_start, stack_NoAdjMean_end, stack_AdjMean_start, stack_AdjMean_end, stackS1_NoAdjMean_start, stackS1_NoAdjMean_end, stackS1_AdjMean_start, stackS1_AdjMean_end, stackAllSingle)
write.csv(stackTable, file="stacktable.csv")

