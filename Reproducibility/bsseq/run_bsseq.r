library(optparse)
library(DSS)
require(bsseq)


option_list <- list(
  make_option(
    c("-a", "--case"),
    type = "character",
    default = NULL,
    help = "A comma-separated paths to case (e.g., 'val1,val2,val3')",
    metavar = "character"
  ),
  make_option(
    c("-o", "--control"),
    type = "character",
    default = NULL,
    help = "A comma-separated paths to control (e.g., 'val1,val2,val3')",
    metavar = "character"
  ),
  make_option(
    c("-s", "--assembly"),
    type = "character",
    default = "hg38",
    help = "reference genome",
    metavar = "character"
  ),
  make_option(
    c("-u", "--output"),
    type = "character",
    default = "",
    help = "Output folder",
    metavar = "character"
  ),
  make_option(
    c("-w", "--cov1"),
    type = "integer",
    default = 5,
    help = "Minimum coverage for case",
    metavar = "integer"
  ),
  make_option(
    c("-x", "--spcov1"),
    type = "integer",
    default = 2,
    help = "Minimum sample coverage for case",
    metavar = "integer"
  ),
  make_option(
    c("-y", "--cov2"),
    type = "integer",
    default = 10,
    help = "Minimum coverage for control",
    metavar = "integer"
  ),
  make_option(
    c("-z", "--spcov2"),
    type = "integer",
    default = 2,
    help = "Minimum sample coverage for control",
    metavar = "integer"
  )

#  make_option(
#    c("-c", "--cutoff"),
#    type = "float",
#    default = 4.6,
#    help = "The cutoff for t-test",
#    metavar = "float"
#  )

)


opt_parser <- OptionParser(option_list = option_list)

opt <- parse_args(opt_parser)

if (!is.null(opt$case)) {
  case <- strsplit(opt$case, ",")[[1]]
} else {
  cat("No values provided.\n")
}


if (!is.null(opt$control)) {
  control <- strsplit(opt$control, ",")[[1]]
} else {
  cat("No values provided.\n")
}



cov1 <- opt$cov1
cov2 <- opt$cov2

spcov1 <- opt$spcov1
spcov2 <- opt$spcov2


case_list <- strsplit(case, ",")
control_list <- strsplit(control, ",")


class_1 <- paste0("case", 1:length(case_list))
class_2 <- paste0("control", 1:length(control_list))

type_1 <- rep("case", length(case_list))
type_2 <- rep("control", length(control_list))

files <- unlist(list(case_list, control_list)) # c(case_list, control_list) # unlist(list(case_list, control_list))


classes <- unlist(list(class_1, class_2))
types <- unlist(list(type_1, type_2))

# treatment <- unlist(list(treatment_1, treatment_2))
#files <- as.list(files)

# print(files)
# print(files)

# class_1 <- as.list(class_1)
# class_2 <- as.list(class_2)
classes <- c(as.list(class_1), as.list(class_2))

# print(unlist(files))
# print(unlist(classes))

# print(classes)
# print(types)

# print(type_1)
# print(type_2)
print(classes)
print(types)




bismarkBSseq <- read.bismark(files = files,
                               rmZeroCov = FALSE,
                               strandCollapse = FALSE,
                               verbose = TRUE)


sampleNames(bismarkBSseq) <- classes

bismarkBSseq.fit <- BSmooth(
    BSseq = bismarkBSseq,
    verbose = TRUE)

bismarkBSseq.fit <- updateObject(bismarkBSseq.fit)

sample_info <- data.frame(
  SampleNames = unlist(classes),
  SampleType = unlist(types)
)

colData(bismarkBSseq.fit) <- DataFrame(SampleType = sample_info$SampleType)

colData(bismarkBSseq) <- DataFrame(SampleType = sample_info$SampleType)

bismarkBSseq.cov <- getCoverage(bismarkBSseq.fit)

keepLoci.ex <- which(rowSums(bismarkBSseq.cov[, bismarkBSseq$SampleType == type_1[0]] >= cov1) >= spcov1 &
                     rowSums(bismarkBSseq.cov[, bismarkBSseq$SampleType == type_2[0]] >= cov2) >= spcov2)


sampleNames(bismarkBSseq.fit) <- sample_info$SampleNames

# print(sampleNames(bismarkBSseq.fit))
print(class_2)
print(class_1)


bismarkBSseq.tstat <- BSmooth.tstat(bismarkBSseq.fit,
                                    group1 = class_2,
                                    group2 = class_1,
                                    estimate.var = "group2",
                                    local.correct = TRUE,
                                    verbose = TRUE,
                                    mc.cores= 18)

print(bismarkBSseq.tstat)

dmrs0 <- dmrFinder(bismarkBSseq.tstat, cutoff = c(-4.6, 4.6))

write.csv(dmrs0, paste(opt$output,"/","new.dsseq.",opt$context ,".",opt$assembly ,".","DMR" , ".csv",  sep="" ), row.names = FALSE)

