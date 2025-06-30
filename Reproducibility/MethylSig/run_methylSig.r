library(optparse)
# library(DSS)
require(bsseq)
library(methylSig)



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
    c("-c", "--cov"),
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
    c("-z", "--spcov2"),
    type = "integer",
    default = 2,
    help = "Minimum sample coverage for control",
    metavar = "integer"
  ),
  make_option(
    c("-i", "--window_size"),
    type = "integer",
    default = 1000,
    help = "Window size if window based is chosen",
    metavar = "integer"
  ),
  make_option(
    c("-w", "--window"),
    type = "logical",
    default = FALSE,
    help = "Window or position based",
    metavar = "logical"
  )
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



cov <- opt$cov

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


bsseq::pData(bismarkBSseq)


bs = filter_loci_by_coverage(bismarkBSseq, min_count = cov)

sample_info <- data.frame(
  SampleNames = unlist(classes),
  SampleType = unlist(types)
)


colData(bs) <- DataFrame(SampleType = sample_info$SampleType)

sample_info <- data.frame(
   Pair = unlist(classes),
   Type = unlist(types)
 )

rownames(sample_info) <- unlist(classes)

# bsseq::pData(bs) <- sample_info



bs = filter_loci_by_group_coverage(
    bs = bs,
    group_column = 'SampleType',
    c("case" = spcov1, 'control' = spcov2))


win <- "non.window"


print(colData(bs))
print("--------------------------------")
print(granges(bs))
print("--------------------------------")

print(nrow(colData(bs)))
print(ncol(colData(bs)))
print(dim(colData(bs)))

pData(bs) <- sample_info

if (opt$window){
	# pData(bs) <- sample_info # DataFrame(SampleType = sample_info$SampleType)
        ws <- opt$window_size
        bs = tile_by_windows(bs = bs, win_size = ws)
        win <- paste("window", as.character(ws), sep=".")
}



diff_gr = diff_methylsig(
    bs = bs,
    group_column = 'Type',
    comparison_groups = c('case' = 'case', 'control' = 'control'),
    disp_groups = c('case' = TRUE, 'control' = TRUE),
    local_window_size = 0,
    t_approx = FALSE, #,TRUE,
    n_cores = 20)

write.csv(diff_gr, paste(opt$output,"/","new.methylSig.",opt$context ,".",opt$assembly ,".",win , ".csv",  sep="" ), row.names = FALSE)

#print(paste(opt$output,"/","new.methylkit",opt$context ,opt$assembly ,win ,"cov",c, ".csv",  sep="" ))

#meth2_2=unite(myobj1, min.per.group=m, destrand=FALSE, mc.cores=20)
#myDiff2=calculateDiffMeth(meth2_2, mc.cores=20)
# myDiff25p2=getMethylDiff(myDiff2,difference=25,qvalue=0.05)


# write.csv(dmls, paste(opt$output,"/","new.dss.",ds ,opt$context ,".",opt$assembly ,"DML" , ".csv",  sep="" ), row.names = FALSE)

# write.csv(dmrs, paste(opt$output,"/","new.dss.",ds ,opt$context ,".",opt$assembly ,"DMR" , ".csv",  sep="" ), row.names = FALSE)
