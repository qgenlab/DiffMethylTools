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
    c("-n", "--context"),
    type = "character",
    default = "CpG",
    help = "CpG or non CpG",
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
    c("-d", "--destranded"),
    type = "logical",
    default = TRUE,
    help = "Destrand the input or not",
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

ds <- ""

if (!opt$destranded){
	ds <- "stranded."
}


case_list <- strsplit(case, ",")
control_list <- strsplit(control, ",")


class_1 <- paste0("case", 1:length(case_list))
class_2 <- paste0("control", 1:length(control_list))


files <- unlist(list(case_list, control_list)) # c(case_list, control_list) # unlist(list(case_list, control_list))


classes <- list(unlist(list(class_1, class_2)))

# treatment <- unlist(list(treatment_1, treatment_2))


files <- lapply(files, file.path)

files <- lapply(files, read.table, header = TRUE)
#files <- as.list(files)

# print(files)
# print(files)

# class_1 <- as.list(class_1)
# class_2 <- as.list(class_2)
classes <- c(as.list(class_1), as.list(class_2))

print(class_1)
print(class_2)

BSobj = makeBSseqData(files, classes )

print(BSobj)

dmlTest = DMLtest(BSobj, group1=class_1, group2=class_2, ncores=20) # ncores=20


dmls = callDML(dmlTest, delta=0.0, p.threshold=1)


dmrs = callDMR(dmlTest, delta=0.1, p.threshold=0.05)
#print(paste(opt$output,"/","new.methylkit",opt$context ,opt$assembly ,win ,"cov",c, ".csv",  sep="" ))

#meth2_2=unite(myobj1, min.per.group=m, destrand=FALSE, mc.cores=20)
#myDiff2=calculateDiffMeth(meth2_2, mc.cores=20)
# myDiff25p2=getMethylDiff(myDiff2,difference=25,qvalue=0.05)


write.csv(dmls, paste(opt$output,"/","new.dss.",ds ,opt$context ,".",opt$assembly ,"DML" , ".csv",  sep="" ), row.names = FALSE)

write.csv(dmrs, paste(opt$output,"/","new.dss.",ds ,opt$context ,".",opt$assembly ,"DMR" , ".csv",  sep="" ), row.names = FALSE)
