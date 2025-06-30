library(optparse)
library(methylKit)

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
    c("-w", "--window"),
    type = "logical",
    default = FALSE,
    help = "Window or position based",
    metavar = "logical"
  ),
  make_option(
    c("-c", "--coverage"),
    type = "integer",
    default = 10L,
    help = "Minimum coverage per position",
    metavar = "integer"
  ),
  make_option(
    c("-m", "--min_group"),
    type = "integer",
    default = 2L,
    help = "Minimum number of samples per group",
    metavar = "integer"
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
    c("-i", "--window_size"),
    type = "integer",
    default = 1000,
    help = "Window size if window based is chosen",
    metavar = "integer"
  ),
  make_option(
    c("-t", "--step"),
    type = "integer",
    default = 500,
    help = "Step size if window based is chosen",
    metavar = "integer"
  ),
  make_option(
    c("-u", "--output"),
    type = "character",
    default = "",
    help = "Output folder",
    metavar = "character"
  ),
  make_option(
    c("-d", "--destrand"),
    type = "logical",
    default = FALSE,
    help = "Destrand the input or not",
    metavar = "logical"
   ),
  make_option(
    c("-b", "--bismark"),
    type = "integer",
    default = 0,
    help = "Is in input file bed (0), bismarkCytosineReport (1) or a bedGraph (any)",
    metavar = "integer"
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


case_list <- strsplit(case, ",")
control_list <- strsplit(control, ",")


class_1 <- rep("case", length(case_list))
class_2 <- rep("control", length(control_list))

treatment_1 <- rep(0, length(case_list))
treatment_2 <- rep(1, length(control_list))

c <- as.integer(opt$coverage)
m <- as.integer(opt$min_group)

files <- unlist(list(case_list, control_list))

#print(files)

# classes <- list(unlist(list(class_1, class_2)))
treatment <- unlist(list(treatment_1, treatment_2))

files <- list(lapply(files, file.path))
#files <- as.list(files)

class_1 <- as.list(class_1)
class_2 <- as.list(class_2)
classes <- c(class_1, class_2)

input_f <- opt$bismark

print("/********************************/ min cov is ")
print(c)
print("/********************************/ min sample is ")
print(m)

if (input_f == 1) {
myobj1=methRead(files[[1]],
           sample.id=classes,
           assembly=opt$assembly,
           treatment=treatment,
           context=opt$context,
                   pipeline = "bismarkCytosineReport",
                   mincov=c
           )
} else if (input_f == 0) {
myobj1=methRead(files[[1]],
           sample.id=classes,
           assembly=opt$assembly,
           treatment=treatment,
                   context=opt$context,
                   pipeline = list(fraction=FALSE,chr.col=1,start.col=2,end.col=2,coverage.col=10,strand.col=6,freqC.col=11),
                   mincov=c
	)
} else {
myobj1=methRead(files[[1]],
           sample.id=classes,
           assembly=opt$assembly,
           treatment=treatment,
                   context=opt$context,
                   pipeline = list(fraction=FALSE,chr.col=2,start.col=3,end.col=3,coverage.col=5,strand.col=4,freqC.col=6),
                   mincov=c
)
}
win <- "non.window"

if (opt$window){
	ws <- opt$window_size
	step <- opt$step
	myobj1 = tileMethylCounts(myobj1,win.size=ws,step.size=step,cov.bases = 10)
	win <- paste("window", as.character(ws),"step", as.character(step), sep=".")
}

destrand <- ""
ds <- opt$destrand

if (ds){
        destrand <- "destranded."
}



print(paste(opt$output,"/","new.methylkit",opt$context ,opt$assembly ,win ,"cov",c, ".csv",  sep="" ))

meth2_2=unite(myobj1, min.per.group=m, destrand=ds, mc.cores=20) # mc.cores=20
print(meth2_2)

#write.csv(meth2_2, "test.csv")
myDiff2=calculateDiffMeth(meth2_2, mc.cores=20)
# myDiff25p2=getMethylDiff(myDiff2,difference=25,qvalue=0.05)


write.csv(myDiff2, paste(opt$output,"/","new.methylkit.",".",destrand,opt$context ,".",opt$assembly ,".",win ,".cov.",c, ".csv",  sep="" ), row.names = FALSE)
