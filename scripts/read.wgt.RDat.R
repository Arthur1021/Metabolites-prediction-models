args = commandArgs(trailingOnly=TRUE);

# test argument number
if (length(args) !=  2) {
   stop("Two and only two arguments are required: input.RDat output.list", call.=FALSE);
} else {
  f1 = args[1];
  f2 = args[2];
}

# read/load RDat to obtain objects in the data: cv.performance, wgt.matrix, etc
load(f1);

#
# > cv.performance
#              blup        lasso         top1         enet
# rsq  5.824996e-02 7.758592e-02 7.059922e-02 7.967104e-02
# pval 1.968176e-34 1.162418e-45 1.404434e-41 6.941667e-47

# > head(wgt.matrix)
#                     blup lasso       top1 enet
# rs12570924 -0.0001140030     0 -1.5895390    0
# rs72772882  0.0008071184     0  1.2651891    0
# rs10795122  0.0001380190     0 -0.4289210    0
# rs10795124 -0.0009524297     0 -0.7788182    0

# get the max of rows, and 'rsq' row
m=(colnames(cv.performance)[apply(cv.performance, 1, which.min)])[2] #lowest p-value, if not unique, will select first.
#m=(colnames(cv.performance)[apply(cv.performance, 1, which.max)])[1]; 

# convert to DF for easy print of IDs
wm = as.data.frame(which(wgt.matrix[,m] != 0));
# wm <- row.names(wgt.matrix)
sink(f2);
cat(rownames(wm), sep="\n");
# cat(wm, sep="\n");
sink();