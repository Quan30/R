
# change here the package name
install.packages("tictoc", lib=my_lib_path)  

# add the path to be recognized by R
my_lib_path <- "~/tools/R/library"
.libPaths(c(my_lib_path, .libPaths()))

# call the package, see if it's working
library(tictoc)