# Show Htslib and Rbcf versions
# load the library
library(rbcf)
#print the version of the associated htslib 
paste("HTSLIB:",htslib.version())
#print the version of rbcf
paste("RBCF:",rcbf.version())

