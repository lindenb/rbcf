## Setting FORMAT in variant genotypes
# load rbcf
library(rbcf)
# vcf input filename
filenamein = "./data/1000G.ALL.2of4intersection.20100804.genotypes.bcf"
# output vcf filename. "-" is standard output
filenameout =  "-"

fp <- bcf.open(filenamein,FALSE)
# error on opening (exit 0 for tests)
if(is.null(fp)) quit(save="no",status=0,runLast=FALSE)

# create a new VCF writer using the header from 'fp'
out <- bcf.new.writer(fp,filenameout)
# error on opening (exit 0 for tests)
if(is.null(out)) quit(save="no",status=0,runLast=FALSE)

n_samples <- bcf.nsamples(fp)

# loop while we can read a variant
while(!is.null(vc<-bcf.next(fp))) {
  # Setting of INT attributes
  old_dp <- variant.genotypes.int.attribute(vc, "DP")
  new_dp <- rnbinom(n = n_samples, size = 5, prob = 0.25)
  variant.genotypes.set.int.attribute(vc, "DP", new_dp)
  stopifnot( all(  new_dp == variant.genotypes.int.attribute(vc, "DP") ))

  old_gq <- variant.genotypes.float.attribute(vc, "GQ")
  new_gq <- rnorm(n = n_samples)
  variant.genotypes.set.float.attribute(vc, "GQ", new_gq)
  stopifnot( all(  new_gq - variant.genotypes.float.attribute(vc, "GQ") <  0.000001 ))
  
  # Updating GT with random allele combinations
  alleles <- c(".", as.character(0:variant.nalleles(vc)))
  all_potential_genotypes = c(".",
    apply(as.matrix(expand.grid(alleles, alleles)), 1, paste, collapse = "/"),
    apply(as.matrix(expand.grid(0:variant.nalleles(vc), 0:variant.nalleles(vc))), 1, paste, collapse = "|")
  )
  
  new_gt <- sample(x = all_potential_genotypes, size = n_samples, replace = TRUE)
  variant.genotypes.set.allele.strings(vc, new_gt)
  new_gt <- sub("-", ".", new_gt)
  subset(data.frame( Expected = new_gt, Seen = variant.genotypes.allele.strings(vc)), Expected != Seen)
  
  stopifnot(all(  new_gt == variant.genotypes.allele.strings(vc) ))
  
  bcf.write.variant(out, vc)
}
# dispose the vcf reader
bcf.close(fp)
# dispose the vcf rwriter
bcf.close(out)
