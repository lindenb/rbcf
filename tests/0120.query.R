##Query the indexed vcf using intervals
# load rbcf
library(rbcf)

# create a function counting variants in a VCF, in some intervals
count.variants<-function(filename,intervals) {
	# open the indexed VCF
	fp <- bcf.open(filename)
	# error on opening
	if(is.null(fp)) return(-1)
	# loop over the intervals
	for(interval in intervals) {
		# try query the interval
		if(bcf.query(fp,interval)) {
			# number of variants
			n<-0
			# loop while we can read a variant
			while(!is.null(vc<-bcf.next(fp))) {
				# increment the count
				n<-n+1
				}
			print(paste("Number of variants in ",basename(filename),"/'",interval,"' :",n,sep=""))
			}
		# query failed
		else {
			print(paste("Cannot query ",basename(filename),"/'",interval,"'",sep=""))
			}
		}
	# dispose the vcf reader
	bcf.close(fp)
}

some_intervals <-c("","RF03","RF03:2000-3000","1:1-10000000","chr1")
count.variants("./data/rotavirus_rf.02.vcf.gz",some_intervals)
count.variants("./data/1000G.ALL.2of4intersection.20100804.genotypes.bcf",some_intervals)

# another way to query is set collect=TRUE to return a vector of variant
fp <- bcf.open("./data/rotavirus_rf.02.vcf.gz")
if(!is.null(fp)) {
	print(paste("Number of variants using collect:",length(bcf.query(fp,"RF03",collect=TRUE))))
	bcf.close(fp)
	}

