\name{bcf.dictionary}
\alias{bcf.dictionary}
\title{bcf.dictionary}
\usage{bcf.dictionary(fp)}
\description{get dictionary from a VCF reader.


fp<-bcf.open("my.vcf.gz",TRUE)
dict<-bcf.dictionary(fp)
print(dict)
chrom size
RF01  RF01 3302
RF02  RF02 2687
RF03  RF03 2592
RF04  RF04 2362
RF05  RF05 1579
RF06  RF06 1356
RF07  RF07 1074
RF08  RF08 1059
RF09  RF09 1062
RF10  RF10  751
RF11  RF11  666
}
\value{@return the vcf dictionary as a table}
\arguments{
	\item{fp}{the vcf reader}
}
