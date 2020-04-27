\name{bcf.open}
\alias{bcf.open}
\title{bcf.open}
\usage{bcf.open(filename,requireIndex=TRUE)}
\description{Open a VCF or a BCF file

fp<-bcf.open("my.vcf.gz",TRUE)
hts.close(fp)
fp<-bcf.open("my.bcf")
hts.close(fp)
}
\value{@return the new VCF reader}
\arguments{
	\item{filename}{the path to the vcf file}
	\item{requireIndex}{load associated vcf index}
}
