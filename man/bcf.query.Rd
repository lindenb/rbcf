\name{bcf.query}
\alias{bcf.query}
\title{bcf.query}
\usage{bcf.query(fp,interval,collect=FALSE)}
\description{prepare the VCF reader for a new vcf iteration over a given interval.
VCF reader must be associated and opened with a valid index.}
\value{@return TRUE on success or a list of variant if collect=TRUE}
\arguments{
	\item{fp}{the vcf reader}
	\item{interval}{the genomic interval to be scanned}
	\item{collect}{if TRUE return a list of variants in the region or NULL on failure}
}
