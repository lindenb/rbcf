\name{bcf.open}
\alias{bcf.open}
\title{bcf.open}
\usage{bcf.open(filename,requireIndex=TRUE)}
\description{Open a VCF or a BCF file}
\value{@return the new VCF reader}
\arguments{
	\item{filename}{the path to the vcf file}
	\item{requireIndex}{load associated vcf index}
}
\examples{

bcf.open("my.vcf.gz",T)
bcf.open("my.bcf")

}
