\name{bcf.close}
\alias{bcf.close}
\title{bcf.close}
\usage{bcf.close(fp)}
\description{Close a VCF reader}
\value{@return true on success}
\arguments{
	\item{fp}{the vcf reader}
}
\examples{

fp<-bcf.open("my.vcf.gz",T)
bcf.close(fp)

}
