\name{bcf.next}
\alias{bcf.next}
\title{bcf.next}
\usage{bcf.next(fp)}
\description{read the next Variant Context in the VCF reader


fp <- bcf.open("in.vcf.gz")
bcf.query(fp,"RF02:1-1000")
while(!is.null(vc<-bcf.next(fp))) {
## do something with vc
}
bcf.close(fp)
}
\value{@return an opaque vcf context or NULL }
\arguments{
	\item{fp}{the vcf reader}
}
