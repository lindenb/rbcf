\name{bcf.new.writer}
\alias{bcf.new.writer}
\title{bcf.new.writer}
\usage{bcf.new.writer(fp,fname="-")}
\description{Open a new VCF writer.
Must be closed with bcf.close
}
\value{@return the writer or NULL on failure}
\arguments{
	\item{fp}{the vcf reader}
	\item{fname}{the name of the output vcf file}
}
