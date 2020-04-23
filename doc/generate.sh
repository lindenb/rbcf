#/bin/bash

set -e

cat << __EOF__ > tmp.tex
\documentclass[12pt, a4paper,notitlepage,onecolumn]{article}
\usepackage[utf8]{inputenc}
\usepackage{listings}
\title{RBcf: An VCF API for R.}

\author{Pierre Lindenbaum / yokofakun/ Institut du Thorax . Nantes.}
\date{\today}
\begin{document}
\maketitle

\section{Abstract}
RBcf uses the Htslib C API for parsing VCF and BCF files.
This API was written by a regular user of the htsjdk library who doesn't like R.

\section{Examples}
__EOF__

ls *.R | sort | while read F
do
	echo -n '\subsection{' >> tmp.tex
	head -n1 "${F}" | cut -c 3- | sed 's/$/}/' >>  tmp.tex
	echo -e 'Code:\n\\begin{verbatim}' >> tmp.tex
	tail -n+2 "${F}" >> tmp.tex
	echo -e '\\end{verbatim}\nOutput:\n\\begin{verbatim}' >> tmp.tex
	tail -n+2 "${F}" | R --no-save --quiet  >> tmp.tex
	echo -e '\\end{verbatim}\n'  >> tmp.tex

done


cat << __EOF__ >> tmp.tex
\end{document}
__EOF__

pdflatex tmp.tex
mv tmp.pdf rbcf_examples.pdf
rm tmp.tex

