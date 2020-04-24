#/bin/bash

set -e

cat << __EOF__ > tmp.tex
\documentclass[12pt, a4paper,notitlepage,onecolumn]{article}
\usepackage[utf8]{inputenc}
\usepackage{xcolor}
\usepackage{hyperref}
\usepackage{listings}
\title{RBcf: An VCF API for R.}

\lstset{
  backgroundcolor = \color{lightgray},
  language = R,
  basicstyle=\ttfamily,
  columns=fullflexible,
  keepspaces=true,
}

\author{Pierre Lindenbaum / @yokofakun/ Institut du Thorax . Nantes.}
\date{\today}
\begin{document}
\maketitle

\section{Abstract}
RBcf uses the \href{https://github.com/samtools/htslib}{Htslib C API} for parsing VCF and BCF files.
This API was written by a regular user of the \href{https://github.com/samtools/htsjdk}{htsjdk} library who doesn't like R.

A list of functions is available at: \href{https://github.com/lindenb/rbcf/blob/master/R/rbcf.R}{https://github.com/lindenb/rbcf/blob/master/R/rbcf.R}

\section{Examples}
__EOF__

ls *.R | sort | while read F
do
	echo -n '\subsection{' >> tmp.tex
	head -n1 "${F}" | cut -c 3- | sed 's/$/}/' >>  tmp.tex
	tail -n+2 "${F}" | grep  "^#'" | cut -c 3-  >> tmp.tex
	echo -e ' \\textbf{Code}:\n\\begin{lstlisting}' >> tmp.tex
	tail -n+2 "${F}" | grep -v "^#'" >> tmp.tex
	echo -e '\\end{lstlisting}\n\\textbf{Output}:\n\\begin{lstlisting}' >> tmp.tex
	tail -n+2 "${F}" | grep -v "^#'" | R --no-save --quiet  >> tmp.tex
	echo -e '\\end{lstlisting}\n'  >> tmp.tex

done


cat << __EOF__ >> tmp.tex
\end{document}
__EOF__

pdflatex tmp.tex
mv tmp.pdf rbcf_examples.pdf
rm tmp.tex

