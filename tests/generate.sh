#/bin/bash

set -e

rm -f tmp.aux  tmp.log  tmp.out  tmp.toc

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
  keepspaces=true
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

awk 'BEGIN{P=0;} /^##[ ]*Examples/ {P=1;print;} {if(P==0) print;}' ../README.md > tmp.md

ls *.R | sort | while read F
do
	echo -n '\addcontentsline{toc}{toc}{' >> tmp.tex
	head -n1 "${F}" | cut -c 3- | sed 's/$/}/' >>  tmp.tex
	echo -n '\subsection{' >> tmp.tex
	head -n1 "${F}" | cut -c 3- | sed 's/$/}/' >>  tmp.tex
	tail -n+2 "${F}" | grep  "^#'" | cut -c 3-  >> tmp.tex
	echo -e ' \\textbf{Code}:\n\\begin{lstlisting}' >> tmp.tex
	tail -n+2 "${F}" | grep -v "^#'" >> tmp.tex
	echo -e '\\end{lstlisting}\n\\textbf{Output}:\n\\begin{lstlisting}' >> tmp.tex
	tail -n+2 "${F}" | grep -v "^#'" | R --no-save --slave --quiet |\
		 awk '{if(length($$0)<100) {print} else {printf("%s(...)\n",substr($$0,1,95));}}'  >> tmp.tex
	echo -e '\\end{lstlisting}\n'  >> tmp.tex
	
	echo -n '### ' >> tmp.md
	head -n1 "${F}" | cut -c 3- >>  tmp.md
	tail -n+2 "${F}" | grep  "^#'" | cut -c 3-  >> tmp.md
	echo -e '\n**Code**:\n\n```' >> tmp.md
	tail -n+2 "${F}" | grep -v "^#'" >> tmp.md
	echo -e '```\n\n**Output**:\n\n```' >> tmp.md
	tail -n+2 "${F}" | grep -v "^#'" | R --no-save --slave --quiet  |\
		 awk '{if(length($$0)<100) {print} else {printf("%s(...)\n",substr($$0,1,95));}}'  >> tmp.md
	echo -e '```\n\n'  >> tmp.md

done


cat << __EOF__ >> tmp.tex
\end{document}
__EOF__

pdflatex tmp.tex
mv -v tmp.pdf ../doc/rbcf_examples.pdf
rm -v tmp.tex

mv -v tmp.md  ../README.md

rm -v -f tmp.log tmp.out tmp.aux
