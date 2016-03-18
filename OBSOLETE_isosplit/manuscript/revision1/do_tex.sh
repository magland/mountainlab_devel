#!/bin/bash
latex isosplit-revision1.tex

if [ $? -neq 0 ]
then
	exit 0
fi

bibtex isosplit-revision1
latex isosplit-revision1.tex
pdflatex isosplit-revision1.tex
#dvipdfm isosplit-revision1.dvi

