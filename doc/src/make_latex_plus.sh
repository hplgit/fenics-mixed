#!/bin/sh
# Basically as make_latex.sh, but make more versions with different
# styles for computer code.

# Version 1: use anslistings.sty (from the FEniCS book) to
# typeset code
doconce format pdflatex fenics-mixed
# Fix anslistings.sty to drop line numbers
doconce subst -m '^(numbers=.+)' '%\g<1>' anslistings.sty
# Turn .p.tex to .tex
#doconce ptex2tex fenics-mixed
#doconce ptex2tex fenics-mixed cpppro=ans:nt cppcod=ans:nt pypro=ans:nt pycod=ans:nt cpro=ans:nt ccod=ans:nt fcod=ans:nt fpro=ans:nt shpro=ans:nt shcod=ans:nt
doconce ptex2tex fenics-mixed envir=ans:nt
pdflatex fenics-mixed
bibtex fenics-mixed
pdflatex fenics-mixed
pdflatex fenics-mixed
cp fenics-mixed.pdf doc

# Version 2: use minted.sty (and the Pygments program) to
# typeset code
doconce ptex2tex fenics-mixed envir=minted
# Add line numbers by fixing \begin{minted}[...,linenos=true,...] headings
doconce replace "linenos=false" "linenos=true" fenics-mixed.tex
pdflatex -shell-escape fenics-mixed
cp fenics-mixed.pdf doc/fenics-mixed-minted.pdf

# Version 3: use Verbatim (fancyvrb) to typeset code
# (problem: get numbered lines also for terminal sessions,
# fixed in version 4 below)
doconce ptex2tex fenics-mixed
# All code environments are now standard \begin{verbatim}
doconce replace ',amsfonts}' ',amsfonts,fancyvrb}' fenics-mixed.tex
# Add line numbers by fixing \begin{minted}[...,linenos=true,...] headings
doconce replace "\begin{quote}\begin{verbatim}" "\begin{Verbatim}[numbers=left,fontsize=\fontsize{9pt}{9pt},tabsize=8,baselinestretch=1.0]" fenics-mixed.tex
doconce replace "\end{verbatim}\end{quote}" "\end{Verbatim}" fenics-mixed.tex
pdflatex fenics-mixed
cp fenics-mixed.pdf doc/fenics-mixed-fancyvrb0.pdf

# Version 4: let sys and sh environments be a fancy Verbatim environment with
# title and horizontal rule and the rest of all environments be Verbatim
# with numbered lines
terminal="\begin{Verbatim}[numbers=none,frame=lines,label=\fbox{{\tiny Terminal}},fontsize=\fontsize{9pt}{9pt},labelposition=topline,framesep=2.5mm,framerule=0.7pt]@\end{Verbatim}"
doconce ptex2tex fenics-mixed "sys=$terminal" "shcod=$terminal" "shpro=$terminal" "envir=\begin{Verbatim}[numbers=left]@\end{Verbatim}"
pdflatex fenics-mixed
cp fenics-mixed.pdf doc/fenics-mixed-fancyvrb.pdf

#doconce clean  # puts temporary files to be cleaned in Trash dir
