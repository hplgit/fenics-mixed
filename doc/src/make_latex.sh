#!/bin/sh
# Compile fenics-mixed.do.txt to PDF via pdflatex.
# Make one version for screen reading and one for printing on paper
# (the difference is the way links are handled).

doconce spellcheck -d .dict4spell.txt fenics-mixed.do.txt
if [ $? -ne 0 ]; then
  echo "Abort due to misspellings!"
  exit 1
fi
rm -rf tmp_stripped

# Make sure publish database is up-to-date if refs.bib is the
# file that is really maintained
rm -f publish_references.pub
publish import refs.bib <<EOF
1
1
EOF

# Version 1: use anslistings.sty (from the FEniCS book) to
# typeset code
doconce format pdflatex fenics-mixed --device=screen --skip_inline_comments
if [ $? -ne 0 ]; then echo 'could not compile'; exit; fi
# Fix *local* anslistings.sty to drop line numbers
doconce subst -m '^(numbers=.+)' '%\g<1>' anslistings.sty
# Turn .p.tex to .tex
doconce ptex2tex fenics-mixed envir=ans:nt -DTODONOTES #-DLATEX_HEADING=traditional
# Run ordinary pdflatex
pdflatex fenics-mixed
bibtex fenics-mixed
pdflatex fenics-mixed
pdflatex fenics-mixed
cp fenics-mixed.pdf ../pub

# Make a version for printing (links appear with URLs in footnotes)
doconce format pdflatex fenics-mixed --device=paper  --skip_inline_comments #-DLATEX_HEADING=traditional
doconce ptex2tex fenics-mixed envir=ans:nt -DTODONOTES
pdflatex fenics-mixed
pdflatex fenics-mixed
cp fenics-mixed.pdf ../pub/fenics-mixed-4paper.pdf



