#!/bin/sh
# Compile fenics-mixed.do.txt to HTML

# plain HTML
doconce format html fenics-mixed --html-style=bloodish #--pygments-html-linenos
if [ $? -ne 0 ]; then echo 'could not compile'; exit; fi
cp fenics-mixed.html ../pub
cp -r fig-fenics-mixed ../pub/
# view ../pub/fenics-mixed.html in a browser

#exit # drop Sphinx

# Sphinx HTML
doconce format sphinx fenics-mixed
doconce split_rst fenics-mixed
theme=fenics # alternative: cbc
doconce sphinx_dir author="K.-A. Mardal and H. P. Langtangen" version=1.0 theme=$theme fenics-mixed
python automake_sphinx.py
rm -rf ../pub/sphinx
cp -r sphinx-rootdir/_build/html ../pub/sphinx
# Build more themes for the fun of it
cd sphinx-rootdir
sh make_themes.sh cbc pyramid basicstrap
cd ..
cp -r sphinx-rootdir/sphinx-* ../pub/
# Display in Google Chrome browser
#firefox ../pub/sphinx/index.html &

# Content
#doconce format html index0 --html-style=bloodish
#cp index0.html ../web/index.html
#doconce sphinx_dir author="K.-A. Mardal and H. P. Langtangen" version=1.0 theme=fenics index0
#python automake_sphinx.py
#rm -rf ../web/sphinx
#cp -r sphinx-rootdir/_build/html ../web/sphinx
