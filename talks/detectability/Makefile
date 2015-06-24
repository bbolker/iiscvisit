%.html: %.rmd
	Rscript -e "library(\"rmarkdown\"); render(\"$<\")"

%.md: %.rmd
	Rscript -e "library(\"knitr\"); knit(\"$<\")"

%.tex: %.md
	pandoc -s -S -t latex -V documentclass=tufte-handout $*.md -o $*.tex

%.pdf: %.tex
	pdflatex --interaction=nonstopmode $*

%.pdf: %.Rnw
	Rscript -e "library(\"knitr\"); knit2pdf(\"$<\")"

clean:
	rm -f *.log *.aux *.md *.out texput.log *.nav *.snm *.toc

