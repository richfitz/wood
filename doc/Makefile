LATEX_EXTENSIONS={aux,bbl,blg,lof,log,out}

all: wood-ms.pdf wood-ms-supporting.pdf

# Main manuscript depends on figure references in supplement:
wood-ms.pdf: wood-ms-supporting.pdf

%.pdf: %.tex wood.bib
	pdflatex -interaction=nonstopmode $<
	bibtex $(basename $<)
	pdflatex -interaction=nonstopmode $<
	pdflatex -interaction=nonstopmode $<
	pdflatex -interaction=nonstopmode $<

clean: tidy
	rm -f wood-ms.pdf wood-ms-supporting.pdf

tidy:
	rm -f wood-ms.${LATEX_EXTENSIONS} wood-ms-supporting.${LATEX_EXTENSIONS}

.PHONY: all clean tidy
