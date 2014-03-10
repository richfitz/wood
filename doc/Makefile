all: wood-ms.pdf

wood-ms.pdf: wood-ms.tex wood.bib
	pdflatex wood-ms
	bibtex wood-ms
	pdflatex wood-ms
	pdflatex wood-ms

clean: tidy
	rm -f wood-ms.pdf

tidy:
	rm -f wood-ms.aux wood-ms.bbl wood-ms.blg wood-ms.log

.PHONY: all clean tidy
