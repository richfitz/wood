all: wood-ms.pdf

wood-ms.pdf: wood-ms.tex wood.bib
	pdflatex -interaction=nonstopmode wood-ms
	bibtex wood-ms
	pdflatex -interaction=nonstopmode wood-ms
	pdflatex -interaction=nonstopmode wood-ms
	pdflatex -interaction=nonstopmode wood-ms

clean: tidy
	rm -f wood-ms.pdf

tidy:
	rm -f wood-ms.aux wood-ms.bbl wood-ms.blg wood-ms.log

.PHONY: all clean tidy
