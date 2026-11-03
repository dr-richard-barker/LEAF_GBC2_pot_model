# Makefile for LEAF GBC1 npj Microgravity LaTeX Manuscript
MAIN = leaf_gbc1_npj_microgravity
PDF = $(MAIN).pdf
ZIP = leaf_gbc1_latex_package.zip
ROOT_DIR = ..

.PHONY: all pdf package clean distclean deploy

all: pdf package deploy

pdf: $(PDF)

$(PDF): $(MAIN).tex references.bib sn-jnl.cls sn-nature.bst
	pdflatex -interaction=nonstopmode $(MAIN).tex
	bibtex $(MAIN)
	pdflatex -interaction=nonstopmode $(MAIN).tex
	pdflatex -interaction=nonstopmode $(MAIN).tex

package: $(PDF)
	@echo "Creating self-contained LaTeX distribution package..."
	rm -f $(ROOT_DIR)/$(ZIP)
	zip -r $(ROOT_DIR)/$(ZIP) \
		sn-jnl.cls sn-nature.bst sn-basic.bst sn-mathphys.bst \
		$(MAIN).tex main.tex references.bib sn-bibliography.bib \
		Makefile compile.sh README.md \
		figures/ \
		$(PDF)
	@echo "Created $(ROOT_DIR)/$(ZIP)"

deploy: $(PDF) package
	@echo "Copying production PDF to web portal root..."
	cp -f $(PDF) $(ROOT_DIR)/$(PDF)
	@echo "Deployment assets refreshed."

clean:
	rm -f *.aux *.bbl *.blg *.log *.out *.toc *.synctex.gz *.fls *.fdb_latexmk *.upa *.upb

distclean: clean
	rm -f $(PDF)
