# LEAF GBC1 npj Microgravity LaTeX Manuscript Package

This directory contains the publication-ready LaTeX source package for the article:

> **Computational fluid dynamics, multi-species canopy gas exchange, and harvest fixation architecture of the LEAF GBC1 lunar plant growth chamber**  
> *Richard Barker (Department of Botany and Plant Pathology, Purdue University, barkerrj@purdue.edu)*  
> Target Journal: **npj Microgravity** (Nature Portfolio / Springer Nature)

---

## Directory Contents

- `leaf_gbc1_npj_microgravity.tex` — Main LaTeX manuscript source formatted in modern double-column *npj Microgravity* style (`iicol`).
- `main.tex` — Symlink to `leaf_gbc1_npj_microgravity.tex` for automatic Overleaf detection.
- `references.bib` / `sn-bibliography.bib` — Full BibTeX database of cited models and spaceflight literature.
- `sn-jnl.cls` — Official Springer Nature unified journal document class (`[pdflatex,sn-nature,iicol]`).
- `sn-nature.bst` — Official Nature numerical superscript citation bibliography style.
- `sn-basic.bst`, `sn-mathphys.bst` — Auxiliary Springer Nature bibliography styles.
- `figures/` — Directory containing all 25 high-resolution figures, blueprints, and 4D mission posters.
- `Makefile` — Build automation for PDF compilation, packaging, and deployment.
- `compile.sh` — Bash compilation script executing `pdflatex` + `bibtex` + `pdflatex` + `pdflatex`.
- `leaf_gbc1_npj_microgravity.pdf` — Compiled production PDF.

---

## Local Compilation

To compile locally using `pdflatex` and `bibtex`:

```bash
cd latex
./compile.sh
# Or using make:
make
```

### Clean Build
```bash
make clean
```

---

## Overleaf Import Instructions

1. Download the bundled ZIP package from the repository root: `leaf_gbc1_latex_package.zip` (or download directly from the GitHub Pages portal).
2. Go to [Overleaf](https://www.overleaf.com/) and click **New Project** $\to$ **Upload Project**.
3. Select `leaf_gbc1_latex_package.zip`.
4. Overleaf automatically detects `main.tex` and compiles using **pdfLaTeX**.
