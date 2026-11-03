#!/usr/bin/env bash
set -e

MAIN="leaf_gbc1_npj_microgravity"
ROOT_DIR=".."

echo "=========================================================="
echo "Compiling LEAF GBC1 Manuscript for npj Microgravity"
echo "=========================================================="

echo ">>> [1/4] Running Initial pdflatex..."
pdflatex -interaction=nonstopmode "${MAIN}.tex"

echo ">>> [2/4] Running BibTeX..."
bibtex "${MAIN}"

echo ">>> [3/4] Running Second pdflatex (resolving citations)..."
pdflatex -interaction=nonstopmode "${MAIN}.tex"

echo ">>> [4/4] Running Final pdflatex (resolving references & labels)..."
pdflatex -interaction=nonstopmode "${MAIN}.tex"

if [ -f "${MAIN}.pdf" ]; then
    echo "=========================================================="
    echo "SUCCESS: ${MAIN}.pdf generated successfully!"
    echo "File size: $(ls -lh ${MAIN}.pdf | awk '{print $5}')"
    echo "=========================================================="

    # Copy to root directory for direct GitHub Pages URL access
    cp -f "${MAIN}.pdf" "${ROOT_DIR}/${MAIN}.pdf"
    echo "Copied ${MAIN}.pdf to root directory for GitHub Pages hosting."

    # Build ZIP distribution archive
    echo "Creating Overleaf / Submission ZIP package..."
    rm -f "${ROOT_DIR}/leaf_gbc1_latex_package.zip"
    zip -r "${ROOT_DIR}/leaf_gbc1_latex_package.zip" \
        sn-jnl.cls sn-nature.bst sn-basic.bst sn-mathphys.bst \
        "${MAIN}.tex" main.tex references.bib sn-bibliography.bib \
        Makefile compile.sh README.md \
        figures/ \
        "${MAIN}.pdf"
    echo "Created ${ROOT_DIR}/leaf_gbc1_latex_package.zip"
else
    echo "ERROR: Compilation failed to generate ${MAIN}.pdf."
    exit 1
fi
