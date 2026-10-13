Model documentation (leaf-conductance.pdf)
==========================================

Sources
  leaf-conductance.tex      LaTeX source
  extra.bib                 sources cited only here (bibliography otherwise ../references.bib)
  make_figures.py           figures (fig/) and numbers (generated/values.tex) from the model code
  matlab_crosscheck.m       MATLAB vs. reference values -> generated/matlab_crosscheck.json

Build (from the repository root)
  matlab -batch "cd docs/manual; matlab_crosscheck('generated/matlab_crosscheck.json')"
  python docs/manual/make_figures.py
  cd docs/manual
  pdflatex leaf-conductance && bibtex leaf-conductance && pdflatex leaf-conductance && pdflatex leaf-conductance

The MATLAB step is optional; without it the MATLAB columns of the cross-check table are missing
and the document does not compile until the macros \MlDev* / \MlVals* / \MatlabVersion are defined.
