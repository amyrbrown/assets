http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1002882

Rockwell G, Guido NJ, Church GM (2013) Redirector: Designing Cell
Factories by Reconstructing the Metabolic Objective. PLoS Comput Biol
9(1): e1002882. doi:10.1371/journal.pcbi.1002882

grockwell@genetics.med.harvard.edu

--------------------------------------------------------------------------------

* The core is built with Python and currently uses the GLPK and SCIP
  solvers. LP optimizations were largely carried out in GLPK because
  of the ability to directly access GLPK functions from Python while
  MILP optimizations are carried out by SCIP for faster solving
  speed. Computation was carried out on the Broad Institute
  computational cluster. The Redirector Package including operational
  software code and metabolic network model files used for this
  publication are available at
  https://github.com//bionomicron/Redirect​or.git.

--------------------------------------------------------------------------------

* Code from https://github.com/bionomicron/Redirector
* File from https://raw.github.com/bionomicron/Redirector/master/core/genetic/ISequenceAnalysis.py
* and https://github.com/bionomicron/Redirector/blob/master/core/model/LPSolver.py

