http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1002957

Park CY, Wong AK, Greene CS, Rowland J, Guan Y, et al. (2013)
Functional Knowledge Transfer for High-accuracy Prediction of
Under-studied Biological Processes. PLoS Comput Biol 9(3):
e1002957. doi:10.1371/journal.pcbi.1002957

ogt@genomics.princeton.edu

--------------------------------------------------------------------------------

* All software used in this study has been implemented in the open
  source and publicly available Sleipnir library [94] available from
  http://libsleipnir.bitbucket.org, which interfaces with the SVMperf
  library [95] for linear kernel SVM classifiers (the error parameter
  C was set to 100 for these experiments through
  cross-validation). L1-regularized logistic regression used the
  LIBLINEAR [28] and Random forest used the MILK (Machine Learning
  Toolkit) python package implementation with 61 decision trees per GO
  term.
* A comprehensive application of FKT to 11,000 biological processes,
  along with the functional relationship networks for all six
  organisms, are available through the IMP web-server portal
  accessible at http://imp.princeton.edu [26].
* Tried to download http://libsleipnir.bitbucket.org/sleipnir-current.tar.gz
  from http://libsleipnir.bitbucket.org: 404
* Resorted to their Bitbucket repository

--------------------------------------------------------------------------------

* Downloaded https://bitbucket.org/libsleipnir/sleipnir/raw/16f7f785d415d136fcff9100c716b6b082bb160e/src/annotationobo.cpp
  (long) and https://bitbucket.org/libsleipnir/sleipnir/raw/16f7f785d415d136fcff9100c716b6b082bb160e/src/clustpivot.cpp (short)
