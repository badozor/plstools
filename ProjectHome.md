## Description: ##
This package contains functions to fit PLS and GPLS models and additional diagnostic elements such as specific summaries, graphical representations, etc ... These functions are partly based on the functions proposed in the packages 'pls' and 'gpls'. However, the functions 'plsone' and 'gplsone' enable the conservation of the missing values (NA) in the explanatory variables in using NIPALS algorithm.

**NOTA:** This version is a first prototype (in construction)!
The function 'gplsone' is always in construction and some thoretical assumptions could be debatable (?). In addition, the computational time related to this one can be relatively high. For the future version, external C functions are in development to increase the computational performance. If you have any comments, please don't hesitate to contact me.

## Unexhaustive References: ##
  * Ding, B.Y. and Gentleman, R. (2003) Classification using generalized partial least squares.
  * Marx, B.D (1996) Iteratively reweighted partial least squares estimation for generalized linear regression. Technometrics 38(4): 374-381.
  * Tenenhaus M.(1998) La Regression PLS. Theorie et pratique. Technip, Paris.
  * Wold H. (1966) Estimation of principal components and related models by iterative least squares. In P. Krishnaiah, editors.Multivariate Analysis, Academic Press, 391-420.
  * Wold H. (1975) Modeling in Complex Situations with Soft Information, Third World Congress of Econometric Society, August 21-26, Toronto, Canada.
  * Wold S., Martens H. & Wold H. (1983) The multivariate calibration problem in chemistry solved by the PLS methods, in Proc. Conf. Matrix Pencils, Ruhe, A. & Kågstrøm, B. (Eds), March 1982, Lecture Notes in Mathematics, Springer Verlag, Heidelberg, pp. 286- 293.
  * R Development Core Team (2008) R: A language and environment for statistical computing. In R Foundation for Statistical Computing. Vienna, Austria.

## Depends: ##
R (>= 2.9.0), ade4
## Suggests: ##
boot, pls, gpls
## License: ##
GPL version 2 or newer
## Date: ##
2009/06/22
## Version: ##
alpha 1.0-6 (in construction)