# NomGen
## Overview:
A burgeoning R package to generate correlated nominal data.

Based on an algorithm by A.J. Lee (1997) [1].

Find above 3 completed implementations (so far):
- CatGen2Var generates 2D categorical vectors with specificed correlation.
- 3VarNominal generated 3D categorical vectors with correlation specified between 2 pairs of variables.
- 3VarNom3Corre does the same as 3VarNominal but accepts a complete correlation structure, returning an error if it is not coherent.

Along with code to generate data presented at 3 minute UROP presentation 12/10/22.

Ultimately, an R package will be compiled which can generate multivariate categorical data of any dimension.

## Notes on usage:
- Owing to inefficiencies, exercise caution when applying algorithm for large numbers of categories (>5), as this can use an extoritionate amount of memory. 

## Sources

[1] Lee, A. , “Some simple methods for generating correlated categorical variates”, Computational Statistics & Data Analysis, vol. 26, 1997
