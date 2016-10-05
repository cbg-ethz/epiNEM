# epiNEM

EpiNEM is a direct extension of the NEM framework. epiNEM focuses on the epistasis between three S-genes, namely two parents and one child. Epistasis occurs, when the single effects of the two parents do not explain their combined effect on the child. EpiNEM models this epistasis with five Boolean logical gates for a triplet A->C, B->C.

- A or B, OR: A has an effect on C, B has an effect on C, A+B has an effect on C
- A and B, AND: A has no effect on C, B has no effect on C, A+B has an effect on C
- (not A and B) or (A and not B), XOR: A has an effect on C, B has an effect on C, A+B has no effect on C
- not A and B, A masks the effect of B: A has no effect on C, B has an effect on C, A+B has no effect on C
- A and not B, B masks the effect of A: A has an effect on C, B has no effect on C, A+B has no effect on C

Install:
--------

Open R and input:

```r
install.packages("devtools")

library(devtools)

install_github("cbg-ethz/epiNEM")

library(epiNEM)
```

Then check out the vignette for working examples.

```r
vignette(package="epiNEM")
```
