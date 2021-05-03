# epiNEM

EpiNEM is a direct extension of the NEM framework. epiNEM focuses on the epistasis between three S-genes, namely two parents and one child. Epistasis occurs, when the single effects of the two parents do not explain their combined effect on the child. EpiNEM models this epistasis with five Boolean logical gates for a triplet A->C, B->C.

- A or B, OR: A has an effect on C, B has an effect on C, A+B has an effect on C
- A and B, AND: A has no effect on C, B has no effect on C, A+B has an effect on C
- (not A and B) or (A and not B), XOR: A has an effect on C, B has an effect on C, A+B has no effect on C
- not A and B, A masks the effect of B: A has no effect on C, B has an effect on C, A+B has no effect on C
- A and not B, B masks the effect of A: A has an effect on C, B has no effect on C, A+B has no effect on C

Install:
--------

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("epiNEM")
```

Most recent (devel) version:

```r
install.packages("devtools")

library(devtools)

install_github("cbg-ethz/epiNEM")

library(epiNEM)
```

Check out the vignette for working examples and the reproduction of the publication.

```r
vignette(package="epiNEM")
```

Toy example:

```r
data <- matrix(sample(c(0,1), 100*4, replace = T), 100, 4)

colnames(data) <- c("A", "A.B", "B", "C")

rownames(data) <- paste("E", 1:100, sep = "_")

res <- epiNEM(data, method = "exhaustive")

plot(res)
```
## References

Pirkl M, Diekmann M, van der Wees M, Beerenwinkel N, Fröhlich H, Markowetz F (2017). "Inferring modulators of genetic interactions with epistatic nested effects models." PLOS Computational Biology, 13, 1-18. doi: 10.1371/journal.pcbi.1005496.
