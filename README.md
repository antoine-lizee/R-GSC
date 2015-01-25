# R-GSC

An *R* implementation of the *Gerstein-Sonnhammer-Chothia algorithm* described in:

> Gerstein M, Sonnhammer ELL, Chothian C (1994) Volume Changes in Protein Evolution. *Journal of Molecular Biology*. [doi:10.1016/0022-2836(94)90012-4](http://dx.doi.org/10.1016/0022-2836(94)90012-4)

The algorithm is described in the paper's appendix titled, "A Method to Weight Protein Sequences to Correct for Unequal Representation." [[Download PDF](https://pdf.yt/d/Sx3jMbr8vANgxAej/download)].

The algorithm weights the leaves of a [dendrogram](https://en.wikipedia.org/wiki/Dendrogram) based on their underrepresentation.

Run the algorithm with the `GSC()` function from the `GSC.R` file. The single input to `GSC()` is an [R dendrogram object](https://stat.ethz.ch/R-manual/R-patched/library/stats/html/dendrogram.html).
