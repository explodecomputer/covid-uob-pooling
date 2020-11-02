# Pooling strategies for detecting Covid in University of Bristol students

The University of Bristol UoB plans to commence the 2020 academic year with a strategy of placing students in living circles, to minimise the spread of SARS-CoV-2 across the student body. To monitor outbreaks it would be costly to simply test everybody on a regular basis, and so whether pooling samples could effectively reduce costs and still adequately detect outbreaks is an open question.

## To run

Obtain xlsx file containing living circle counts for halls of residence, place in `docs`, then run

```
snakemake
```

The simulations take about 1 hour on 16 cores.



---

## Notes

Papers on pooling

- https://www.medrxiv.org/content/10.1101/2020.10.02.20204859v1.full.pdf
- https://www.medrxiv.org/content/10.1101/2020.09.02.20183830v1.full.pdf
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3983103/
- https://www.nejm.org/doi/full/10.1056/NEJMp2025631


