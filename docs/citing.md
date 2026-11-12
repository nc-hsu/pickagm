# Citing pickagm

If you use pickagm in published work, please cite the software and give
the exact version you used, for example "pickagm 0.1.0".

The citation metadata is in
[`CITATION.cff`](https://github.com/nc-hsu/pickagm/blob/main/CITATION.cff)
in the repository; on GitHub, the "Cite this repository" button turns it
into APA or BibTeX format.

```text
Clemett, N. pickagm (version 0.1.0) [Computer software].
https://github.com/nc-hsu/pickagm
```

## Version 0.1.0 and the paper's record selection

The code in `src/pickagm/` of version 0.1.0 is identical to commit
[`9175fe6`](https://github.com/nc-hsu/pickagm/commit/9175fe6c825feb42087765537bc4b377a1ef5f43)
(tag `record-selection-run`), the commit that the record selection for the
accompanying paper ran on. To reproduce that selection, install
`pickagm==0.1.0`.

## Data and models

If your work depends on the bundled data, please also cite its sources
(see {doc}`data`):

- the ESM flatfile: Lanzano et al. (2018),
  <https://doi.org/10.13127/esm/flatfile.1.0>, and Lanzano et al. (2019),
  <https://doi.org/10.1007/s10518-018-0480-z>;
- the correlation models: Clemett and Gündel (paper in preparation).

Please also cite the ground-motion models you use, as given in their
OpenQuake documentation.
