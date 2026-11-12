# pickagm

:::{warning}
**pickagm is under heavy development and may change at any time.**
Functions, arguments and results can change between releases without
notice. If you need exactly the same results later, pin the version, for
example `pickagm==0.1.0`.
:::

pickagm is a Python library for selecting earthquake ground-motion records
for nonlinear structural analysis.

When engineers test a structural model under earthquake shaking, they run it
with a set of recorded ground motions. The set should represent the shaking
that the seismic hazard at the site predicts: not just one "average"
spectrum, but the spread of spectral accelerations, durations and other
properties that the earthquakes behind the hazard would produce. A set like
that is called *hazard-consistent*.

pickagm selects such sets with the **generalised conditional intensity
measure (GCIM)** approach. You can think of it like casting actors for a
crowd scene: first you describe the crowd you need (the target
distributions), then you write a short description of each person in it
(random realisations of the targets), and finally you go through the
records you have and pick the closest match for each description.

pickagm builds on [OpenQuake](https://github.com/gem/oq-engine): it uses
OpenQuake's intensity measure types and ground-motion models, and its
rupture contexts are NumPy record arrays in the form OpenQuake models
expect. The [ESM flatfile](https://doi.org/10.13127/esm/flatfile.1.0) is
included as a record database.

## Where to start

- {doc}`installation`: how to install pickagm and OpenQuake.
- {doc}`workflow`: the steps of a selection and the functions for each.
- {doc}`example`: a complete, small selection you can run.
- {doc}`api/index`: every module, function and class.

```{toctree}
:maxdepth: 2
:caption: User guide
:hidden:

installation
workflow
conventions
correlation-models
data
limitations
example
```

```{toctree}
:maxdepth: 2
:caption: Reference
:hidden:

api/index
```

```{toctree}
:maxdepth: 1
:caption: About
:hidden:

citing
changelog
licence
```
