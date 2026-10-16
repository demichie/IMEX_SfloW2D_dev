# IMEX_SfloW2D

**A depth-averaged numerical model for geophysical and volcanic flows.**

IMEX_SfloW2D combines physical flow models with high-resolution spatial discretizations and IMEX time integration. This website provides the documentation, example configurations and Fortran source reference for the development repository.

[Explore examples](examples/index.md){ .md-button .md-button--primary }
[Getting started](getting-started/index.md){ .md-button }
[Source code](https://github.com/demichie/IMEX_SfloW2D_dev){ .md-button }

## Documentation

- **[Getting started](getting-started/index.md)** — build prerequisites and first simulation.
- **[Model](model/index.md)** — physical and numerical formulation overview.
- **[Examples](examples/index.md)** — automatically discovered simulation cases, with input parameters and execution instructions.
- **[Development](development/index.md)** — testing and development guidance.
- **[Doxygen](https://demichie.github.io/IMEX_SfloW2D_dev/doxygen/)** — generated Fortran source-code documentation.
- **[References](references/index.md)** — publication and repository links.

!!! info "Automatically published releases"
    The documentation is built on demand or when a tag or release is published. Example pages are regenerated from the README.md files in `EXAMPLES/` at the selected revision.

!!! warning "GitHub-generated source snapshots"
    GitHub's automatic **Source code (zip/tar.gz)** downloads are snapshots of
    the tagged repository, not the prepared Autotools source archive; they may
    omit generated files such as `configure`. When available, use the named
    `IMEX_SfloW2D-<version>.tar.gz` and
    `IMEX_SfloW2D-<version>-examples.tar.gz` assets on the
    [release page](https://github.com/demichie/IMEX_SfloW2D_dev/releases).
