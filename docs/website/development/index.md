# Development and validation

IMEX-SfloW2D is maintained as a scientific Fortran solver with explicit ownership of canonical state and physical configurations. The integrated single- and two-layer paths are tested together.

## Current workflows

The active GitHub Actions workflows are **manual** (`workflow_dispatch`). Ordinary commits, PRs and tags do not automatically compile, validate, publish documentation or release a container.

- [Validation and acceptance](validation.md): optimized/strict builds, physical and numerical regressions, solver/OpenMP integration, and scientific acceptance scope.
- [Documentation and Docker publication](publishing.md): generate MkDocs and Doxygen, deploy Pages and publish a GHCR container intentionally.
- [Contributing guidance](https://github.com/demichie/IMEX_SfloW2D_dev/blob/master/CONTRIBUTING.md)
- [Current development roadmap](https://github.com/demichie/IMEX_SfloW2D_dev/blob/master/docs/IMEX_SfloW2D_Development_Roadmap.md)

The GitHub [Actions dashboard](https://github.com/demichie/IMEX_SfloW2D_dev/actions) records results for specific commit SHAs; a green run on another commit is not evidence that the current branch has been validated.

Source documentation is generated with [Doxygen](https://demichie.github.io/IMEX_SfloW2D_dev/doxygen/) from the exact revision selected for Pages deployment.
