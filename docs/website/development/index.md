# Development and validation

The development repository maintains automated solver validation, Intel `ifx` portability and Docker workflows. See [GitHub Actions](https://github.com/demichie/IMEX_SfloW2D_dev/actions) for the current results.

The [N9 scientific baseline and M12c build acceptance records](https://github.com/demichie/IMEX_SfloW2D_dev/tree/master/docs/acceptance) document tested revisions, compiler configurations and limitations. Those acceptance records must not be interpreted as automatic approval of subsequent numerical changes.

## Documentation workflow

The site is rebuilt manually or on tag/release events. Each example's `README.md` is authoritative: its metadata drives the catalogue, and the complete Markdown body is published as the example page. Source-code API pages are created from `template.doxy` in the same build.
