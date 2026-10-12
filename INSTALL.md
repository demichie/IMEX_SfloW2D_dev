# IMEX_SfloW2D documentation site integration

Extract this ZIP **at the repository root**, not in `EXAMPLES/`.

```bash
unzip IMEX_SfloW2D_GitHub_Pages_site.zip -d /path/to/IMEX_SfloW2D_dev
```

The archive contains `mkdocs.yml`, `requirements-docs.txt`, `docs/website/`, `tools/docs/` and `.github/workflows/publish-documentation.yml`. It does not change the Fortran solver, `template.doxy`, or the README files already in `EXAMPLES/`.

## First publish

1. Commit these files to the repository's default branch (`master`).
2. Go to GitHub **Settings > Pages > Build and deployment > Source > GitHub Actions**.
3. Open **Actions > Publish documentation (MkDocs + Doxygen) > Run workflow**, selecting `master`.
4. After successful deployment, open https://demichie.github.io/IMEX_SfloW2D_dev/.

The workflow also triggers on tag push (`*`) and release publication. Both events can occur for the same version; runs are serialized, so the final deployment represents whichever run finishes last. The checkout for a release uses its tag. No normal branch push or pull request runs this workflow.

Note: GitHub may require permissions for Pages deployment in repository settings. Publishing a release referring to an older tag documents that older revision and replaces the currently published Pages website.

## Maintenance

- Add a new `EXAMPLES/**/README.md` with front matter (`title`, `summary`, `category`, `tags`) and the standard content sections; the site generator discovers it recursively.
- Preserve older `README.txt` files until historical audit completes.
- Edit static pages in `docs/website/`.
- Edit `template.doxy` for general Doxygen settings; the workflow overrides the output folder so the result is under `/doxygen/`.
- The build fails for missing or invalid README metadata; missing content sections are listed as warnings in the generated audit file.
- Run locally: `python tools/docs/build_examples.py --repo . --output build/docs-generated --revision master`, then `mkdocs build --strict -f build/docs-generated/mkdocs.generated.yml -d "$PWD/build/site"`.

## Caveats

The site does not itself validate numerical input values or execute examples. The project is configured to display the **currently deployed revision**, not to retain documentation for multiple versions. Markdown links to non-documentation files are rewritten to GitHub links pinned to the source commit.
