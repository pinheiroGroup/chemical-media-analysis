# chemical-media-analysis

Analysis pipeline and outputs for the chemical-media growth screen reported in
the GUIbiont manuscript. The published dataset contains 13,608 original
*E. coli* BW25113 growth-curve records; the manuscript pipeline retains 13,400
curves with sufficient numeric measurements, spanning 1,026 defined media,
across seven experimental rounds. The exact nine-workbook input snapshot is
maintained in
[`chemical-media-dataset`](https://github.com/pinheiroGroup/chemical-media-dataset)
(seven growth workbooks, the curve-to-medium mapping and the 44-compound medium
composition table). The data from the associated
[Scientific Data article](https://doi.org/10.1038/s41597-025-05356-3)
are converted to GUIbiont experiment folders, fitted with GUIbiont's log-linear
sliding-window estimator, and analysed against a per-medium matrix of 44
compound concentrations to relate medium composition to maximum growth rate and
saturation OD.

## Reproduction

See [`GUIBIONT_REPRODUCTION.md`](GUIBIONT_REPRODUCTION.md) for the full
step-by-step reproduction: prerequisites, the five pipeline commands, the
expected counts at each stage, and the list of canonical outputs. To verify an
existing checkout without rerunning the pipeline:

```bash
python scripts/validate_results.py
```

## License

Repository-authored software is available under the MIT License; see
[`LICENSE`](LICENSE). The underlying growth-curve data originate from the
Scientific Data article cited above and are distributed under CC BY 4.0; see
[`DATA_LICENSE`](DATA_LICENSE).
