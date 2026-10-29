# SNP-panel pipeline — interface reserved, not implemented

Future independent program:

```text
one/multiple species-consistent VCF/BCF + population_labels.tsv
    -> sample/variant QC -> per-population statistics -> SNP ranking
    -> LD filtering -> candidate panel sizes -> classifier
    -> nested cross-validation -> smallest **tested** qualifying panel
```

Inputs: `sample`, `population` TSV columns and one or more VCF/BCF.
`balanced_accuracy >= 0.95` is an example, **not** a default target.
If no evaluated size meets the threshold, return an explicit "not achieved"
status; do not report a non-qualifying panel as successful.

For unbiased estimates, fit data-dependent QC thresholds, FST ranking,
LD filter, panel-size choice and classifier hyperparameters within each
relevant training split. See `docs/data_contracts.md`.

There is intentionally no executable here yet and no import from WGS modules.
