# Symmetric simulation: original estimates and regenerated summaries

The parameter CSV contains all 15,000 original paired estimates: 5,000 at
each of n = 500, 1000 and 2000. Its original MD5 matches the recorded input
manifest. Every pair is recorded as successful.

The supplied original raw-risk export contained 1,260,000 rows, with 84
cost ratios per pair. Every risk value was checked against its corresponding
fitted parameters. The compact risk summary was aggregated from those
original rows and reproduces the manuscript's Table 1, equal-cost intervals
and RMSE values.

The 188,751,829-byte raw-risk CSV is omitted because the fitted parameters
and original cost grid are sufficient to regenerate it without sampling
or fitting. From the repository root:

```bash
python 09_rebuild_symmetric_risks.py .
python 09_rebuild_symmetric_risks.py . --write-raw
```

The first command verifies/rebuilds the summary; the second additionally
writes the large raw export. The original raw export's MD5 is recorded in
`validation/symmetric_validation.json`. Numerical equivalence does not
imply byte-identical CSV formatting.
