# archive

Nothing here is part of the pipeline. Kept for reference only.

**Do not use the tier scripts.** `classify_confidence.py` and
`tier_assign_202609.py` are earlier re-implementations of the tier rules in
`gene_prio_utils.R`. Checked against the published assignment they disagree with
it on 7.5% and 10.3% of genes respectively, so running them will not reproduce
the released tiers. The authoritative assignment is made by the build itself and
carried in the `top_confidence` column of the release.

`tier_assign_202609.py` also only prints counts; it writes nothing.

`dev/` holds ad-hoc scripts used to inspect intermediate tables during
development. `jobs/` holds cluster submission scripts specific to the BU SCC.
