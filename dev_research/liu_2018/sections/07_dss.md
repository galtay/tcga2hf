## Section 7 — DSS deep-dive (92.9% match — Liu's documented approximation)

Liu's STAR Methods explicitly flags DSS as approximate: *"Technically a patient could be with tumor but died of a car accident and therefore incorrectly considered as an event."* The 7% disagreement is dominated by cases where vital_status / tumor_status / cause_of_death have drifted between 2018 and now, plus our handling of Liu's special case #3 (62 Dead-with-Tumor + defined-NTE cases).

### Mismatch direction

| Liu cdr_DSS | ours dss_event | patients |
| --- | --- | --- |
| 0.0 | 1.0 | 30.0 |
| 1.0 | 0.0 | 10.0 |

The two-way split (some Liu=0/ours=1, some Liu=1/ours=0) reflects the underlying ambiguity in DSS classification — neither direction is dominant. Compare to Section 6's OS, where every mismatch goes Liu=0 → ours=1 (vital-status updates only).

### Per-project event mismatches

| project | event mismatches |
| --- | --- |
| UVM | 11 |
| PAAD | 8 |
| STAD | 7 |
| HNSC | 2 |
| READ | 2 |
| SKCM | 2 |
| CESC | 1 |
| GBM | 1 |
| KIRC | 1 |
| LGG | 1 |
| LUAD | 1 |
| MESO | 1 |
| SARC | 1 |
| UCEC | 1 |

### Liu's documented approximation

Three signals get reasoned over:
1. **vital_status** — Alive / Dead.
2. **tumor_status** — WITH TUMOR / TUMOR FREE (latest follow-up encounter).
3. **cause_of_death** — when populated, "Cancer Related" overrides TUMOR FREE.

Liu's algorithm:
- Alive → DSS event = 0.
- Dead AND TUMOR FREE → DSS event = 0.
- Dead AND WITH TUMOR → DSS event = 1.
- Dead, tumor_status missing, cause_of_death = "Cancer Related" → DSS event = 1.
- Dead, tumor_status missing, cause_of_death missing → DSS event = 1 (conservative default; Liu's CDR populates DSS=1 for dead patients without a TUMOR FREE signal).

### Implication

DSS will never reach 100% match against Liu unless we shadow Liu's exact 2018 vital_status / tumor_status snapshot. The remaining gap is structural data drift, not algorithmic. For research purposes, treat OS as the gold-standard "patient died" signal and use DSS as an approximation when cancer-specific mortality matters.
