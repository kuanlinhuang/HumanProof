# The grid of experiments, results and features of HumanProof

This page says what each model of the gene safety pipeline asks, how many genes stand behind each number, what each cut of those genes shows, and which inputs carry the signal.
Every number on it is computed from files of this repo by `scripts/grid_doc.py`, with no model retrained, no network and no paid service.
The numbers are also written to `docs/grid/numbers.csv`, and `tests/test_grid_doc.py` fails when this page or that file is stale.
A row of every table is a gene (an Ensembl id of the Open Targets file) and never a drug, a trial or a patient.

```
uv run --with pandas --with numpy --with pyarrow --with scikit-learn python scripts/grid_doc.py
uv run --with pandas --with numpy --with pyarrow --with scikit-learn python scripts/grid_doc.py --check
uv run --with pandas --with numpy --with pyarrow --with scikit-learn --with pytest python -m pytest tests/test_grid_doc.py -q
```

The first command writes the page and the CSV, the second fails when either is stale, and the third runs those checks and the invariants of the counts.

## At a glance

### Summary

The safety model asks whether Open Targets lists a known safety event for a gene, learning from its Open Targets columns, loss of function genetics, constraint and cell type expression, with a drugged-only classifier and a DR+PU pipeline that scores all 17,745 genes.

### Unit

One row is one gene, an Ensembl id of the Open Targets file (never a drug, a trial or a patient), and every file read here holds one row per id with no duplicate id and no duplicate symbol.

### Label

A positive is a gene whose `hasSafetyEvent` is `-1` (945 of the 78,725 genes of the file, 415 of them drugged), an assumed negative is a drugged gene without that flag (1,149), and a gene with neither a phase nor a flag is unlabeled (76,631) and carries no label, so no gene is confirmed free of a safety event.

### Population

The drugged-only model reads the 1,564 genes that have a maximum phase (415 flagged), and the DR+PU pipeline reads 17,745 genes (1,506 drugged, 409 of them flagged) and fits its outcome model on the 1,506 drugged ones.

### Test

Stratified 5-fold cross-validation with shuffling and seed 42 scores each gene once out of fold (XGBoost, early stopping on the scored fold), with 95% percentile intervals from 2,000 bootstrap draws over genes at seed 42, no split in time and no test set held back.

### Metric

AUROC and average precision (AUPRC), each with a 95% bootstrap interval over genes, beside the number of genes, the number of flagged genes and the base rate, which is the AUPRC of a random ranking.

### Baseline

The repo compares each model only with chance (AUROC 0.5, AUPRC equal to the base rate) and fits no baseline model, and for context this page reads each of the 14 Open Targets columns alone, the best being `hasSmallMoleculeBinder` at AUROC 0.636 [0.617, 0.655] in its best direction, against 0.717 for the drugged-only model.

### Not comparable

These are out-of-fold AUROCs against assumed negatives at one snapshot with random folds, so they cannot be compared with the forecast of a later event, with a score on confirmed outcomes, with a deployed score that was fitted to its own labels, or with a test that groups paralogs or splits in time.

### How many

The drugged-only model scores 1,564 genes (415 flagged, base rate 0.265) and the DR+PU pipeline 17,745 genes (1,506 drugged, 409 of them flagged), and the files name no drug and no trial, so drug and trial counts are not in this repo's files.

### Cuts

Every file holds one row per gene, cancer only is not run because no file holds an indication (95 of the 1,564 drugged genes carry a cancer driver flag), combination does not apply to genes, and the stage cut by maximum phase holds 28, 296, 222 and 1,018 drugged genes at the codes 0.25, 0.5, 0.75 and 1.0.

### Features that matter

In the drugged-only model the 14 Open Targets columns carry 95.4% of the mean absolute SHAP (mouseKOScore 0.7683, geneticConstraint 0.5274 and mouseOrthologMaxIdentityPercentage 0.5171 lead), probably because the other 57 columns are filled for few genes, and the DR+PU importance is saved only as a picture.

### Language model only

Not run: the labels are Open Targets 25.12 safety flags, which carry no date and were all published before the 2026-04-24 release of DeepSeek V4 Flash, and well known safety targets are public (8 of 8 listed ones carry the flag), so a language model would be graded on recall and not on a forecast.

## Status of the evidence

All results are out-of-fold numbers on the labels of one Open Targets snapshot, release 25.12.
No result is blind: the snapshot carries no dates, so there is no split in time, and no test set was held back.
The intervals are 95% percentile intervals of 2,000 bootstrap draws over genes at seed 42, so they do not count the choice of settings or a repeat of the folds.
XGBoost stops adding trees when the AUPRC of the scored fold stops improving (line 542 of `data/safety_model.py` and line 224 of `data/safety_model_dr.py`), so the tree count is picked on the genes that are then scored.
The out-of-fold numbers are therefore a little optimistic, by an amount that these files do not show.

## What the page finds

1. The drugged-only model reads AUROC 0.717 [0.689, 0.747] and AUPRC 0.470 [0.422, 0.525] out of fold on 1,564 drugged genes (415 flagged, base rate 0.265).
2. The outcome model of the DR+PU pipeline reads AUROC 0.727 [0.701, 0.754] and AUPRC 0.454 [0.411, 0.509] on 1,506 drugged genes, and its AUROC minus that of the drugged-only model on the 1,506 genes that both score is +0.013 [-0.008, +0.035], an interval that holds 0.
3. By stage the drugged-only AUROC reads 0.577 at code 0.5, 0.699 at code 0.75 and 0.752 at code 1.0, and the cell of code 0.25 holds 28 genes with 8 flagged in the drugged-only model and 27 genes with 8 flagged in step 3, so it shows counts only.
4. The deployed DR+PU score cannot be evaluated on drugged genes: every flagged drugged gene scores above every unflagged one in sample (AUROC 1.000), because the pseudo-outcome it is fitted to holds the label.
5. On the 16,239 undrugged genes, whose labels no step read, the deployed score reads AUROC 0.636 [0.611, 0.662] and the outcome model 0.722 [0.698, 0.744] for the 508 flagged genes (base rate 0.031), against unlabeled genes that are unknown and not safe.
6. 10 of the 16 numbers that the written docs state about the outcome model, the prior, the pseudo-outcome and the labels differ from the saved files and the code (section 8).
7. Cancer only is not run, combination does not apply to genes, the files name no drug and no trial, and no DeepSeek V4 Flash arm exists (sections 3 and 6).

## 1. The experiments

The scripts hold a drugged-only model, an all-labeled model and the DR+PU pipeline, and the all-labeled model has no saved output.

| Reading | Question | Label | Genes | Columns | Evaluation | Saved in the repo |
|---|---|---|---:|---|---|---|
| Drugged-only model (`train_drugged_only`, Model B) | Which genes that entered a trial carry a safety flag? | flag against no flag | 1,564 | 71 saved columns: 14 Open Targets, 1 LOEUF, 18 pLoF and 38 expression | 5-fold out of fold, seed 42 | `predictions_oof.csv`, `shap_summary.csv`, `shap_values.npy` |
| All-labeled model (`train_and_evaluate`, Model A) | Which targets carry a flag, learned from flagged genes and drugged genes without a flag? | flag, or drugged without a flag; the rest is left out | 2,094 trained, 78,725 scored | built by the script, with the same four groups | 5-fold, printed only | none: its four output files are not in the repo |
| DR+PU step 1, propensity | Which genes enter a trial? | entered a trial against not | 17,745 | 191 columns, as `README.md` states | 5-fold cross-fitted | `pi_hat` of `dr/predictions.csv` |
| DR+PU step 2, prior | What share of genes carry a liability? | none | none | none | a constant, 0.40, not an estimate | none |
| DR+PU step 3, outcome model | Which drugged genes carry a flag, with weights for selection and for unlabeled negatives? | flag against no flag | 1,506 | the same 191 | 5-fold cross-fitted | `m_hat` of `dr/predictions.csv` |
| DR+PU step 5, deployed score | A score for any gene, fitted to the pseudo-outcome of step 4 | the pseudo-outcome, which holds the label of a drugged gene | 17,745 | the same 191 | none out of fold | `safety_score_dr` of `dr/predictions.csv` |
| Single Open Targets columns (made by this page, no model) | How far does one column go on its own? | flag against no flag | 1,564 | 14 columns, one at a time | no fitting, so no folds | `docs/grid/numbers.csv` |

All learners are XGBoost.
The drugged-only classifier grows up to 500 trees of depth 5 with learning rate 0.05, a row and a column subsample of 0.8 and 0.8, and stops after 30 rounds without a gain in AUPRC.
The all-labeled classifier differs in its depth, 6.
The DR+PU classifiers share the settings of the drugged-only classifier, and its final regressor has 300 trees.
All of them use the seed 42 and 5 stratified folds with shuffling.
This page adds one reading that no script holds, the single Open Targets columns of section 5, and it fits no model.

## 2. What is counted

| Word | What it counts |
|---|---|
| Gene | One Ensembl id of the Open Targets file, which is the row of every table of this page. |
| Flag | `hasSafetyEvent` is `-1`: Open Targets lists a known safety event for the target. No other value occurs, so a gene without a flag is not confirmed free of one. |
| Drugged | A gene with a value in `maxClinicalTrialPhase`, the maximum phase of its record. The code reads a gene with no value as not drugged. |
| positive | A flagged gene, whether or not it entered a trial (the label `positive` of the saved predictions). |
| drugged_safe | A drugged gene without a flag. It is an assumed negative, not a confirmed one. |
| unlabeled | A gene with no phase and no flag. Its label is unknown, and it is not safe by default. |
| Events | The genes whose label is the flag. The base rate is the events over the genes, which is also the AUPRC of a random ranking. |
| Out of fold | A gene scored by a model that did not train on it. |

### The labels

The code builds the labels from two columns of the Open Targets file, `hasSafetyEvent` and `maxClinicalTrialPhase`.
The first holds `-1` for a flagged gene and an empty cell for every other gene, so the file never says that a gene is safe.
The table counts the genes of each universe by label.

| Universe | Genes | positive | of which drugged | of which no phase | drugged_safe | unlabeled |
|---|---:|---:|---:|---:|---:|---:|
| Open Targets file | 78,725 | 945 | 415 | 530 | 1,149 | 76,631 |
| Drugged-only rows | 1,564 | 415 | 415 | 0 | 1,149 | 0 |
| DR+PU file | 17,745 | 917 | 409 | 508 | 1,097 | 15,731 |

The drugged-only model keeps the genes that have a phase, so it holds no unlabeled gene.
The all-labeled model trains on the flagged genes, including those that never entered a trial, and on the drugged genes without a flag, and it scores every target.
The saved DR+PU file calls a gene positive when it is flagged, whether or not it entered a trial, so 508 of its 917 positives never entered a trial.

### Genes used to train and to score

| Reading | Genes in the universe | Trained on | Scored | Per fold, trained | Per fold, scored | Events | Base rate |
|---|---:|---:|---:|---:|---:|---:|---:|
| Drugged-only model | 78,725 | 1,564 | 1,564 | 1,251 to 1,252 | 312 to 313 | 415 | 0.265 |
| All-labeled model (no saved result) | 78,725 | 2,094 | 78,725 | 1,675 to 1,676 | 418 to 419 | 945 | 0.451 |
| DR+PU step 1, propensity | 17,745 | 17,745 | 17,745 | 14,196 | 3,549 | 1,506 | 0.085 |
| DR+PU step 3, outcome model | 17,745 | 1,506 | 1,506 | 1,204 to 1,205 | 301 to 302 | 409 | 0.272 |
| DR+PU step 5, deployed score | 17,745 | 17,745 | 17,745 | no folds | no folds | no events | none |

The folds are not stored in the files.
Their sizes come from the fold rule of the code (class counts, 5 folds, seed 42), which fixes how many genes each fold holds.
The final drugged-only model, fitted to all its genes, is not saved either, but its SHAP values are.

### Drugs and trials

The Open Targets file has 17 columns: an Ensembl id, 14 features, the safety flag and the maximum phase.
None of them is a drug, a trial, an indication or a date, and no column of the saved prediction files is one either.
The data files that git tracks are this file, the saved predictions and SHAP values, and two small lists under `data/demo` (gene symbols, and organs with their cell types).
The number of distinct drugs or trials behind the labels is therefore not in this repo's files, and this page states none.

### Deduplication

| File | Rows | Unique Ensembl ids | Duplicate ids | Rows with a symbol | Unique symbols | Duplicate symbols |
|---|---:|---:|---:|---:|---:|---:|
| Open Targets parquet | 78,725 | 78,725 | 0 | not in the file | not in the file | not in the file |
| Drugged-only predictions | 1,564 | 1,564 | 0 | 73 | 73 | 0 |
| DR+PU predictions | 17,745 | 17,745 | 0 | 17,745 | 17,745 | 0 |
| data/demo/genes.json | 19,155 | not in the file | not in the file | 19,155 | 19,155 | 0 |

Every file has one row per Ensembl id, and the files that carry symbols have one row per symbol.
No gene family is collapsed, and the folds are not grouped by family (the scripts split with `StratifiedKFold` only), so a gene and a paralog can fall in different folds.
`paralogMaxIdentityPercentage` is a feature of the models, and the files hold no family id to group by.
The list `data/demo/genes.json` holds 19,155 unique symbols, which is the count that `Methodology.md` writes as ~19,155 (its section 1), and all 17,745 symbols of the DR+PU file are in it.

## 3. The cuts

A cut keeps genes by a fact that the files hold, and no model reads it.

| Cut | What it keeps | In these files |
|---|---|---|
| each gene once | one row per Ensembl id | holds in every file (section 2) |
| each drug once | a drug counts once | not applicable: the files name no drug |
| cancer only, not cancer | by the disease of the trial | not run: the files hold no indication (section 3) |
| single agent, combination, mixed | by the regimen | not applicable: a row is a gene, and a gene has no regimen |
| stage | by the maximum phase | counts for every universe and results for the 1,564 drugged genes (sections 3 and 4) |
| driver gene | by `isCancerDriverGene` | counts only, and it is no cancer cut (section 3) |

### Stage

The stage of a gene is `maxClinicalTrialPhase`, the maximum phase of its record.
`Methodology.md` (line 161) says that the codes 0.25, 0.5, 0.75 and 1.0 stand for Phases I to IV, and an empty cell means that the gene never entered a trial.
The first table counts the genes of the Open Targets file, and its four coded rows are also the 1,564 genes of the drugged-only model.
The second table counts the genes of the DR+PU file.

The Open Targets file:

| Max phase code | Phase | Genes | positive | drugged_safe | unlabeled | Flag rate |
|---|---|---:|---:|---:|---:|---:|
| 0.25 | I | 28 | 8 | 20 | 0 | 0.286 |
| 0.5 | II | 296 | 59 | 237 | 0 | 0.199 |
| 0.75 | III | 222 | 61 | 161 | 0 | 0.275 |
| 1.0 | IV | 1,018 | 287 | 731 | 0 | 0.282 |
| none | not drugged | 77,161 | 530 | 0 | 76,631 | 0.007 |
| all | all | 78,725 | 945 | 1,149 | 76,631 | 0.012 |

The DR+PU file:

| Max phase code | Phase | Genes | positive | drugged_safe | unlabeled | Flag rate |
|---|---|---:|---:|---:|---:|---:|
| 0.25 | I | 27 | 8 | 19 | 0 | 0.296 |
| 0.5 | II | 268 | 56 | 212 | 0 | 0.209 |
| 0.75 | III | 212 | 61 | 151 | 0 | 0.288 |
| 1.0 | IV | 999 | 284 | 715 | 0 | 0.284 |
| none | not drugged | 16,239 | 508 | 0 | 15,731 | 0.031 |
| all | all | 17,745 | 917 | 1,097 | 15,731 | 0.052 |

### Cancer only

Cancer only is not run.
It needs a field that says in which disease each drug was tested, such as an indication or an area for each target, and the files hold none.
The one cancer field is `isCancerDriverGene`, a property of the gene, so a cut on it asks whether the target is a driver and not whether the drug treats a cancer.
In the Open Targets file 334 genes carry the flag, 239 of them never entered a trial, and the drugged genes split as follows.

| Cancer driver flag | Drugged genes | positive | drugged_safe | Flag rate |
|---|---:|---:|---:|---:|
| flagged as driver | 95 | 44 | 51 | 0.463 |
| not flagged as driver | 1,469 | 371 | 1,098 | 0.253 |
| all drugged genes | 1,564 | 415 | 1,149 | 0.265 |

### Combination therapy

Single agent against combination does not apply, because a row is a gene.
A gene has no regimen, and the files name no drug.

## 4. The results

Each row shows the number of genes, the number of flagged genes (the events) and the base rate, next to the AUROC and the AUPRC with their 95% intervals.
A cell with fewer than 30 genes, or with fewer than 10 flagged or unflagged genes, shows its counts only.

| Reading | Genes | Events | Base rate | AUROC [95% interval] | AUPRC [95% interval] |
|---|---:|---:|---:|---|---|
| Drugged-only model, out of fold | 1,564 | 415 | 0.265 | 0.717 [0.689, 0.747] | 0.470 [0.422, 0.525] |
| DR+PU step 3 outcome model, out of fold | 1,506 | 409 | 0.272 | 0.727 [0.701, 0.754] | 0.454 [0.411, 0.509] |
| Difference, step 3 minus drugged-only, on the 1,506 genes that both score (paired) | 1,506 | 409 | 0.272 | +0.013 [-0.008, +0.035] | -0.019 [-0.062, +0.026] |
| DR+PU step 1 propensity (entered a trial or not), out of fold | 17,745 | 1,506 | 0.085 | 0.915 [0.907, 0.923] | 0.600 [0.574, 0.625] |
| DR+PU outcome model on undrugged genes, which it never trained on | 16,239 | 508 | 0.031 | 0.722 [0.698, 0.744] | 0.088 [0.076, 0.106] |
| DR+PU deployed score on undrugged genes | 16,239 | 508 | 0.031 | 0.636 [0.611, 0.662] | 0.073 [0.061, 0.089] |
| DR+PU deployed score on drugged genes, in sample (not a result) | 1,506 | 409 | 0.272 | 1.000 in sample, no interval | 1.000 in sample, no interval |

The rows read as follows.

1. The drugged-only row scores every drugged gene once out of fold, from the saved `predictions_oof.csv`.
2. The step 3 row is the outcome model of the DR+PU pipeline on its 1,506 drugged genes, from the column `m_hat` of `predictions.csv`, which holds the out-of-fold score of a drugged gene.
3. The paired row repeats both readings on the 1,506 genes that both files score, on the same bootstrap draws.
   The two models differ in their columns, their weights and their folds, so the difference is not the effect of one change.
4. The step 1 row is the propensity model, which predicts whether a gene entered a trial.
   Its target is not a safety label, and it is shown because the weights of the pipeline come from it.
5. The two undrugged rows are the one check of the pipeline in which no step read the label of the genes that are scored.
   The 16,239 undrugged genes were not used to fit the outcome model, and the labels of the 508 that carry a flag were read by no step, because the prior of step 2 is fixed in the code (line 61 of `data/safety_model_dr.py`) and step 4 gives an undrugged gene the score of the outcome model as its target.
   The deployed regressor of step 5 is still fitted on the features of these genes, so the check is blind to their labels and not to the genes.
   The 15,731 other undrugged genes are unlabeled, which means unknown and not safe, so every unknown flagged gene counts as a miss and the check understates a perfect ranking.
   The deployed score is clipped to `[0, 1]` (line 432 of `data/safety_model_dr.py`), so 9,650 of the 16,239 undrugged genes tie at `0` (205 of them flagged), and the AUROC counts a tie as half.
6. The last row is not a result.
   The deployed regressor is fitted to the pseudo-outcome of step 4, which adds the label of a drugged gene to its score, so every flagged drugged gene scores above every unflagged one in sample (the lowest flagged score is 0.483 and the highest unflagged is 0.351).
   The saved deployed score of a drugged gene therefore carries its own label, and it cannot be evaluated out of fold.

### By stage

| Reading | Max phase code | Phase | Genes | Events | Base rate | AUROC [95% interval] | AUPRC [95% interval] |
|---|---|---|---:|---:|---:|---|---|
| Drugged-only model, out of fold | 0.25 | I | 28 | 8 | 0.286 | counts only | counts only |
| Drugged-only model, out of fold | 0.5 | II | 296 | 59 | 0.199 | 0.577 [0.494, 0.662] | 0.263 [0.190, 0.376] |
| Drugged-only model, out of fold | 0.75 | III | 222 | 61 | 0.275 | 0.699 [0.622, 0.772] | 0.445 [0.344, 0.584] |
| Drugged-only model, out of fold | 1.0 | IV | 1,018 | 287 | 0.282 | 0.752 [0.721, 0.784] | 0.535 [0.477, 0.596] |
| DR+PU step 3 outcome model, out of fold | 0.25 | I | 27 | 8 | 0.296 | counts only | counts only |
| DR+PU step 3 outcome model, out of fold | 0.5 | II | 268 | 56 | 0.209 | 0.615 [0.536, 0.693] | 0.276 [0.210, 0.385] |
| DR+PU step 3 outcome model, out of fold | 0.75 | III | 212 | 61 | 0.288 | 0.685 [0.604, 0.761] | 0.440 [0.334, 0.563] |
| DR+PU step 3 outcome model, out of fold | 1.0 | IV | 999 | 284 | 0.284 | 0.762 [0.731, 0.794] | 0.518 [0.464, 0.582] |

The cell of code 0.25 holds 28 genes with 8 flagged in the drugged-only model and 27 genes with 8 flagged in step 3, so it shows counts only.
Each stage is read on fewer genes than the whole set, so its intervals are wider.

## 5. The features that matter

### The drugged-only model

The saved `shap_summary.csv` holds the mean absolute SHAP value of each of the 71 columns of the model, and `shap_values.npy` holds the matrix behind it (1,564 genes by 71 columns), which reproduces the summary.

| Group | Columns | Of which | Sum of mean absolute SHAP | Share | Largest column |
|---|---:|---|---:|---:|---|
| Open Targets | 14 | 14 columns | 4.1067 | 95.4% | `mouseKOScore` 0.7683 |
| LOEUF | 1 | 1 column | 0.0088 | 0.2% | `loeuf_score` 0.0088 |
| pLoF | 18 | 4 global and 14 by organ system | 0.0775 | 1.8% | `plof_respiratory_min_p` 0.0166 |
| Expression | 38 | 24 cell type (12 mean, 12 percent expressing) and 14 by organ | 0.1105 | 2.6% | `expr_organ_Heart` 0.0163 |

The Open Targets columns carry 95.4% of the total, and the 57 other columns together carry 4.6%.
The scripts fill LOEUF, pLoF and expression columns only for the 108 genes of the fixed list `GENES` (line 74 of `data/safety_model.py`), and give a symbol only to those genes (line 305).
The saved predictions carry a symbol for 73 rows, and these are exactly the 73 drugged genes of that list.
So the three groups other than Open Targets were very likely empty for the other 1,491 genes, although the saved files do not store the feature values to confirm it.
Of the 38 expression columns of the saved model, 24 carry names that the current `get_feature_cols` does not select (it takes only the prefixes `expr_ct_` and `expr_organ_`), so these files come from an earlier version of the script.

| Rank | Column | Group | Mean absolute SHAP | Share of the total |
|---:|---|---|---:|---:|
| 1 | `mouseKOScore` | Open Targets | 0.7683 | 17.9% |
| 2 | `geneticConstraint` | Open Targets | 0.5274 | 12.3% |
| 3 | `mouseOrthologMaxIdentityPercentage` | Open Targets | 0.5171 | 12.0% |
| 4 | `hasSmallMoleculeBinder` | Open Targets | 0.4403 | 10.2% |
| 5 | `hasHighQualityChemicalProbes` | Open Targets | 0.3291 | 7.6% |
| 6 | `paralogMaxIdentityPercentage` | Open Targets | 0.3249 | 7.5% |
| 7 | `hasLigand` | Open Targets | 0.3152 | 7.3% |
| 8 | `tissueSpecificity` | Open Targets | 0.2230 | 5.2% |
| 9 | `hasPocket` | Open Targets | 0.2124 | 4.9% |
| 10 | `tissueDistribution` | Open Targets | 0.1743 | 4.1% |
| 11 | `isInMembrane` | Open Targets | 0.1655 | 3.8% |
| 12 | `isSecreted` | Open Targets | 0.0889 | 2.1% |
| 13 | `plof_respiratory_min_p` | pLoF | 0.0166 | 0.4% |
| 14 | `expr_organ_Heart` | Expression | 0.0163 | 0.4% |
| 15 | `isCancerDriverGene` | Open Targets | 0.0159 | 0.4% |

The SHAP values explain the final model, which was fitted to all 1,564 genes and reads its own labels almost perfectly in sample (AUROC 0.9996 from the summed SHAP values of each gene).
The rows of the matrix follow the order of the drugged genes in the parquet, which the rank correlation of 0.637 between the summed SHAP values and the out-of-fold scores supports (a shuffled copy gives -0.008).
The table therefore says which columns the model uses, and it is not an out-of-fold measure of signal.

### One column alone

For each Open Targets column of the drugged-only model, the next table ranks the drugged genes by that column alone and reads the AUROC against the flag, with no model.
A column that holds one value only (`isCancerDriverGene` and `hasTEP`) is read as present against empty, and any other column leaves out the genes with no value.
A column with fewer than 30 present or ranked genes, or with fewer than 10 flagged or unflagged ones, shows its counts only.

| Column | Basis | Genes | Flagged | AUROC alone [95% interval] | Reading | SHAP rank | Mean absolute SHAP |
|---|---|---:|---:|---|---|---:|---:|
| `mouseKOScore` | genes with a value | 1,300 | 388 | 0.377 [0.346, 0.410] | lower value goes with a flag | 1 | 0.7683 |
| `geneticConstraint` | genes with a value | 1,470 | 401 | 0.451 [0.417, 0.483] | lower value goes with a flag | 2 | 0.5274 |
| `mouseOrthologMaxIdentityPercentage` | genes with a value | 1,485 | 406 | 0.523 [0.493, 0.553] | no ranking shown | 3 | 0.5171 |
| `hasSmallMoleculeBinder` | genes with a value | 1,545 | 413 | 0.636 [0.617, 0.655] | higher value goes with a flag | 4 | 0.4403 |
| `hasHighQualityChemicalProbes` | genes with a value | 502 | 222 | 0.520 [0.480, 0.562] | no ranking shown | 5 | 0.3291 |
| `paralogMaxIdentityPercentage` | genes with a value | 1,355 | 389 | 0.490 [0.462, 0.516] | no ranking shown | 6 | 0.3249 |
| `hasLigand` | genes with a value | 1,545 | 413 | 0.635 [0.612, 0.656] | higher value goes with a flag | 7 | 0.3152 |
| `tissueSpecificity` | genes with a value | 1,523 | 411 | 0.521 [0.490, 0.551] | no ranking shown | 8 | 0.2230 |
| `hasPocket` | genes with a value | 1,545 | 413 | 0.591 [0.567, 0.617] | higher value goes with a flag | 9 | 0.2124 |
| `tissueDistribution` | genes with a value | 1,523 | 411 | 0.536 [0.507, 0.567] | higher value goes with a flag | 10 | 0.1743 |
| `isInMembrane` | genes with a value | 1,539 | 409 | 0.560 [0.531, 0.587] | higher value goes with a flag | 11 | 0.1655 |
| `isSecreted` | genes with a value | 1,539 | 409 | 0.478 [0.456, 0.499] | lower value goes with a flag | 12 | 0.0889 |
| `isCancerDriverGene` | present against empty | 1,564 (95 present) | 415 | 0.531 [0.515, 0.547] | present goes with a flag | 15 | 0.0159 |
| `hasTEP` | present against empty | 1,564 (8 present) | 415 | counts only | counts only | 26 | 0.0043 |

The best of the 14 columns is `hasSmallMoleculeBinder`, which reads AUROC 0.636 [0.617, 0.655] on its own in its best direction, against 0.717 [0.689, 0.747] for the drugged-only model.
The best column is picked from 14 by its point value, so its interval is a little optimistic.

### The DR+PU pipeline

The importance of the DR+PU model is saved only as pictures (`data/safety_model_output/dr/figures/04_feature_importance.png` and two more), not as a table, and its model files and its SHAP file are not in the repo.
`README.md` (line 425) states 191 columns in five groups (Open Targets platform features 14, LOEUF constraint score 1, pLoF Genebass aggregates 42, cell-type expression 120 and organ-level max expression 14), which add up to that total.
This page cannot recount them, because the matrix is built from `backend/humanproof.db`, which is not in the repo.

## 6. The language model only arm: not run

The arm would ask DeepSeek V4 Flash about a gene alone, with no model of this repo, and would grade it only on outcomes that come after its release date, 2026-04-24.
That date is an input of this page and is not read from a file of this repo.
No such arm exists in this repo, and no model was called to make this page.
The table gives the reasons from the files.

| Reason | What the files give | What an arm would need |
|---|---|---|
| The labels have no date | `hasSafetyEvent` is a flag of release 25.12, and none of the 17 columns of the file is a date | flags dated after 2026-04-24 |
| The snapshot ends before the release | release 25.12 is named for December 2025, which ends before 2026-04-24 | a later Open Targets release, which is not in this repo |
| The textbook liabilities are public | 8 of 8 well known safety targets carry the flag (`KCNH2`, `HTR2B`, `SCN5A`, `PTGS1`, `PTGS2`, `ADRB2`, `DRD2`, `HRH1`) | an outcome that a language model has not read |
| No forward list exists | no forecast of this repo was frozen before a label was read | forecasts frozen before the later release is read, graded once on the genes it newly flags |
| Nothing paid was run | no model was called to make this page | a recorded spend |

The label is a flag in one release of Open Targets, not an outcome that happens after a forecast day.
A forecast of it would test what a language model recalls about known safety targets, and the list above shows that these targets are public.
The list is a short fixed list written into the script, an illustration and not a sample.
An arm would need a later Open Targets release, forecasts frozen before that release is read, and a grade on the genes that the release newly flags.

## 7. The cells that were not run

Each row is a cell that a reader may expect and that this page did not fill, with its reason.
A file counts as absent when git does not track it, and every absent file named on this page is excluded by a rule of `.gitignore`, which keeps databases, pickles, compressed tables, model files and the outputs of the all-labeled model out of the repo.
A maintainer may therefore hold some of them locally, and this page does not see them.

| Cell | Why it is not run | What it would need |
|---|---|---|
| Result of the all-labeled model (Model A) | its output files are not in the repo (4 of 4 absent), and no run of it is saved | the saved output, or a run of the script |
| Retraining any model here | 7 of the 7 input files that the scripts read are absent, among them `backend/humanproof.db` | those files, which `README.md` says how to fetch and build |
| Feature group ablations (drop Open Targets, genetics or expression) | no script holds one, and the saved model has genetics and expression values for few genes | the feature matrix, built from `backend/humanproof.db`, and XGBoost |
| A split in time or a held-back test | the labels carry no dates, and every number is out of fold | dated labels, such as two Open Targets releases, and a gene list fixed before the labels are read |
| Cancer only and not cancer | the files hold no indication of any trial | a field that names the disease of the trial, or an area for each target |
| Single agent, combination, mixed | a row is a gene, and the files name no drug | a table of the drugs and their regimens for each target |
| Counts of distinct drugs and trials | no column of the files names a drug or a trial | a table of the drugs and the trials behind each flag and each phase |
| Folds that group paralogs | no family id is in the files, and the scripts split with `StratifiedKFold` only | a family id for each gene |
| Out-of-fold evaluation of the deployed score on drugged genes | the regressor is fitted to the pseudo-outcome, which holds the label of these genes | a cross-fitted step 5 |
| The prior of step 2 as an estimate | the code fixes it at 0.40 | a run with `PI_OVERRIDE` set to `None` |
| A table of the DR+PU feature importance | it is saved as pictures only, and 3 of 3 model and SHAP files are absent | `gene_shap_dr.json`, or a run of the script |
| The DeepSeek V4 Flash only arm | see section 6 | flags dated after 2026-04-24 and forecasts frozen before they are read |

## 8. Where the written docs, the code and the saved files differ

The table sets 16 numbers that the written docs state about the DR+PU pipeline against what the saved files and the code give.
A stated number agrees when it lies within one and a half units of its last digit of the recomputed one.

| Statement | Where | Stated | Saved files and code give | Verdict |
|---|---|---|---|---|
| Step 3 out-of-fold AUROC and AUPRC | `README.md` line 440 | 0.691 and 0.426 | 0.727 and 0.454 | differs |
| Step 3 out-of-fold AUROC and AUPRC | `Methodology.md` line 382 | 0.691 and 0.426 | 0.727 and 0.454 | differs |
| Step 3 out-of-fold AUROC and AUPRC | `CONTRIBUTING.md` line 146 | 0.691 and 0.426 | 0.727 and 0.454 | differs |
| Step 3 out-of-fold AUROC and AUPRC | `frontend/src/app/methodology/page.tsx` line 50 | 0.727 and 0.454 | 0.727 and 0.454 | agrees |
| Prior of step 2, the share of genes with a liability | `Methodology.md` line 350 | 0.090 | 0.40 fixed in the code (line 61), not stored in the saved files | differs |
| Prior of step 2, the share of genes with a liability | `README.md` line 439 | 9% | 0.40 fixed in the code (line 61), not stored in the saved files | differs |
| Prior of step 2, the share of genes with a liability | `CONTRIBUTING.md` line 145 | 9% | 0.40 fixed in the code (line 61), not stored in the saved files | differs |
| Prior of step 2, the share of genes with a liability | `frontend/src/app/methodology/page.tsx` line 37 | 40% | 0.40 fixed in the code (line 61), not stored in the saved files | agrees |
| Mean and standard deviation of the pseudo-outcome, drugged genes | `Methodology.md` line 419 | 0.17 and 1.00 | 0.01 and 1.02 | differs |
| Mean and standard deviation of the pseudo-outcome, undrugged genes | `Methodology.md` line 420 | 0.02 and 0.06 | 0.04 and 0.10 | differs |
| Mean and standard deviation of the pseudo-outcome, drugged genes | `frontend/src/app/methodology/page.tsx` line 64 | 0.01 and 1.02 | 0.01 and 1.02 | agrees |
| Mean and standard deviation of the pseudo-outcome, undrugged genes | `frontend/src/app/methodology/page.tsx` line 64 | 0.04 and 0.10 | 0.04 and 0.10 | agrees |
| Genes labeled positive among the DB genes | `Methodology.md` line 491 | 409 | 917 | differs |
| Genes labeled drugged_safe among the DB genes | `Methodology.md` line 492 | 1,097 | 1,097 | agrees |
| Genes labeled unlabeled among the DB genes | `Methodology.md` line 493 | 16,239 | 15,731 | differs |
| Propensity AUROC | `Methodology.md` line 321 | 0.91 | 0.915 | agrees |

Of these 16 numbers, 10 differ.
The saved outcome model reads 0.727 and 0.454, which are the numbers that the app page (`frontend/src/app/methodology/page.tsx`) states, and the saved summaries of the pseudo-outcome also match that page and not `Methodology.md`.
The code fixes the prior at 0.40 (line 61 of `data/safety_model_dr.py`) and the app page says so, while `Methodology.md` (line 350) and `README.md` still describe the Elkan-Noto estimate, which `Methodology.md` writes as 0.090.
The comment of the script (lines 57 and 58 of `data/safety_model_dr.py`) and the app page (line 35) call the fixed prior a monotone shift that leaves AUROC and AUPRC unchanged, and the app page also says that it raised the out-of-fold AUROC from 0.691 to 0.727.
Both cannot hold for a monotone shift, and in the code the prior is a sample weight of step 3 (line 208), so another prior refits the trees and can change the ranking.
The saved files hold the run with the fixed prior only, so they cannot show what the prior changed.
This page edits neither the scripts nor the written docs.
