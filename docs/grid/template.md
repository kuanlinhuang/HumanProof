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

{s_summary}

### Unit

{s_unit}

### Label

{s_label}

### Population

{s_population}

### Test

{s_test}

### Metric

{s_metric}

### Baseline

{s_baseline}

### Not comparable

{s_not_comparable}

### How many

{s_how_many}

### Cuts

{s_cuts}

### Features that matter

{s_features}

### Language model only

{s_llm}

## Status of the evidence

All results are out-of-fold numbers on the labels of one Open Targets snapshot, release {ot_release}.
No result is blind: the snapshot carries no dates, so there is no split in time, and no test set was held back.
The intervals are 95% percentile intervals of {draws} bootstrap draws over genes at seed {seed}, so they do not count the choice of settings or a repeat of the folds.
XGBoost stops adding trees when the AUPRC of the scored fold stops improving (line {e1_stop_line} of `data/safety_model.py` and line {dr_stop_line} of `data/safety_model_dr.py`), so the tree count is picked on the genes that are then scored.
The out-of-fold numbers are therefore a little optimistic, by an amount that these files do not show.

## What the page finds

1. {f1}
2. {f2}
3. {f3}
4. {f4}
5. {f5}
6. {f6}
7. {f7}

## 1. The experiments

The scripts hold a drugged-only model, an all-labeled model and the DR+PU pipeline, and the all-labeled model has no saved output.

{table_experiments}

All learners are XGBoost.
The drugged-only classifier grows up to {e1_trees} trees of depth {e1_depth} with learning rate {e1_lr}, a row and a column subsample of {e1_sub} and {e1_col}, and stops after {e1_stop} rounds without a gain in AUPRC.
The all-labeled classifier differs in its depth, {e2_depth}.
The DR+PU classifiers share the settings of the drugged-only classifier, and its final regressor has {dr_reg_trees} trees.
All of them use the seed {seed} and {folds} stratified folds with shuffling.
This page adds one reading that no script holds, the single Open Targets columns of section 5, and it fits no model.

## 2. What is counted

{table_words}

### The labels

The code builds the labels from two columns of the Open Targets file, `hasSafetyEvent` and `maxClinicalTrialPhase`.
The first holds `-1` for a flagged gene and an empty cell for every other gene, so the file never says that a gene is safe.
The table counts the genes of each universe by label.

{table_labels}

The drugged-only model keeps the genes that have a phase, so it holds no unlabeled gene.
The all-labeled model trains on the flagged genes, including those that never entered a trial, and on the drugged genes without a flag, and it scores every target.
The saved DR+PU file calls a gene positive when it is flagged, whether or not it entered a trial, so {dr_und_flag} of its {dr_flag} positives never entered a trial.

### Genes used to train and to score

{table_setup}

The folds are not stored in the files.
Their sizes come from the fold rule of the code (class counts, {folds} folds, seed {seed}), which fixes how many genes each fold holds.
The final drugged-only model, fitted to all its genes, is not saved either, but its SHAP values are.

### Drugs and trials

The Open Targets file has {n_cols} columns: an Ensembl id, {n_features} features, the safety flag and the maximum phase.
None of them is a drug, a trial, an indication or a date, and no column of the saved prediction files is one either.
The data files that git tracks are this file, the saved predictions and SHAP values, and two small lists under `data/demo` (gene symbols, and organs with their cell types).
The number of distinct drugs or trials behind the labels is therefore not in this repo's files, and this page states none.

### Deduplication

{table_dedup}

Every file has one row per Ensembl id, and the files that carry symbols have one row per symbol.
No gene family is collapsed, and the folds are not grouped by family (the scripts split with `StratifiedKFold` only), so a gene and a paralog can fall in different folds.
`paralogMaxIdentityPercentage` is a feature of the models, and the files hold no family id to group by.
The list `data/demo/genes.json` holds {gene_list_n} unique symbols, which is the count that `Methodology.md` writes as {meth_universe} (its section 1), and all {dr_n} symbols of the DR+PU file are in it.

## 3. The cuts

A cut keeps genes by a fact that the files hold, and no model reads it.

{table_cut_words}

### Stage

The stage of a gene is `maxClinicalTrialPhase`, the maximum phase of its record.
`Methodology.md` (line {phase_line}) says that the codes {phase_codes} stand for Phases I to IV, and an empty cell means that the gene never entered a trial.
The first table counts the genes of the Open Targets file, and its four coded rows are also the {e1_n} genes of the drugged-only model.
The second table counts the genes of the DR+PU file.

The Open Targets file:

{table_stage_ot}

The DR+PU file:

{table_stage_dr}

### Cancer only

Cancer only is not run.
It needs a field that says in which disease each drug was tested, such as an indication or an area for each target, and the files hold none.
The one cancer field is `isCancerDriverGene`, a property of the gene, so a cut on it asks whether the target is a driver and not whether the drug treats a cancer.
In the Open Targets file {drivers_all} genes carry the flag, {drivers_undrugged} of them never entered a trial, and the drugged genes split as follows.

{table_driver}

### Combination therapy

Single agent against combination does not apply, because a row is a gene.
A gene has no regimen, and the files name no drug.

## 4. The results

Each row shows the number of genes, the number of flagged genes (the events) and the base rate, next to the AUROC and the AUPRC with their 95% intervals.
A cell with fewer than {min_genes} genes, or with fewer than {min_class} flagged or unflagged genes, shows its counts only.

{table_results}

The rows read as follows.

1. The drugged-only row scores every drugged gene once out of fold, from the saved `predictions_oof.csv`.
2. The step 3 row is the outcome model of the DR+PU pipeline on its {dr3_n} drugged genes, from the column `m_hat` of `predictions.csv`, which holds the out-of-fold score of a drugged gene.
3. The paired row repeats both readings on the {shared_n} genes that both files score, on the same bootstrap draws.
   The two models differ in their columns, their weights and their folds, so the difference is not the effect of one change.
4. The step 1 row is the propensity model, which predicts whether a gene entered a trial.
   Its target is not a safety label, and it is shown because the weights of the pipeline come from it.
5. The two undrugged rows are the one check of the pipeline in which no step read the label of the genes that are scored.
   The {und_n} undrugged genes were not used to fit the outcome model, and the labels of the {und_pos} that carry a flag were read by no step, because the prior of step 2 is fixed in the code (line {pi_line} of `data/safety_model_dr.py`) and step 4 gives an undrugged gene the score of the outcome model as its target.
   The deployed regressor of step 5 is still fitted on the features of these genes, so the check is blind to their labels and not to the genes.
   The {und_neg} other undrugged genes are unlabeled, which means unknown and not safe, so every unknown flagged gene counts as a miss and the check understates a perfect ranking.
   The deployed score is clipped to `{clip_range}` (line {clip_line} of `data/safety_model_dr.py`), so {und_zero} of the {und_n} undrugged genes tie at `{clip_lo}` ({und_flag_zero} of them flagged), and the AUROC counts a tie as half.
6. The last row is not a result.
   The deployed regressor is fitted to the pseudo-outcome of step 4, which adds the label of a drugged gene to its score, so every flagged drugged gene scores above every unflagged one in sample (the lowest flagged score is {min_flag_score} and the highest unflagged is {max_safe_score}).
   The saved deployed score of a drugged gene therefore carries its own label, and it cannot be evaluated out of fold.

### By stage

{table_stage_results}

The cell of {small_text}, so it shows counts only.
Each stage is read on fewer genes than the whole set, so its intervals are wider.

## 5. The features that matter

### The drugged-only model

The saved `shap_summary.csv` holds the mean absolute SHAP value of each of the {n_saved} columns of the model, and `shap_values.npy` holds the matrix behind it ({shap_rows} genes by {shap_cols} columns), which reproduces the summary.

{table_groups}

The Open Targets columns carry {ot_share} of the total, and the {n_nonot} other columns together carry {nonot_share}.
The scripts fill LOEUF, pLoF and expression columns only for the {genes_list_n108} genes of the fixed list `GENES` (line {genes_line} of `data/safety_model.py`), and give a symbol only to those genes (line {symbol_line}).
The saved predictions carry a symbol for {n_sym} rows, and these are exactly the {n_sym_drugged} drugged genes of that list.
So the three groups other than Open Targets were very likely empty for the other {n_nosym} genes, although the saved files do not store the feature values to confirm it.
Of the {n_expr} expression columns of the saved model, {n_expr_old} carry names that the current `get_feature_cols` does not select (it takes only the prefixes `expr_ct_` and `expr_organ_`), so these files come from an earlier version of the script.

{table_shap_top}

The SHAP values explain the final model, which was fitted to all {e1_n} genes and reads its own labels almost perfectly in sample (AUROC {shap_insample} from the summed SHAP values of each gene).
The rows of the matrix follow the order of the drugged genes in the parquet, which the rank correlation of {rho} between the summed SHAP values and the out-of-fold scores supports (a shuffled copy gives {rho_shuffled}).
The table therefore says which columns the model uses, and it is not an out-of-fold measure of signal.

### One column alone

For each Open Targets column of the drugged-only model, the next table ranks the drugged genes by that column alone and reads the AUROC against the flag, with no model.
A column that holds one value only (`isCancerDriverGene` and `hasTEP`) is read as present against empty, and any other column leaves out the genes with no value.
A column with fewer than {min_genes} present or ranked genes, or with fewer than {min_class} flagged or unflagged ones, shows its counts only.

{table_columns}

The best of the {n_alone} columns is `{best_col}`, which reads AUROC {best_auc} on its own in its best direction, against {e1_auc} for the drugged-only model.
The best column is picked from {n_alone} by its point value, so its interval is a little optimistic.

### The DR+PU pipeline

The importance of the DR+PU model is saved only as pictures (`{dr_picture}` and two more), not as a table, and its model files and its SHAP file are not in the repo.
`README.md` (line {readme_feat_line}) states {dr_feat_total} columns in five groups ({dr_feat_groups}), which add up to that total.
This page cannot recount them, because the matrix is built from `backend/humanproof.db`, which is not in the repo.

## 6. The language model only arm: not run

The arm would ask {llm} about a gene alone, with no model of this repo, and would grade it only on outcomes that come after its release date, {release}.
That date is an input of this page and is not read from a file of this repo.
No such arm exists in this repo, and no model was called to make this page.
The table gives the reasons from the files.

{table_llm}

The label is a flag in one release of Open Targets, not an outcome that happens after a forecast day.
A forecast of it would test what a language model recalls about known safety targets, and the list above shows that these targets are public.
The list is a short fixed list written into the script, an illustration and not a sample.
An arm would need a later Open Targets release, forecasts frozen before that release is read, and a grade on the genes that the release newly flags.

## 7. The cells that were not run

Each row is a cell that a reader may expect and that this page did not fill, with its reason.
A file counts as absent when git does not track it, and every absent file named on this page is excluded by a rule of `.gitignore`, which keeps databases, pickles, compressed tables, model files and the outputs of the all-labeled model out of the repo.
A maintainer may therefore hold some of them locally, and this page does not see them.

{table_not_run}

## 8. Where the written docs, the code and the saved files differ

The table sets {n_claims} numbers that the written docs state about the DR+PU pipeline against what the saved files and the code give.
A stated number agrees when it lies within one and a half units of its last digit of the recomputed one.

{table_claims}

Of these {n_claims} numbers, {n_differ} differ.
The saved outcome model reads {dr3_auc_point} and {dr3_ap_point}, which are the numbers that the app page (`frontend/src/app/methodology/page.tsx`) states, and the saved summaries of the pseudo-outcome also match that page and not `Methodology.md`.
The code fixes the prior at {pi} (line {pi_line} of `data/safety_model_dr.py`) and the app page says so, while `Methodology.md` (line {pi_doc_line}) and `README.md` still describe the Elkan-Noto estimate, which `Methodology.md` writes as {pi_doc}.
The comment of the script (lines {shift_comment_lines} of `data/safety_model_dr.py`) and the app page (line {app_shift_line}) call the fixed prior a monotone shift that leaves AUROC and AUPRC unchanged, and the app page also says that it raised the out-of-fold AUROC from {shift_from} to {shift_to}.
Both cannot hold for a monotone shift, and in the code the prior is a sample weight of step 3 (line {weight_line}), so another prior refits the trees and can change the ranking.
The saved files hold the run with the fixed prior only, so they cannot show what the prior changed.
This page edits neither the scripts nor the written docs.
