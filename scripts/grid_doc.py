"""Make docs/GRID.md, the grid of experiments, results and features of HumanProof, from the saved files of the safety model.

    uv run --with pandas --with numpy --with pyarrow --with scikit-learn python scripts/grid_doc.py [--check]

The text is `docs/grid/template.md`, and every figure of the page is computed here from files of this repo:
the Open Targets parquet, the saved predictions, the saved SHAP files, the constants of the scripts under `data/` and a few lines of the written docs.
The numbers behind the page are written to `docs/grid/numbers.csv`.
`--check` fails when the committed page or CSV differs from what the files give, which is also a test (`tests/test_grid_doc.py`).
No model is retrained, and no network or paid service is used.
A claim of the text that the files do not hold stops the script: `values` asserts each one before its sentence is made.
A file is in the repo when git tracks it, so the page is the same in a fresh clone and in a checkout that also holds the files that `.gitignore` excludes (it needs a git checkout, and only reads from it).
"""
import argparse
import ast
import calendar
import csv
import datetime
import io
import json
import os
import re
import subprocess
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DOC = os.path.join(ROOT, 'docs', 'GRID.md')
TEMPLATE = os.path.join(ROOT, 'docs', 'grid', 'template.md')
NUMBERS = os.path.join(ROOT, 'docs', 'grid', 'numbers.csv')

SEED = 42                       # the seed of every bootstrap of this page
DRAWS = 2000                    # bootstrap draws over genes
CHUNK = 100                     # draws made at once
MIN_GENES = 30                  # a cell with fewer genes shows its counts only
MIN_CLASS = 10                  # a cell with fewer flagged or fewer unflagged genes shows its counts only
LIMIT = 380                     # the longest line of the section "At a glance"
RELEASE = '2026-04-24'          # the release date of the language model: an input of this page, not read from a file of the repo
LLM = 'DeepSeek V4 Flash'
PHASES = (0.25, 0.5, 0.75, 1.0)
ROMAN = {0.25: 'I', 0.5: 'II', 0.75: 'III', 1.0: 'IV'}      # as `Methodology.md` states (the page asserts that line)
TEXTBOOK = ('KCNH2', 'HTR2B', 'SCN5A', 'PTGS1', 'PTGS2', 'ADRB2', 'DRD2', 'HRH1')   # a fixed list of well known safety targets: an illustration, not a sample

OT_FILE = 'data/opentargets_target_prioritisation.parquet'     # paths of the repo are written with slashes, as git lists them
OUT = 'data/safety_model_output'
DOCS = {
    'README.md': ('README.md',),
    'Methodology.md': ('Methodology.md',),
    'CONTRIBUTING.md': ('CONTRIBUTING.md',),
    'frontend/src/app/methodology/page.tsx': ('frontend', 'src', 'app', 'methodology', 'page.tsx'),
}
APP_PAGE = 'frontend/src/app/methodology/page.tsx'
MODEL_A_FILES = tuple(f'{OUT}/{n}' for n in ('predictions.csv', 'shap_summary.csv', 'shap_values.npy', 'model.json'))
DRUGGED_MODEL_FILE = f'{OUT}/drugged_only/model.json'
DR_FILES = tuple(f'{OUT}/dr/{n}' for n in ('model_final.json', 'model_outcome.json', 'gene_shap_dr.json'))
DR_PICTURES = (f'{OUT}/dr/figures/04_feature_importance.png', 'figures/fig3_feature_importance.png', 'figures/fig8_feature_importance_full.png')
INPUT_FILES = ('backend/humanproof.db', 'data/LOEUF_scores.csv.gz', 'data/genebass_pLoF_filtered.pkl', 'data/cellxgene/celltype_log1p_mean_expression.csv.gz',
               'data/cellxgene/celltype_fraction_expressing.csv.gz', 'data/cellxgene/celltype_metadata.csv', 'data/cellxgene/gene_metadata.csv')
ABSENT = MODEL_A_FILES + (DRUGGED_MODEL_FILE,) + DR_FILES + INPUT_FILES      # the files that the page calls absent, each in a sentence or a count
DATA_EXPECTED = tuple(sorted(('data/demo/genes.json', 'data/demo/tissues.json', OT_FILE, f'{OUT}/dr/predictions.csv', f'{OUT}/drugged_only/predictions_oof.csv',
                              f'{OUT}/drugged_only/shap_summary.csv', f'{OUT}/drugged_only/shap_values.npy')))      # the data files of the repo: the page reads all of them but one small list
DATA_SUFFIXES = ('.parquet', '.csv', '.tsv', '.pkl', '.pickle', '.db', '.gz', '.npy', '.npz', '.feather', '.xlsx', '.h5ad', '.sqlite', '.json')
APART = ('docs/', 'scripts/', 'tests/', 'frontend/', '.claude/')       # this page and its checks, the app and the settings of the editor: neither data nor code of a model
GITIGNORE = ('*.db', '*.pkl', '*.csv.gz', 'data/cellxgene/*.csv', 'data/safety_model_output/**/*.json', 'data/safety_model_output/*.json', 'data/safety_model_output/*.csv',
             'data/safety_model_output/*.npy')                         # the rules of `.gitignore` that keep the absent files out
NAMING = re.compile(r'drug|trial|nct|chembl|indication|disease|cancer|oncolog|tumou?r|regimen|combination|sponsor|compound|efo', re.I)
EXPLAINED = ('maxClinicalTrialPhase', 'isCancerDriverGene', 'is_drugged')      # the columns that match NAMING: a phase, a driver flag and the drugged indicator, all of them properties of a gene

GROUP_WORD = {'Open Targets': 'Open Targets', 'LOEUF': 'LOEUF', 'pLoF': 'pLoF', 'Expression': 'expression'}
READING_LABEL = {
    'drugged_only': 'Drugged-only model, out of fold',
    'dr_step3': 'DR+PU step 3 outcome model, out of fold',
    'dr_step1': 'DR+PU step 1 propensity (entered a trial or not), out of fold',
    'dr_und_outcome': 'DR+PU outcome model on undrugged genes, which it never trained on',
    'dr_und_deployed': 'DR+PU deployed score on undrugged genes',
    'dr_deployed_drugged': 'DR+PU deployed score on drugged genes, in sample',
}


# ---------------------------------------------------------------------------------------------------------------------------------------------
# formats
# ---------------------------------------------------------------------------------------------------------------------------------------------
def num(x):
    return f'{int(round(float(x))):,}'


def f3(x):
    return f'{float(x):.3f}'


def f4(x):
    return f'{float(x):.4f}'


def pct(x):
    return f'{100 * float(x):.1f}%'


def interval(v, lo, hi):
    return f'{v:.3f} [{lo:.3f}, {hi:.3f}]'


def signed(v, lo, hi):
    return f'{v:+.3f} [{lo:+.3f}, {hi:+.3f}]'


def span(lo, hi):
    return num(lo) if lo == hi else f'{num(lo)} to {num(hi)}'


def joined(items):
    items = list(items)
    if not items:
        return 'none'
    return items[0] if len(items) == 1 else ', '.join(items[:-1]) + ' and ' + items[-1]


def md(head, rows, right=()):
    """A markdown table; a cell never holds a pipe, which would split it."""
    sep = ['---:' if i in right else '---' for i in range(len(head))]
    lines = ['| ' + ' | '.join(head) + ' |', '|' + '|'.join(sep) + '|']
    for r in rows:
        cells = [str(c) for c in r]
        assert len(cells) == len(head), (head, cells)
        assert not any('|' in c for c in cells), cells
        lines.append('| ' + ' | '.join(cells) + ' |')
    return '\n'.join(lines)


def read_text(path):
    with open(path, encoding='utf-8') as f:
        return f.read()


# ---------------------------------------------------------------------------------------------------------------------------------------------
# the bootstrap over genes: the AUROC and the average precision of every draw, made at once from the weights of the draw
# ---------------------------------------------------------------------------------------------------------------------------------------------
def _prepare(y, s):
    y = np.asarray(y).astype(np.int64)
    s = np.asarray(s, dtype=float)
    assert len(y) == len(s) and not np.isnan(s).any()
    order = np.argsort(s, kind='mergesort')
    ys, ss = y[order], s[order]
    starts = np.flatnonzero(np.r_[True, ss[1:] != ss[:-1]])        # a block of tied scores is one group
    return order, ys, starts


def _metrics(ys, starts, w):
    """AUROC and average precision of each row of `w` (draws by genes, in the order of ascending score); a draw with no flagged or no unflagged gene gives nan."""
    pos, neg = w * ys, w * (1 - ys)
    gp, gn = np.add.reduceat(pos, starts, axis=1), np.add.reduceat(neg, starts, axis=1)
    n_pos, n_neg = gp.sum(axis=1), gn.sum(axis=1)
    below = np.cumsum(gn, axis=1) - gn
    gp_d, gn_d = gp[:, ::-1], gn[:, ::-1]
    tp, fp = np.cumsum(gp_d, axis=1), np.cumsum(gn_d, axis=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        auc = (gp * (below + 0.5 * gn)).sum(axis=1) / (n_pos * n_neg)
        precision = np.where(tp + fp > 0, tp / (tp + fp), 0.0)
        ap = (gp_d / n_pos[:, None] * precision).sum(axis=1)
    bad = (n_pos == 0) | (n_neg == 0)
    auc[bad], ap[bad] = np.nan, np.nan
    return auc, ap


def _draws(n, draws, seed):
    rng = np.random.default_rng(seed)
    for k in range(0, draws, CHUNK):
        yield rng.multinomial(n, np.full(n, 1.0 / n), size=min(CHUNK, draws - k))      # how often each gene enters the draw


def _percentiles(values):
    values = values[~np.isnan(values)]
    return tuple(float(x) for x in np.percentile(values, [2.5, 97.5]))


def point(y, s):
    """The AUROC and the average precision of the genes as they are."""
    order, ys, starts = _prepare(y, s)
    a, p = _metrics(ys, starts, np.ones((1, len(ys))))
    return float(a[0]), float(p[0])


def bootstrap(y, s, draws=DRAWS, seed=SEED):
    """The AUROC and the average precision with 95% percentile intervals over genes."""
    order, ys, starts = _prepare(y, s)
    a0, p0 = point(y, s)
    auc, ap = [], []
    for counts in _draws(len(ys), draws, seed):
        a, p = _metrics(ys, starts, counts[:, order].astype(float))
        auc.append(a)
        ap.append(p)
    auc, ap = np.concatenate(auc), np.concatenate(ap)
    lo_a, hi_a = _percentiles(auc)
    lo_p, hi_p = _percentiles(ap)
    return dict(auroc=a0, auroc_lo=lo_a, auroc_hi=hi_a, auprc=p0, auprc_lo=lo_p, auprc_hi=hi_p, draws_used=int((~np.isnan(auc)).sum()))


def bootstrap_pair(y, s1, s2, draws=DRAWS, seed=SEED):
    """The difference s2 minus s1 of the AUROC and of the average precision, on the same draws, with 95% percentile intervals."""
    o1, y1, st1 = _prepare(y, s1)
    o2, y2, st2 = _prepare(y, s2)
    a1, p1 = point(y, s1)
    a2, p2 = point(y, s2)
    da, dp = [], []
    for counts in _draws(len(y1), draws, seed):
        x1, q1 = _metrics(y1, st1, counts[:, o1].astype(float))
        x2, q2 = _metrics(y2, st2, counts[:, o2].astype(float))
        da.append(x2 - x1)
        dp.append(q2 - q1)
    da, dp = np.concatenate(da), np.concatenate(dp)
    lo_a, hi_a = _percentiles(da)
    lo_p, hi_p = _percentiles(dp)
    return dict(auroc_1=a1, auroc_2=a2, auroc_diff=a2 - a1, auroc_diff_lo=lo_a, auroc_diff_hi=hi_a,
                auprc_1=p1, auprc_2=p2, auprc_diff=p2 - p1, auprc_diff_lo=lo_p, auprc_diff_hi=hi_p, draws_used=int((~np.isnan(da)).sum()))


def fold_sizes(n_pos, n_neg, n_splits, seed):
    """The sizes of the folds of `StratifiedKFold(shuffle=True)`, which depend on the class counts only (the folds themselves are not stored)."""
    y = np.r_[np.ones(n_pos, dtype=int), np.zeros(n_neg, dtype=int)]
    sizes = [(len(tr), len(va)) for tr, va in StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed).split(np.zeros((len(y), 1)), y)]
    train, score = [a for a, b in sizes], [b for a, b in sizes]
    assert all(a + b == len(y) for a, b in sizes) and sum(score) == len(y)
    return min(train), max(train), min(score), max(score)


# ---------------------------------------------------------------------------------------------------------------------------------------------
# facts of the code, read from the syntax tree of the scripts, so that no constant is typed here
# ---------------------------------------------------------------------------------------------------------------------------------------------
UNKNOWN = object()


def evaluate(node, env):
    """The value of a literal, a name of the module, a list, a dict or a call of dict(); anything else is UNKNOWN."""
    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, ast.Name):
        return env.get(node.id, UNKNOWN)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        x = evaluate(node.operand, env)
        return UNKNOWN if x is UNKNOWN else -x
    if isinstance(node, (ast.List, ast.Tuple)):
        items = [evaluate(e, env) for e in node.elts]
        return UNKNOWN if any(i is UNKNOWN for i in items) else items
    if isinstance(node, ast.Dict):
        out = {}
        for k, v in zip(node.keys, node.values):
            key = evaluate(k, env) if k is not None else UNKNOWN
            if key is UNKNOWN:
                return UNKNOWN
            out[key] = evaluate(v, env)
        return out
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == 'dict':
        out = {}
        for kw in node.keywords:
            val = evaluate(kw.value, env)
            if kw.arg is None:
                if isinstance(val, dict):
                    out.update(val)
            else:
                out[kw.arg] = val
        return out
    return UNKNOWN


def module_env(tree):
    env = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            name, value = node.targets[0].id, node.value
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name) and node.value is not None:
            name, value = node.target.id, node.value
        else:
            continue
        got = evaluate(value, env)
        if got is not UNKNOWN:
            env[name] = got
    return env


def line_of_assign(tree, name):
    for node in tree.body:
        target = node.targets[0] if isinstance(node, ast.Assign) and len(node.targets) == 1 else getattr(node, 'target', None)
        if isinstance(target, ast.Name) and target.id == name:
            return node.lineno
    raise KeyError(name)


def function_node(tree, name):
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return node
    raise KeyError(name)


def function_var(tree, func, var, env):
    for node in ast.walk(function_node(tree, func)):
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name) and node.targets[0].id == var:
            return evaluate(node.value, env)
    raise KeyError(var)


def function_call(tree, func, callee, env):
    for node in ast.walk(function_node(tree, func)):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == callee:
            return {kw.arg: evaluate(kw.value, env) for kw in node.keywords}
    raise KeyError(callee)


def keyword_line(tree, func, attr, keyword):
    """The line of the keyword `keyword` of the call `.attr(...)` inside the function `func`."""
    for node in ast.walk(function_node(tree, func)):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == attr:
            for kw in node.keywords:
                if kw.arg == keyword:
                    return kw.value.lineno
    raise KeyError((func, attr, keyword))


def code_facts(texts):
    sm, dr = ast.parse(texts['safety_model.py']), ast.parse(texts['safety_model_dr.py'])
    env_sm, env_dr = module_env(sm), module_env(dr)
    lines_sm = texts['safety_model.py'].splitlines()
    symbol_line = [i + 1 for i, line in enumerate(lines_sm) if 'df["targetId"].map(ENSG_TO_SYMBOL)' in line]
    assert len(symbol_line) == 1
    lines_dr = texts['safety_model_dr.py'].splitlines()
    clip = [(i + 1, m) for i, line in enumerate(lines_dr) if (m := re.search(r'np\.clip\(final_model\.predict\([^)]*\), *([\d.]+), *([\d.]+)\)', line))]
    assert len(clip) == 1, 'the deployed score is no longer clipped in one line'
    assert '**{k: v for k, v in _XGB_COMMON.items() if k != "n_estimators"}' in texts['safety_model_dr.py'], 'the final regressor no longer shares the settings of the classifiers'
    shift = [i + 1 for i, line in enumerate(lines_dr) if 'monotone logit-shift' in line]
    weight = [i + 1 for i, line in enumerate(lines_dr) if 'pu_w = np.where(Y == 1, pi_p, 1.0 - pi_p)' in line]
    assert len(shift) == 1 and 'AUROC and AUPRC exactly' in lines_dr[shift[0]], 'the comment on the prior no longer calls it a monotone shift that keeps AUROC and AUPRC'
    assert len(weight) == 1 and 'sample_weight=sw' in texts['safety_model_dr.py'], 'the prior is no longer a sample weight of step 3'
    return dict(
        genes=env_sm['GENES'], ot_features=env_sm['OT_FEATURES'], loeuf_features=env_sm['LOEUF_FEATURES'], plof_base=env_sm['PLOF_BASE_FEATURES'],
        genes_line=line_of_assign(sm, 'GENES'), symbol_line=symbol_line[0], clip_line=clip[0][0], clip_range=(float(clip[0][1][1]), float(clip[0][1][2])),
        shift_comment_lines=f'{shift[0]} and {shift[0] + 1}', weight_line=weight[0],
        pi_override=env_dr['PI_OVERRIDE'], pi_line=line_of_assign(dr, 'PI_OVERRIDE'), n_folds=env_dr['N_FOLDS'], seed=env_dr['SEED'], max_ipw=env_dr['MAX_IPW'],
        e1_params=function_var(sm, 'train_drugged_only', 'model_params', env_sm), e1_cv=function_call(sm, 'train_drugged_only', 'StratifiedKFold', env_sm),
        e1_stop_line=keyword_line(sm, 'train_drugged_only', 'fit', 'eval_set'),
        e2_params=function_var(sm, 'train_and_evaluate', 'model_params', env_sm), e2_cv=function_call(sm, 'train_and_evaluate', 'StratifiedKFold', env_sm),
        dr_clf=env_dr['_XGB_CLF'], dr_reg=env_dr['_XGB_REG'], dr_cv=function_call(dr, 'fit_outcome_crossfit', 'StratifiedKFold', env_dr),
        dr_prop_cv=function_call(dr, 'fit_propensity_crossfit', 'StratifiedKFold', env_dr),
        dr_stop_line=keyword_line(dr, 'fit_outcome_crossfit', 'fit', 'eval_set'),
    )


LEARNER = re.compile(r'\b(XGB[A-Za-z]+|Logistic[A-Za-z]+|RandomForest[A-Za-z]+|Dummy[A-Za-z]+|GradientBoosting[A-Za-z]+|HistGradient[A-Za-z]+|SVC|Ridge[A-Za-z]*|Lasso[A-Za-z]*|ElasticNet)\(')
SPLITTER = re.compile(r'\b(StratifiedKFold|KFold|GroupKFold|StratifiedGroupKFold|TimeSeriesSplit|train_test_split)\(')


def names_in(texts, pattern):
    """The distinct names that a pattern finds in the texts of the scripts."""
    return set().union(*[set(pattern.findall(t)) for t in texts.values()])


# ---------------------------------------------------------------------------------------------------------------------------------------------
# the sources
# ---------------------------------------------------------------------------------------------------------------------------------------------
def run_git(root, *args, stdin=None):
    """The return code and the output of a git command that only reads the repo."""
    try:
        done = subprocess.run(['git', *args], cwd=root, input=stdin, capture_output=True, text=True)
    except OSError as e:
        raise RuntimeError(f'git is needed to read which files the repo tracks: {e}') from e
    return done.returncode, done.stdout, done.stderr


def tracked_files(root):
    """The files that git tracks under the root, written with slashes."""
    code, out, err = run_git(root, 'ls-files', '-z')
    if code != 0:
        raise RuntimeError(f'`git ls-files` failed in {root}, and this page needs a git checkout: {err.strip()}')
    return {p for p in out.split('\0') if p}


def ignored_among(root, paths):
    """The paths that a rule of `.gitignore` excludes, whether or not the file is there."""
    if not paths:
        return set()
    code, out, err = run_git(root, 'check-ignore', '--stdin', stdin='\n'.join(paths) + '\n')
    if code not in (0, 1):
        raise RuntimeError(f'`git check-ignore` failed in {root}: {err.strip()}')
    return set(out.split('\n')) - {''}


class Sources:
    """The files of the repo that the page reads, and nothing else."""

    def __init__(self, root=ROOT):
        self.root = root
        join = lambda *p: os.path.join(root, *'/'.join(p).split('/'))              # noqa: E731
        self.ot = pd.read_parquet(join(OT_FILE))
        self.oof = pd.read_csv(join(OUT, 'drugged_only', 'predictions_oof.csv'))
        self.dr = pd.read_csv(join(OUT, 'dr', 'predictions.csv'))
        self.shap = pd.read_csv(join(OUT, 'drugged_only', 'shap_summary.csv'))
        self.shap_values = np.load(join(OUT, 'drugged_only', 'shap_values.npy'), allow_pickle=False)
        with open(join('data', 'demo', 'genes.json'), encoding='utf-8') as f:
            self.gene_list = json.load(f)
        with open(join('data', 'demo', 'tissues.json'), encoding='utf-8') as f:
            self.tissues = json.load(f)
        self.lines = {name: read_text(join(*path)).splitlines() for name, path in DOCS.items()}
        self.tracked = tracked_files(root)
        self.gitignore = [x.strip() for x in read_text(join('.gitignore')).splitlines()]
        scripts = {p: read_text(join(p)) for p in sorted(self.tracked) if p.endswith('.py') and not p.startswith(APART)}
        self.scripts = scripts
        self.facts = code_facts({n: scripts[f'data/{n}'] for n in ('safety_model.py', 'safety_model_dr.py')})
        self.all_learners = names_in(scripts, LEARNER)
        self.all_splitters = names_in(scripts, SPLITTER)
        self.data_files = sorted(p for p in self.tracked if p.endswith(DATA_SUFFIXES) and not p.startswith(APART))
        paths = ABSENT + DR_PICTURES
        self.present = {p: p in self.tracked for p in paths}
        self.ignored = ignored_among(root, [p for p in paths if not self.present[p]])

    def grab(self, file, pattern):
        """The one line of a written doc that matches the pattern: its number and its match."""
        hits = [(i + 1, m) for i, line in enumerate(self.lines[file]) if (m := re.search(pattern, line))]
        assert len(hits) == 1, f'{file}: {len(hits)} lines match {pattern!r}'
        return hits[0]


# ---------------------------------------------------------------------------------------------------------------------------------------------
# the numbers
# ---------------------------------------------------------------------------------------------------------------------------------------------
def flag_of(df):
    return df['hasSafetyEvent'] == -1.0


def drugged_of(df):
    return df['maxClinicalTrialPhase'].notna()


def verify(src):
    """The files agree with one another, so that every count below is a count of the same genes."""
    ot, oof, dr, k = src.ot, src.oof, src.dr, src.facts
    assert ot.shape[1] == 1 + len(k['ot_features']) + 2 and set(ot.columns) == {'targetId', 'hasSafetyEvent', 'maxClinicalTrialPhase'} | set(k['ot_features'])
    assert set(ot['hasSafetyEvent'].dropna().unique()) == {-1.0}, 'the flag holds a value other than -1'
    assert set(ot['maxClinicalTrialPhase'].dropna().unique()) == set(PHASES)
    assert ot['targetId'].is_unique
    assert set(oof['targetId']) == set(ot.loc[drugged_of(ot), 'targetId']) and oof['targetId'].is_unique
    merged = oof.merge(ot, on='targetId', suffixes=('', '_ot'))
    assert (merged['true_label'] == (merged['hasSafetyEvent_ot'] == -1.0).astype(int)).all()
    assert (merged['maxClinicalTrialPhase'] == merged['maxClinicalTrialPhase_ot']).all()
    assert dr['targetId'].is_unique and set(dr['targetId']) <= set(ot['targetId'])
    merged = dr.merge(ot, on='targetId', suffixes=('', '_ot'))
    assert (merged['maxClinicalTrialPhase'].fillna(-1) == merged['maxClinicalTrialPhase_ot'].fillna(-1)).all()
    assert (merged['hasSafetyEvent'].fillna(0) == merged['hasSafetyEvent_ot'].fillna(0)).all()
    assert (dr['is_drugged'] == drugged_of(dr).astype(int)).all()
    expected = np.where(flag_of(dr), 'positive', np.where(drugged_of(dr), 'drugged_safe', 'unlabeled'))
    assert (dr['safety_label'].values == expected).all(), 'the saved labels of the DR+PU file are not the rule of the code'
    assert src.shap_values.shape == (len(oof), len(src.shap)) and src.shap['mean_abs_shap'].is_monotonic_decreasing
    assert np.allclose(np.sort(np.abs(src.shap_values).mean(axis=0))[::-1], src.shap['mean_abs_shap'].values, atol=1e-6), 'the SHAP summary is not the SHAP matrix'
    assert src.all_learners == {'XGBClassifier', 'XGBRegressor'}, f'the scripts hold other learners: {src.all_learners}'
    assert src.all_splitters == {'StratifiedKFold'}, f'the scripts hold other splitters: {src.all_splitters}'
    # the page says that a model is compared with chance only: the figures draw the chance lines, and no other learner exists
    assert any('Random classifier' in t for t in src.scripts.values()) and any('Random baseline' in t for t in src.scripts.values()), 'no script draws the chance lines that the page names'
    # the page says that the files name no drug and no trial: the data files of the repo are these, and no column of a table names one
    assert src.data_files == list(DATA_EXPECTED), f'the repo tracks a data file that this page does not know, or lacks one: {sorted(set(src.data_files) ^ set(DATA_EXPECTED))}'
    named = sorted({c for df in (ot, oof, dr) for c in df.columns if NAMING.search(c)} - set(EXPLAINED))
    assert not named, f'a column names a drug, a trial or a disease: {named}'
    assert isinstance(src.gene_list, list) and all(isinstance(g, str) for g in src.gene_list), 'the gene list is not a list of symbols'
    assert isinstance(src.tissues, dict) and all(isinstance(v, dict) and set(v) == {'organ', 'cell_types'} for v in src.tissues.values()), 'the tissue list is not organs with cell types'


def label_row(name, df):
    flag, drugged = flag_of(df), drugged_of(df)
    return dict(universe=name, genes=len(df), positive=int(flag.sum()), positive_drugged=int((flag & drugged).sum()), positive_no_phase=int((flag & ~drugged).sum()),
                drugged_safe=int((~flag & drugged).sum()), unlabeled=int((~flag & ~drugged).sum()))


def stage_rows(name, df):
    flag, drugged = flag_of(df), drugged_of(df)
    rows = []
    for code in PHASES + (None,):
        sel = df['maxClinicalTrialPhase'] == code if code is not None else ~drugged
        n = int(sel.sum())
        rows.append(dict(universe=name, phase='none' if code is None else str(code), genes=n, positive=int((flag & sel).sum()), drugged_safe=int((~flag & drugged & sel).sum()),
                         unlabeled=int((~flag & ~drugged & sel).sum()), flag_rate=float((flag & sel).sum() / n)))
    n = len(df)
    rows.append(dict(universe=name, phase='all', genes=n, positive=int(flag.sum()), drugged_safe=int((~flag & drugged).sum()), unlabeled=int((~flag & ~drugged).sum()),
                     flag_rate=float(flag.sum() / n)))
    assert sum(r['genes'] for r in rows[:-1]) == n
    return rows


def read_result(reading, stage, y, s):
    """One reading: the counts, and the metrics with their intervals when the cell is large enough."""
    y = np.asarray(y).astype(int)
    n, e = len(y), int(y.sum())
    row = dict(reading=reading, stage=stage, genes=n, events=e, base_rate=e / n, shown=bool(n >= MIN_GENES and e >= MIN_CLASS and n - e >= MIN_CLASS))
    keys = ('auroc', 'auroc_lo', 'auroc_hi', 'auprc', 'auprc_lo', 'auprc_hi', 'draws_used')
    row.update(dict.fromkeys(keys, np.nan))
    if row['shown']:
        row.update(bootstrap(y, s))
    return row


def compute_counts(src, F, K):
    ot, oof, dr, k = src.ot, src.oof, src.dr, src.facts
    drugged_ot = ot[drugged_of(ot)]
    F['labels'] = pd.DataFrame([label_row('Open Targets file', ot), label_row('Drugged-only rows', drugged_ot), label_row('DR+PU file', dr)])
    lab = F['labels'].set_index('universe')
    assert (lab.loc['Open Targets file', ['positive', 'drugged_safe', 'unlabeled']].sum() == len(ot))
    saved = dr['safety_label'].value_counts()
    assert (saved['positive'], saved['drugged_safe'], saved['unlabeled']) == tuple(lab.loc['DR+PU file', ['positive', 'drugged_safe', 'unlabeled']])
    K.update(ot_n=len(ot), ot_flag=int(lab.loc['Open Targets file', 'positive']), ot_safe=int(lab.loc['Open Targets file', 'drugged_safe']), ot_unl=int(lab.loc['Open Targets file', 'unlabeled']),
             ot_no_phase_flag=int(lab.loc['Open Targets file', 'positive_no_phase']), ot_flag_drugged=int(lab.loc['Open Targets file', 'positive_drugged']),
             e1_n=len(oof), e1_events=int(oof['true_label'].sum()), dr_n=len(dr), dr_drugged=int(dr['is_drugged'].sum()), dr_und=int((dr['is_drugged'] == 0).sum()),
             dr_flag=int(lab.loc['DR+PU file', 'positive']), dr_flag_drugged=int(lab.loc['DR+PU file', 'positive_drugged']), dr_flag_no_phase=int(lab.loc['DR+PU file', 'positive_no_phase']),
             dr_safe=int(lab.loc['DR+PU file', 'drugged_safe']), dr_unl=int(lab.loc['DR+PU file', 'unlabeled']))
    assert K['e1_n'] == int(lab.loc['Drugged-only rows', 'genes']) and K['e1_events'] == int(lab.loc['Drugged-only rows', 'positive'])
    assert int(lab.loc['Drugged-only rows', 'unlabeled']) == 0 and int(lab.loc['Drugged-only rows', 'positive_no_phase']) == 0
    # the genes that each model trains on and scores
    labeled_a = K['ot_flag'] + K['ot_safe']
    spec = [
        ('drugged_only', K['ot_n'], K['e1_n'], K['e1_n'], K['e1_events'], k['e1_cv']),
        ('all_labeled', K['ot_n'], labeled_a, K['ot_n'], K['ot_flag'], k['e2_cv']),
        ('dr_step1', K['dr_n'], K['dr_n'], K['dr_n'], K['dr_drugged'], k['dr_prop_cv']),
        ('dr_step3', K['dr_n'], K['dr_drugged'], K['dr_drugged'], K['dr_flag_drugged'], k['dr_cv']),
    ]
    rows = []
    for name, universe, trained, scored, events, cv in spec:
        assert cv['shuffle'] is True
        tr_lo, tr_hi, sc_lo, sc_hi = fold_sizes(events, trained - events, cv['n_splits'], cv['random_state'])
        rows.append(dict(reading=name, universe=universe, trained=trained, scored=scored, fold_train_min=tr_lo, fold_train_max=tr_hi, fold_score_min=sc_lo, fold_score_max=sc_hi,
                         events=events, base_rate=events / trained, folds=cv['n_splits'], seed=cv['random_state']))
    rows.append(dict(reading='dr_step5', universe=K['dr_n'], trained=K['dr_n'], scored=K['dr_n'], fold_train_min=np.nan, fold_train_max=np.nan, fold_score_min=np.nan, fold_score_max=np.nan,
                     events=np.nan, base_rate=np.nan, folds=np.nan, seed=np.nan))
    F['setup'] = pd.DataFrame(rows)
    assert len({(r['folds'], r['seed']) for r in rows[:4]}) == 1, 'the models do not share the fold rule'
    K['folds'], K['seed_code'] = int(rows[0]['folds']), int(rows[0]['seed'])
    assert K['seed_code'] == SEED
    # one row per gene
    ids = lambda df: (len(df), int(df['targetId'].nunique()))                         # noqa: E731
    dedup = []
    for name, df, sym in (('Open Targets parquet', ot, None), ('Drugged-only predictions', oof, oof['gene_symbol'].dropna()), ('DR+PU predictions', dr, dr['gene_symbol'].dropna())):
        rows_n, uid = ids(df)
        dedup.append(dict(file=name, rows=rows_n, unique_ids=uid, duplicate_ids=rows_n - uid, rows_with_symbol=np.nan if sym is None else len(sym),
                          unique_symbols=np.nan if sym is None else int(sym.nunique()), duplicate_symbols=np.nan if sym is None else len(sym) - int(sym.nunique())))
    genes = pd.Series(src.gene_list)
    dedup.append(dict(file='data/demo/genes.json', rows=len(genes), unique_ids=np.nan, duplicate_ids=np.nan, rows_with_symbol=len(genes), unique_symbols=int(genes.nunique()),
                      duplicate_symbols=len(genes) - int(genes.nunique())))
    F['dedup'] = pd.DataFrame(dedup)
    assert (F['dedup']['duplicate_ids'].dropna() == 0).all() and (F['dedup']['duplicate_symbols'].dropna() == 0).all(), 'a file holds a duplicate gene'
    K.update(gene_list_n=len(genes), dr_symbols_in_list=int(dr['gene_symbol'].isin(set(genes)).sum()), e1_sym=int(oof['gene_symbol'].notna().sum()), e1_nosym=int(oof['gene_symbol'].isna().sum()))
    # stage and cancer driver
    F['stage'] = pd.DataFrame(stage_rows('Open Targets file', ot) + stage_rows('DR+PU file', dr))
    driver = ot['isCancerDriverGene'] == -1.0
    d = ot[drugged_of(ot)]
    rows = []
    for name, sel in (('flagged as driver', driver[d.index]), ('not flagged as driver', ~driver[d.index])):
        sub = d[sel]
        e = int(flag_of(sub).sum())
        rows.append(dict(group=name, genes=len(sub), flagged=e, not_flagged=len(sub) - e, flag_rate=e / len(sub)))
    e = int(flag_of(d).sum())
    rows.append(dict(group='all drugged genes', genes=len(d), flagged=e, not_flagged=len(d) - e, flag_rate=e / len(d)))
    F['driver'] = pd.DataFrame(rows)
    K.update(drivers_all=int(driver.sum()), drivers_undrugged=int((driver & ~drugged_of(ot)).sum()), drivers_drugged=int(rows[0]['genes']))
    assert rows[0]['genes'] + rows[1]['genes'] == len(d)


def compute_results(src, F, K):
    oof, dr = src.oof, src.dr
    d, u = dr[dr['is_drugged'] == 1], dr[dr['is_drugged'] == 0]
    y_d, y_u = flag_of(d).astype(int).values, flag_of(u).astype(int).values
    readings = {
        'drugged_only': (oof['true_label'].values, oof['safety_score_oof'].values, oof['maxClinicalTrialPhase'].values),
        'dr_step3': (y_d, d['m_hat'].values, d['maxClinicalTrialPhase'].values),
        'dr_step1': (dr['is_drugged'].values, dr['pi_hat'].values, None),
        'dr_und_outcome': (y_u, u['m_hat'].values, None),
        'dr_und_deployed': (y_u, u['safety_score_dr'].values, None),
    }
    rows = [read_result(name, 'all', y, s) for name, (y, s, ph) in readings.items()]
    for name in ('drugged_only', 'dr_step3'):
        y, s, ph = readings[name]
        for code in PHASES:
            sel = ph == code
            rows.append(read_result(name, str(code), y[sel], s[sel]))
    # step 4 gives an undrugged gene the score of the outcome model as its target, and the deployed score stays inside the clip of the code
    assert (u['Y_tilde'] == u['m_hat']).all(), 'step 4 no longer gives an undrugged gene the score of the outcome model'
    lo, hi = src.facts['clip_range']
    assert dr['safety_score_dr'].min() == lo and dr['safety_score_dr'].max() <= hi, 'the saved deployed score leaves the clip of the code'
    # the deployed score on the drugged genes is fitted to these labels: its in-sample values are shown for the record only
    a, p = point(y_d, d['safety_score_dr'].values)
    rows.append(dict(reading='dr_deployed_drugged', stage='all', genes=len(y_d), events=int(y_d.sum()), base_rate=float(y_d.mean()), shown=False, auroc=a, auroc_lo=np.nan, auroc_hi=np.nan,
                     auprc=p, auprc_lo=np.nan, auprc_hi=np.nan, draws_used=np.nan))
    F['results'] = pd.DataFrame(rows)
    res = F['results'].set_index(['reading', 'stage'])
    assert res.loc[('drugged_only', 'all'), 'genes'] == K['e1_n'] and res.loc[('dr_step3', 'all'), 'genes'] == K['dr_drugged']
    assert res.loc[('dr_step3', 'all'), 'events'] == K['dr_flag_drugged'] and res.loc[('dr_und_outcome', 'all'), 'events'] == K['dr_flag_no_phase']
    assert (res['draws_used'].dropna() >= DRAWS - 5).all(), 'many draws have no flagged or no unflagged gene'
    # the in-sample claim: every flagged drugged gene scores above every unflagged one
    s_d = d['safety_score_dr'].values
    K.update(min_flag_score=float(s_d[y_d == 1].min()), max_safe_score=float(s_d[y_d == 0].max()),
             und_zero=int((u['safety_score_dr'] == lo).sum()), und_flag_zero=int(((u['safety_score_dr'] == lo) & (y_u == 1)).sum()),
             und_n=len(u), und_pos=int(y_u.sum()), und_neg=int((y_u == 0).sum()), deployed_min=float(dr['safety_score_dr'].min()), deployed_max=float(dr['safety_score_dr'].max()))
    # the paired difference on the genes that both files score
    m = d[['targetId', 'm_hat']].merge(oof[['targetId', 'safety_score_oof', 'true_label']], on='targetId', how='inner')
    assert len(m) == K['dr_drugged'] and set(m['true_label']) == {0, 1}
    pair = bootstrap_pair(m['true_label'].values, m['safety_score_oof'].values, m['m_hat'].values)
    pair.update(genes=len(m), events=int(m['true_label'].sum()), base_rate=float(m['true_label'].mean()))
    F['paired'] = pd.DataFrame([pair])
    K['shared_n'] = len(m)
    K['not_in_dr'] = len(oof) - len(m)
    K['not_in_dr_flag'] = int(oof['true_label'].sum() - m['true_label'].sum())


def feature_group(name, k):
    if name in k['ot_features']:
        return 'Open Targets'
    if name in k['loeuf_features']:
        return 'LOEUF'
    if name.startswith('plof_'):
        return 'pLoF'
    if name.startswith('expr_'):
        return 'Expression'
    raise ValueError(name)


def column_alone(v, flag):
    """Rank the drugged genes by one column and read the AUROC against the flag, or give the counts only."""
    indicator = v.dropna().nunique() == 1
    if indicator:
        s, y = v.notna().astype(float).values, flag.values
        used, present = len(y), int(v.notna().sum())
        basis = 'present against empty'
        sel = v.notna().values
        big = present >= MIN_GENES and int(y[sel].sum()) >= MIN_CLASS and int((1 - y[sel]).sum()) >= MIN_CLASS
    else:
        keep = v.notna().values
        s, y = v.values[keep], flag.values[keep]
        used, present = int(keep.sum()), np.nan
        basis = 'genes with a value'
        big = used >= MIN_GENES and int(y.sum()) >= MIN_CLASS and int((1 - y).sum()) >= MIN_CLASS
    out = dict(basis=basis, genes=used, present=present, flagged=int(np.sum(y)), shown=bool(big))
    keys = ('auroc', 'auroc_lo', 'auroc_hi')
    out.update(dict.fromkeys(keys, np.nan))
    if big:
        b = bootstrap(y, s)
        out.update(auroc=b['auroc'], auroc_lo=b['auroc_lo'], auroc_hi=b['auroc_hi'])
    return out


def compute_features(src, F, K):
    ot, oof, k = src.ot, src.oof, src.facts
    sh = src.shap.copy()
    sh['rank'] = np.arange(1, len(sh) + 1)
    sh['group'] = sh['feature'].map(lambda f: feature_group(f, k))
    sh['share'] = sh['mean_abs_shap'] / sh['mean_abs_shap'].sum()
    F['shap'] = sh[['rank', 'feature', 'group', 'mean_abs_shap', 'share']]
    rows = []
    for g in ('Open Targets', 'LOEUF', 'pLoF', 'Expression'):
        sub = sh[sh['group'] == g]
        if g == 'pLoF':
            base = int(sub['feature'].isin(k['plof_base']).sum())
            detail = f'{base} global and {len(sub) - base} by organ system'
        elif g == 'Expression':
            organ = int(sub['feature'].str.startswith('expr_organ_').sum())
            pct_cols = int(sub['feature'].str.endswith('_pct').sum())
            detail = f'{len(sub) - organ} cell type ({len(sub) - organ - pct_cols} mean, {pct_cols} percent expressing) and {organ} by organ'
        else:
            detail = f'{len(sub)} column' + ('' if len(sub) == 1 else 's')
        top = sub.iloc[0]
        rows.append(dict(group=g, columns=len(sub), detail=detail, mean_abs_shap=float(sub['mean_abs_shap'].sum()), share=float(sub['share'].sum()), top_feature=top['feature'], top_value=float(top['mean_abs_shap'])))
    F['groups'] = pd.DataFrame(rows)
    assert F['groups']['columns'].sum() == len(sh) and abs(F['groups']['share'].sum() - 1) < 1e-9
    g = F['groups'].set_index('group')
    expr = sh[sh['group'] == 'Expression']['feature']
    K.update(n_saved=len(sh), ot_share=float(g.loc['Open Targets', 'share']), n_ot=int(g.loc['Open Targets', 'columns']), n_nonot=len(sh) - int(g.loc['Open Targets', 'columns']),
             n_expr=len(expr), n_expr_old=int((~(expr.str.startswith('expr_ct_') | expr.str.startswith('expr_organ_'))).sum()), shap_rows=src.shap_values.shape[0], shap_cols=src.shap_values.shape[1],
             genes_list_n=len(k['genes']))
    symbols = set(oof['gene_symbol'].dropna())
    listed_drugged = {s for s, e in k['genes'].items() if e in set(oof['targetId'])}
    assert symbols == listed_drugged, 'the symbols of the drugged-only file are not the drugged genes of the fixed list'
    K['n_sym_drugged'] = len(listed_drugged)
    assert K['e1_sym'] == K['n_sym_drugged']
    # the SHAP matrix follows the order of the drugged genes of the parquet: checked against the out-of-fold scores
    drugged_ot = ot[drugged_of(ot)].reset_index(drop=True)
    rowsum = src.shap_values.sum(axis=1)
    aligned = drugged_ot['targetId'].map(oof.set_index('targetId')['safety_score_oof']).values
    K['shap_rho'] = float(pd.Series(rowsum).corr(pd.Series(aligned), method='spearman'))
    K['shap_rho_shuffled'] = float(pd.Series(np.random.default_rng(SEED).permutation(rowsum)).corr(pd.Series(aligned), method='spearman'))
    K['shap_insample_auc'] = float(roc_auc_score(flag_of(drugged_ot).astype(int).values, rowsum))
    # one Open Targets column alone
    d = ot[drugged_of(ot)]
    flag = flag_of(d).astype(int)
    rank = dict(zip(sh['feature'], sh['rank']))
    value = dict(zip(sh['feature'], sh['mean_abs_shap']))
    rows = []
    for col in sorted(k['ot_features'], key=lambda c: rank[c]):
        r = dict(column=col, shap_rank=rank[col], mean_abs_shap=value[col])
        r.update(column_alone(d[col].reset_index(drop=True), flag.reset_index(drop=True)))
        rows.append(r)
    F['columns'] = pd.DataFrame(rows)
    shown = F['columns'][F['columns']['shown']].copy()
    shown['directed'] = shown['auroc'].where(shown['auroc'] >= 0.5, 1 - shown['auroc'])
    best = shown.loc[(shown['auroc'] - 0.5).abs().idxmax()]
    flipped = best['auroc'] < 0.5
    K.update(best_col=str(best['column']), best_flipped=bool(flipped), n_alone=len(F['columns']), n_alone_shown=len(shown),
             best_auc=float(1 - best['auroc'] if flipped else best['auroc']), best_lo=float(1 - best['auroc_hi'] if flipped else best['auroc_lo']),
             best_hi=float(1 - best['auroc_lo'] if flipped else best['auroc_hi']))


def readme_feature_groups(src):
    """The five groups of columns that README.md states for the DR+PU model, with their counts."""
    lines = src.lines['README.md']
    start = [i for i, line in enumerate(lines) if re.search(r'\*\*Features\*\* \((\d+) total', line)]
    assert len(start) == 1
    total = int(re.search(r'\((\d+) total', lines[start[0]])[1])
    groups = []
    for line in lines[start[0] + 1:]:
        m = re.match(r'  - (.+?)(?::| \()', line)
        c = re.search(r'\((\d+)\)\s*$', line)
        if not (m and c):
            break
        groups.append((m[1], int(c[1])))
    assert len(groups) == 5 and sum(c for _, c in groups) == total, 'the groups of README.md do not add up to its total'
    return start[0] + 1, total, groups


def agrees(stated, computed, decimals):
    """A stated number agrees when it is within one and a half units of its last stated digit."""
    return abs(stated - computed) <= 1.5 * 10 ** (-decimals)


def decimals_of(text):
    return len(text.split('.')[1]) if '.' in text else 0


def compute_claims(src, F, K):
    """The numbers that the written docs state about the DR+PU pipeline, set against what the saved files and the code give."""
    dr, k = src.dr, src.facts
    res = F['results'].set_index(['reading', 'stage'])
    st3, st1 = res.loc[('dr_step3', 'all')], res.loc[('dr_step1', 'all')]
    d, u = dr[dr['is_drugged'] == 1], dr[dr['is_drugged'] == 0]
    mean_std = {'drugged': (float(d['Y_tilde'].mean()), float(d['Y_tilde'].std())), 'undrugged': (float(u['Y_tilde'].mean()), float(u['Y_tilde'].std()))}
    lab = F['labels'].set_index('universe').loc['DR+PU file']
    rows = []

    def add(item, file, line, stated, gives, verdict):
        rows.append(dict(item=item, file=file, line=line, stated=stated, files_give=gives, verdict='agrees' if verdict else 'differs'))

    for file, pattern in (('README.md', r'OOF AUROC (\d\.\d+), AUPRC (\d\.\d+)'), ('Methodology.md', r'AUROC = (\d\.\d+), AUPRC = (\d\.\d+)\.'),
                          ('CONTRIBUTING.md', r'OOF AUROC (\d\.\d+), AUPRC (\d\.\d+)'), (APP_PAGE, r'OOF AUROC (\d\.\d+), AUPRC (\d\.\d+)')):
        ln, m = src.grab(file, pattern)
        ok = agrees(float(m[1]), st3['auroc'], decimals_of(m[1])) and agrees(float(m[2]), st3['auprc'], decimals_of(m[2]))
        add('Step 3 out-of-fold AUROC and AUPRC', file, ln, f'{m[1]} and {m[2]}', f"{f3(st3['auroc'])} and {f3(st3['auprc'])}", ok)
    pi = k['pi_override']
    for file, pattern, scale in (('Methodology.md', r'π_p ≈ (0\.\d+)', 1.0), ('README.md', r'π_p ≈ (\d+)%', 0.01), ('CONTRIBUTING.md', r'π_p ≈ (\d+)%', 0.01), (APP_PAGE, r'metric: "π_p = (\d+)%"', 0.01)):
        ln, m = src.grab(file, pattern)
        stated = float(m[1]) * scale
        add('Prior of step 2, the share of genes with a liability', file, ln, m[1] if scale == 1.0 else f'{m[1]}%',
            f'{pi:.2f} fixed in the code (line {k["pi_line"]}), not stored in the saved files', agrees(stated, pi, 2))
    for file, pattern, key in (('Methodology.md', r'Drugged genes .{1,3}: mean ≈ (\d\.\d+), std ≈ (\d\.\d+)', 'drugged'), ('Methodology.md', r'Undrugged genes .{1,3}: mean ≈ (\d\.\d+), std ≈ (\d\.\d+)', 'undrugged')):
        ln, m = src.grab(file, pattern)
        mean, std = mean_std[key]
        ok = agrees(float(m[1]), mean, decimals_of(m[1])) and agrees(float(m[2]), std, decimals_of(m[2]))
        add(f'Mean and standard deviation of the pseudo-outcome, {key} genes', file, ln, f'{m[1]} and {m[2]}', f'{mean:.2f} and {std:.2f}', ok)
    ln, m = src.grab(APP_PAGE, r'Drugged gene .{1,3}: mean (\d\.\d+), std (\d\.\d+); Undrugged: mean (\d\.\d+), std (\d\.\d+)')
    for key, i in (('drugged', 1), ('undrugged', 3)):
        mean, std = mean_std[key]
        ok = agrees(float(m[i]), mean, decimals_of(m[i])) and agrees(float(m[i + 1]), std, decimals_of(m[i + 1]))
        add(f'Mean and standard deviation of the pseudo-outcome, {key} genes', APP_PAGE, ln, f'{m[i]} and {m[i + 1]}', f'{mean:.2f} and {std:.2f}', ok)
    for label, column in (('positive', 'positive'), ('drugged_safe', 'drugged_safe'), ('unlabeled', 'unlabeled')):
        ln, m = src.grab('Methodology.md', rf'^\| `{label}` \|.*~([\d,]+)')
        stated = int(m[1].replace(',', ''))
        add(f'Genes labeled {label} among the DB genes', 'Methodology.md', ln, m[1], num(lab[column]), abs(stated - lab[column]) <= 0.005 * lab[column])
    ln, m = src.grab('Methodology.md', r'Propensity AUROC ~(0\.\d+)')
    add('Propensity AUROC', 'Methodology.md', ln, m[1], f3(st1['auroc']), agrees(float(m[1]), st1['auroc'], decimals_of(m[1])))
    F['claims'] = pd.DataFrame(rows)
    K['n_claims'] = len(rows)
    K['n_differ'] = int((F['claims']['verdict'] == 'differs').sum())
    K['pi_doc'] = src.grab('Methodology.md', r'π_p ≈ (0\.\d+)')[1][1]
    K['pi_doc_line'] = src.grab('Methodology.md', r'π_p ≈ (0\.\d+)')[0]
    K['pi'] = pi
    # the app page calls the fixed prior a monotone shift that leaves AUROC and AUPRC unchanged, and says that it raised the AUROC
    ln, m = src.grab(APP_PAGE, r'monotone logit-shift of all scores.*unchanged by theory.*improved OOF AUROC from (\d\.\d+) to (\d\.\d+)')
    assert agrees(float(m[2]), st3['auroc'], decimals_of(m[2])) and float(m[1]) < float(m[2]), 'the AUROC that the app page gives for the fixed prior is not the saved one'
    K.update(app_shift_line=ln, shift_from=m[1], shift_to=m[2])
    app = F['claims'][F['claims']['file'] == APP_PAGE]
    assert (app['verdict'] == 'agrees').all(), 'the app page no longer agrees with the saved files'
    assert (F['claims'][F['claims']['file'].isin(['README.md', 'Methodology.md', 'CONTRIBUTING.md']) & F['claims']['item'].str.startswith('Step 3')]['verdict'] == 'differs').all()


def compute_llm(src, F, K):
    ot, dr = src.ot, src.dr
    releases = set()
    for f in ('README.md', 'Methodology.md'):
        found = re.findall(r'Open Targets Platform (?:release )?(\d{2}\.\d{2})', ' '.join(src.lines[f]))
        assert found, f'{f} names no Open Targets release'
        releases |= set(found)
    assert len(releases) == 1, releases
    release = releases.pop()
    year, month = 2000 + int(release[:2]), int(release[3:])
    end = datetime.date(year, month, calendar.monthrange(year, month)[1])
    assert end < datetime.date.fromisoformat(RELEASE), 'the Open Targets release ends after the release date of the language model'
    dated = [c for c in ot.columns if re.search(r'date|year|time|day', c, re.I)]
    assert not dated, f'the parquet holds a column that may be a date: {dated}'
    ln, _ = src.grab('Methodology.md', r'`0\.25/0\.50/0\.75/1\.0` encoding for Phases I.IV')
    k = dr.set_index('gene_symbol')
    assert set(TEXTBOOK) <= set(k.index)
    listed = k.loc[list(TEXTBOOK)]
    K.update(ot_release=release, ot_month=f'{calendar.month_name[month]} {year}', n_cols=ot.shape[1], n_date_cols=len(dated), phase_line=ln,
             textbook_n=len(TEXTBOOK), textbook_flagged=int(flag_of(listed).sum()), ot_flag_total=int(flag_of(ot).sum()))
    F['facts_llm'] = pd.DataFrame([dict(item=g, label=listed.loc[g, 'safety_label']) for g in TEXTBOOK])


def compute_dr_features(src, F, K):
    """The five groups of columns that README.md states for the DR+PU model: they cannot be recounted, the matrix is built from a database that is not in the repo."""
    line, total, groups = readme_feature_groups(src)
    word = lambda n: n if n.split()[0] in ('Open', 'LOEUF', 'pLoF') else n[0].lower() + n[1:]                    # noqa: E731
    K.update(readme_feat_line=line, dr_feat_total=total, dr_feat_groups=joined(f'{word(name)} {num(c)}' for name, c in groups))


def compute(src):
    verify(src)
    F, K = {}, {}
    compute_counts(src, F, K)
    compute_results(src, F, K)
    compute_features(src, F, K)
    compute_dr_features(src, F, K)
    compute_claims(src, F, K)
    compute_llm(src, F, K)
    K['present'] = src.present
    K['splitters'] = sorted(src.all_splitters)
    K['learners'] = sorted(src.all_learners)
    return F, K


# ---------------------------------------------------------------------------------------------------------------------------------------------
# the page: tables
# ---------------------------------------------------------------------------------------------------------------------------------------------
def result_cells(r):
    if not r['shown']:
        return ['counts only', 'counts only']
    return [interval(r['auroc'], r['auroc_lo'], r['auroc_hi']), interval(r['auprc'], r['auprc_lo'], r['auprc_hi'])]


def experiments_table(F, K, k):
    e1 = k['e1_params']
    rows = [
        ['Drugged-only model (`train_drugged_only`, Model B)', 'Which genes that entered a trial carry a safety flag?', 'flag against no flag', num(K['e1_n']),
         f"{K['n_saved']} saved columns: {joined(f'{int(r.columns)} {GROUP_WORD[r.group]}' for r in F['groups'].itertuples())}", f"{K['folds']}-fold out of fold, seed {K['seed_code']}",
         '`predictions_oof.csv`, `shap_summary.csv`, `shap_values.npy`'],
        ['All-labeled model (`train_and_evaluate`, Model A)', 'Which targets carry a flag, learned from flagged genes and drugged genes without a flag?', 'flag, or drugged without a flag; the rest is left out',
         f"{num(K['ot_flag'] + K['ot_safe'])} trained, {num(K['ot_n'])} scored", 'built by the script, with the same four groups', f"{K['folds']}-fold, printed only", 'none: its four output files are not in the repo'],
        ['DR+PU step 1, propensity', 'Which genes enter a trial?', 'entered a trial against not', num(K['dr_n']), f"{K['dr_feat_total']} columns, as `README.md` states", f"{K['folds']}-fold cross-fitted", '`pi_hat` of `dr/predictions.csv`'],
        ['DR+PU step 2, prior', 'What share of genes carry a liability?', 'none', 'none', 'none', f"a constant, {K['pi']:.2f}, not an estimate", 'none'],
        ['DR+PU step 3, outcome model', 'Which drugged genes carry a flag, with weights for selection and for unlabeled negatives?', 'flag against no flag', num(K['dr_drugged']),
         f"the same {K['dr_feat_total']}", f"{K['folds']}-fold cross-fitted", '`m_hat` of `dr/predictions.csv`'],
        ['DR+PU step 5, deployed score', 'A score for any gene, fitted to the pseudo-outcome of step 4', 'the pseudo-outcome, which holds the label of a drugged gene', num(K['dr_n']),
         f"the same {K['dr_feat_total']}", 'none out of fold', '`safety_score_dr` of `dr/predictions.csv`'],
        ['Single Open Targets columns (made by this page, no model)', 'How far does one column go on its own?', 'flag against no flag', num(K['e1_n']), f"{K['n_ot']} columns, one at a time", 'no fitting, so no folds', '`docs/grid/numbers.csv`'],
    ]
    return md(['Reading', 'Question', 'Label', 'Genes', 'Columns', 'Evaluation', 'Saved in the repo'], rows, right=(3,))


def words_table():
    rows = [
        ['Gene', 'One Ensembl id of the Open Targets file, which is the row of every table of this page.'],
        ['Flag', '`hasSafetyEvent` is `-1`: Open Targets lists a known safety event for the target. No other value occurs, so a gene without a flag is not confirmed free of one.'],
        ['Drugged', 'A gene with a value in `maxClinicalTrialPhase`, the maximum phase of its record. The code reads a gene with no value as not drugged.'],
        ['positive', 'A flagged gene, whether or not it entered a trial (the label `positive` of the saved predictions).'],
        ['drugged_safe', 'A drugged gene without a flag. It is an assumed negative, not a confirmed one.'],
        ['unlabeled', 'A gene with no phase and no flag. Its label is unknown, and it is not safe by default.'],
        ['Events', 'The genes whose label is the flag. The base rate is the events over the genes, which is also the AUPRC of a random ranking.'],
        ['Out of fold', 'A gene scored by a model that did not train on it.'],
    ]
    return md(['Word', 'What it counts'], rows)


def labels_table(F):
    rows = [[r.universe, num(r.genes), num(r.positive), num(r.positive_drugged), num(r.positive_no_phase), num(r.drugged_safe), num(r.unlabeled)] for r in F['labels'].itertuples()]
    return md(['Universe', 'Genes', 'positive', 'of which drugged', 'of which no phase', 'drugged_safe', 'unlabeled'], rows, right=(1, 2, 3, 4, 5, 6))


def setup_table(F):
    names = {'drugged_only': 'Drugged-only model', 'all_labeled': 'All-labeled model (no saved result)', 'dr_step1': 'DR+PU step 1, propensity', 'dr_step3': 'DR+PU step 3, outcome model',
             'dr_step5': 'DR+PU step 5, deployed score'}
    rows = []
    for r in F['setup'].itertuples():
        if r.reading == 'dr_step5':
            rows.append([names[r.reading], num(r.universe), num(r.trained), num(r.scored), 'no folds', 'no folds', 'no events', 'none'])
            continue
        rows.append([names[r.reading], num(r.universe), num(r.trained), num(r.scored), span(r.fold_train_min, r.fold_train_max), span(r.fold_score_min, r.fold_score_max), num(r.events), f3(r.base_rate)])
    return md(['Reading', 'Genes in the universe', 'Trained on', 'Scored', 'Per fold, trained', 'Per fold, scored', 'Events', 'Base rate'], rows, right=(1, 2, 3, 4, 5, 6, 7))


def dedup_table(F):
    def cell(x):
        return 'not in the file' if pd.isna(x) else num(x)
    rows = [[r.file, cell(r.rows), cell(r.unique_ids), cell(r.duplicate_ids), cell(r.rows_with_symbol), cell(r.unique_symbols), cell(r.duplicate_symbols)] for r in F['dedup'].itertuples()]
    return md(['File', 'Rows', 'Unique Ensembl ids', 'Duplicate ids', 'Rows with a symbol', 'Unique symbols', 'Duplicate symbols'], rows, right=(1, 2, 3, 4, 5, 6))


def cut_words_table(K):
    rows = [
        ['each gene once', 'one row per Ensembl id', 'holds in every file (section 2)'],
        ['each drug once', 'a drug counts once', 'not applicable: the files name no drug'],
        ['cancer only, not cancer', 'by the disease of the trial', 'not run: the files hold no indication (section 3)'],
        ['single agent, combination, mixed', 'by the regimen', 'not applicable: a row is a gene, and a gene has no regimen'],
        ['stage', 'by the maximum phase', f"counts for every universe and results for the {num(K['e1_n'])} drugged genes (sections 3 and 4)"],
        ['driver gene', 'by `isCancerDriverGene`', 'counts only, and it is no cancer cut (section 3)'],
    ]
    return md(['Cut', 'What it keeps', 'In these files'], rows)


def stage_table(F, universe):
    sub = F['stage'][F['stage']['universe'] == universe]
    rows = []
    for r in sub.itertuples():
        code = float(r.phase) if r.phase not in ('none', 'all') else None
        rows.append([r.phase, ROMAN[code] if code is not None else ('not drugged' if r.phase == 'none' else 'all'), num(r.genes), num(r.positive), num(r.drugged_safe), num(r.unlabeled), f3(r.flag_rate)])
    return md(['Max phase code', 'Phase', 'Genes', 'positive', 'drugged_safe', 'unlabeled', 'Flag rate'], rows, right=(2, 3, 4, 5, 6))


def driver_table(F):
    rows = [[r.group, num(r.genes), num(r.flagged), num(r.not_flagged), f3(r.flag_rate)] for r in F['driver'].itertuples()]
    return md(['Cancer driver flag', 'Drugged genes', 'positive', 'drugged_safe', 'Flag rate'], rows, right=(1, 2, 3, 4))


def results_table(F, K):
    r = F['results'].set_index(['reading', 'stage'])
    rows = []

    def line(key):
        x = r.loc[(key, 'all')]
        return [READING_LABEL[key], num(x['genes']), num(x['events']), f3(x['base_rate'])] + result_cells(x)
    rows += [line('drugged_only'), line('dr_step3')]
    p = F['paired'].iloc[0]
    rows.append([f"Difference, step 3 minus drugged-only, on the {num(p['genes'])} genes that both score (paired)", num(p['genes']), num(p['events']), f3(p['base_rate']),
                 signed(p['auroc_diff'], p['auroc_diff_lo'], p['auroc_diff_hi']), signed(p['auprc_diff'], p['auprc_diff_lo'], p['auprc_diff_hi'])])
    rows += [line('dr_step1'), line('dr_und_outcome'), line('dr_und_deployed')]
    x = r.loc[('dr_deployed_drugged', 'all')]
    rows.append([READING_LABEL['dr_deployed_drugged'] + ' (not a result)', num(x['genes']), num(x['events']), f3(x['base_rate']), f"{f3(x['auroc'])} in sample, no interval", f"{f3(x['auprc'])} in sample, no interval"])
    return md(['Reading', 'Genes', 'Events', 'Base rate', 'AUROC [95% interval]', 'AUPRC [95% interval]'], rows, right=(1, 2, 3))


def stage_results_table(F):
    rows = []
    for key in ('drugged_only', 'dr_step3'):
        sub = F['results'][(F['results']['reading'] == key) & (F['results']['stage'] != 'all')]
        for r in sub.itertuples():
            rows.append([READING_LABEL[key], r.stage, ROMAN[float(r.stage)], num(r.genes), num(r.events), f3(r.base_rate)] + result_cells(r._asdict()))
    return md(['Reading', 'Max phase code', 'Phase', 'Genes', 'Events', 'Base rate', 'AUROC [95% interval]', 'AUPRC [95% interval]'], rows, right=(3, 4, 5))


def groups_table(F):
    rows = [[r.group, num(r.columns), r.detail, f4(r.mean_abs_shap), pct(r.share), f'`{r.top_feature}` {f4(r.top_value)}'] for r in F['groups'].itertuples()]
    return md(['Group', 'Columns', 'Of which', 'Sum of mean absolute SHAP', 'Share', 'Largest column'], rows, right=(1, 3, 4))


def shap_top_table(F, n=15):
    rows = [[r['rank'], f"`{r['feature']}`", r['group'], f4(r['mean_abs_shap']), pct(r['share'])] for _, r in F['shap'].head(n).iterrows()]
    return md(['Rank', 'Column', 'Group', 'Mean absolute SHAP', 'Share of the total'], rows, right=(0, 3, 4))


def columns_table(F):
    rows = []
    for r in F['columns'].itertuples():
        genes = num(r.genes) if pd.isna(r.present) else f'{num(r.genes)} ({num(r.present)} present)'
        if r.shown:
            reading = ('present goes with a flag' if r.auroc_lo > 0.5 else 'present goes with no flag' if r.auroc_hi < 0.5 else 'no ranking shown') if r.basis.startswith('present') else (
                'higher value goes with a flag' if r.auroc_lo > 0.5 else 'lower value goes with a flag' if r.auroc_hi < 0.5 else 'no ranking shown')
            cell = interval(r.auroc, r.auroc_lo, r.auroc_hi)
        else:
            reading, cell = 'counts only', 'counts only'
        rows.append([f'`{r.column}`', r.basis, genes, num(r.flagged), cell, reading, r.shap_rank, f4(r.mean_abs_shap)])
    return md(['Column', 'Basis', 'Genes', 'Flagged', 'AUROC alone [95% interval]', 'Reading', 'SHAP rank', 'Mean absolute SHAP'], rows, right=(2, 3, 6, 7))


def llm_table(K):
    rows = [
        ['The labels have no date', f"`hasSafetyEvent` is a flag of release {K['ot_release']}, and none of the {K['n_cols']} columns of the file is a date", f'flags dated after {RELEASE}'],
        ['The snapshot ends before the release', f"release {K['ot_release']} is named for {K['ot_month']}, which ends before {RELEASE}", 'a later Open Targets release, which is not in this repo'],
        ['The textbook liabilities are public', f"{K['textbook_flagged']} of {K['textbook_n']} well known safety targets carry the flag ({', '.join(f'`{g}`' for g in TEXTBOOK)})", 'an outcome that a language model has not read'],
        ['No forward list exists', 'no forecast of this repo was frozen before a label was read', 'forecasts frozen before the later release is read, graded once on the genes it newly flags'],
        ['Nothing paid was run', 'no model was called to make this page', 'a recorded spend'],
    ]
    return md(['Reason', 'What the files give', 'What an arm would need'], rows)


def claims_table(F):
    rows = [[r.item, f'`{r.file}` line {r.line}', r.stated, r.files_give, r.verdict] for r in F['claims'].itertuples()]
    return md(['Statement', 'Where', 'Stated', 'Saved files and code give', 'Verdict'], rows)


def not_run_table(K, k):
    absent = lambda paths: [p for p in paths if not K['present'][p]]                   # noqa: E731
    a_missing, in_missing = absent(MODEL_A_FILES), absent(INPUT_FILES)
    rows = [
        ['Result of the all-labeled model (Model A)', f'its output files are not in the repo ({len(a_missing)} of {len(MODEL_A_FILES)} absent), and no run of it is saved', 'the saved output, or a run of the script'],
        ['Retraining any model here', f'{len(in_missing)} of the {len(INPUT_FILES)} input files that the scripts read are absent, among them `backend/humanproof.db`', 'those files, which `README.md` says how to fetch and build'],
        ['Feature group ablations (drop Open Targets, genetics or expression)', 'no script holds one, and the saved model has genetics and expression values for few genes', 'the feature matrix, built from `backend/humanproof.db`, and XGBoost'],
        ['A split in time or a held-back test', 'the labels carry no dates, and every number is out of fold', 'dated labels, such as two Open Targets releases, and a gene list fixed before the labels are read'],
        ['Cancer only and not cancer', 'the files hold no indication of any trial', 'a field that names the disease of the trial, or an area for each target'],
        ['Single agent, combination, mixed', 'a row is a gene, and the files name no drug', 'a table of the drugs and their regimens for each target'],
        ['Counts of distinct drugs and trials', 'no column of the files names a drug or a trial', 'a table of the drugs and the trials behind each flag and each phase'],
        ['Folds that group paralogs', f"no family id is in the files, and the scripts split with `{joined(K['splitters'])}` only", 'a family id for each gene'],
        ['Out-of-fold evaluation of the deployed score on drugged genes', 'the regressor is fitted to the pseudo-outcome, which holds the label of these genes', 'a cross-fitted step 5'],
        ['The prior of step 2 as an estimate', f"the code fixes it at {K['pi']:.2f}", 'a run with `PI_OVERRIDE` set to `None`'],
        ['A table of the DR+PU feature importance', f"it is saved as pictures only, and {len(absent(DR_FILES))} of {len(DR_FILES)} model and SHAP files are absent", '`gene_shap_dr.json`, or a run of the script'],
        [f'The {LLM} only arm', 'see section 6', f'flags dated after {RELEASE} and forecasts frozen before they are read'],
    ]
    return md(['Cell', 'Why it is not run', 'What it would need'], rows)


# ---------------------------------------------------------------------------------------------------------------------------------------------
# the page: figures and sentences
# ---------------------------------------------------------------------------------------------------------------------------------------------
def values(F, K, src):
    """Every figure of the text, each asserted against the tables before its sentence is made."""
    k = src.facts
    res = F['results'].set_index(['reading', 'stage'])
    v = {}
    e1, st3, st1 = res.loc[('drugged_only', 'all')], res.loc[('dr_step3', 'all')], res.loc[('dr_step1', 'all')]
    uo, ud = res.loc[('dr_und_outcome', 'all')], res.loc[('dr_und_deployed', 'all')]
    pair = F['paired'].iloc[0]
    # claims of the text that the tables must hold
    assert pair['auroc_diff_lo'] < 0 < pair['auroc_diff_hi'], 'the paired AUROC difference no longer holds 0'
    assert K['min_flag_score'] > K['max_safe_score'], 'a flagged drugged gene scores below an unflagged one in sample'
    assert float(res.loc[('dr_deployed_drugged', 'all'), 'auroc']) == 1.0
    assert K['pi'] != float(K['pi_doc']), 'the code and Methodology.md now state the same prior'
    p1, p2 = k['e1_params'], k['dr_clf']
    same = ('n_estimators', 'max_depth', 'learning_rate', 'subsample', 'colsample_bytree', 'eval_metric', 'early_stopping_rounds', 'random_state', 'tree_method')
    assert all(p1[x] == p2[x] for x in same), 'the DR+PU classifiers no longer share the settings of the drugged-only model'
    assert k['e2_params']['max_depth'] != p1['max_depth'] and all(k['e2_params'][x] == p1[x] for x in same if x != 'max_depth')
    assert k['e1_params']['early_stopping_rounds'] and k['dr_clf']['early_stopping_rounds']
    assert abs(K['ot_share'] + F['groups'].set_index('group').loc[['LOEUF', 'pLoF', 'Expression'], 'share'].sum() - 1) < 1e-9
    assert K['ot_share'] > 0.9 and K['shap_rho'] > 0.5 and abs(K['shap_rho_shuffled']) < 0.1, 'the SHAP matrix is not aligned with the drugged genes'
    assert K['textbook_flagged'] >= 0 and K['n_date_cols'] == 0
    assert K['dr_symbols_in_list'] == K['dr_n'] and K['gene_list_n'] == int(src.grab('Methodology.md', r'~([\d,]+) unique protein-coding genes')[1][1].replace(',', ''))
    assert all(not K['present'][p] and p in src.ignored for p in ABSENT), 'a file that the page calls absent is tracked, or no rule of .gitignore excludes it'
    assert all(K['present'][p] for p in DR_PICTURES) and set(GITIGNORE) <= set(src.gitignore), 'a picture of the DR+PU importance is gone, or a rule of .gitignore that the page names'
    assert K['n_expr_old'] > 0 and K['drivers_all'] > K['drivers_drugged']
    # the figures
    v.update(draws=num(DRAWS), seed=str(SEED), release=RELEASE, llm=LLM, ot_release=K['ot_release'], ot_month=K['ot_month'], min_genes=num(MIN_GENES), min_class=num(MIN_CLASS), folds=str(K['folds']))
    v.update(e1_stop_line=str(k['e1_stop_line']), dr_stop_line=str(k['dr_stop_line']), pi_line=str(k['pi_line']), genes_line=str(k['genes_line']), symbol_line=str(k['symbol_line']), phase_line=str(K['phase_line']))
    v.update(e1_trees=num(p1['n_estimators']), e1_depth=str(p1['max_depth']), e1_lr=str(p1['learning_rate']), e1_sub=str(p1['subsample']), e1_col=str(p1['colsample_bytree']), e1_stop=str(p1['early_stopping_rounds']),
             e2_depth=str(k['e2_params']['max_depth']), dr_reg_trees=num(k['dr_reg']['n_estimators']))
    v.update(e1_n=num(K['e1_n']), e1_events=num(K['e1_events']), e1_base=f3(e1['base_rate']), e1_auc=interval(e1['auroc'], e1['auroc_lo'], e1['auroc_hi']), e1_ap=interval(e1['auprc'], e1['auprc_lo'], e1['auprc_hi']),
             dr3_n=num(st3['genes']), dr3_auc=interval(st3['auroc'], st3['auroc_lo'], st3['auroc_hi']), dr3_ap=interval(st3['auprc'], st3['auprc_lo'], st3['auprc_hi']),
             dr3_auc_point=f3(st3['auroc']), dr3_ap_point=f3(st3['auprc']), e1_auc_point=f3(e1['auroc']),
             shared_n=num(K['shared_n']), paired_auc=signed(pair['auroc_diff'], pair['auroc_diff_lo'], pair['auroc_diff_hi']),
             und_n=num(K['und_n']), und_pos=num(K['und_pos']), und_neg=num(K['und_neg']), und_base=f3(uo['base_rate']),
             und_out_auc=interval(uo['auroc'], uo['auroc_lo'], uo['auroc_hi']), und_dep_auc=interval(ud['auroc'], ud['auroc_lo'], ud['auroc_hi']),
             und_zero=num(K['und_zero']), und_flag_zero=num(K['und_flag_zero']), min_flag_score=f3(K['min_flag_score']), max_safe_score=f3(K['max_safe_score']),
             insample_auc=f3(res.loc[('dr_deployed_drugged', 'all'), 'auroc']),
             ot_n=num(K['ot_n']), ot_flag=num(K['ot_flag']), ot_flag_drugged=num(K['ot_flag_drugged']), ot_safe=num(K['ot_safe']), ot_unl=num(K['ot_unl']), dr_n=num(K['dr_n']), dr_drugged=num(K['dr_drugged']), dr_events=num(K['dr_flag_drugged']),
             dr_flag=num(K['dr_flag']), dr_und_flag=num(K['dr_flag_no_phase']), n_cols=num(K['n_cols']), n_features=num(K['n_cols'] - 3), gene_list_n=num(K['gene_list_n']),
             meth_universe=src.grab('Methodology.md', r'(~[\d,]+) unique protein-coding genes')[1][1],
             drivers_all=num(K['drivers_all']), drivers_undrugged=num(K['drivers_undrugged']), drivers_drugged=num(K['drivers_drugged']), n_sym=num(K['e1_sym']), n_nosym=num(K['e1_nosym']),
             n_sym_drugged=num(K['n_sym_drugged']), genes_list_n108=num(K['genes_list_n']), n_saved=num(K['n_saved']), n_ot=num(K['n_ot']), n_nonot=num(K['n_nonot']), n_expr=num(K['n_expr']), n_expr_old=num(K['n_expr_old']),
             shap_rows=num(K['shap_rows']), shap_cols=num(K['shap_cols']), ot_share=pct(K['ot_share']), nonot_share=pct(1 - K['ot_share']), shap_insample=f4(K['shap_insample_auc']),
             rho=f3(K['shap_rho']), rho_shuffled=f3(K['shap_rho_shuffled']), best_col=K['best_col'], best_auc=interval(K['best_auc'], K['best_lo'], K['best_hi']), n_alone=num(K['n_alone']),
             textbook_n=num(K['textbook_n']), textbook_flagged=num(K['textbook_flagged']), n_claims=num(K['n_claims']), n_differ=num(K['n_differ']), pi=f"{K['pi']:.2f}", pi_doc=K['pi_doc'])
    v.update(readme_feat_line=str(K['readme_feat_line']), dr_feat_total=num(K['dr_feat_total']), dr_feat_groups=K['dr_feat_groups'], dr_picture=DR_PICTURES[0], pi_doc_line=str(K['pi_doc_line']))
    v.update(clip_range=f"[{k['clip_range'][0]:g}, {k['clip_range'][1]:g}]", clip_lo=f"{k['clip_range'][0]:g}", clip_line=str(k['clip_line']))
    v.update(shift_comment_lines=k['shift_comment_lines'], weight_line=str(k['weight_line']), app_shift_line=str(K['app_shift_line']), shift_from=K['shift_from'], shift_to=K['shift_to'])
    codes = [str(c) for c in PHASES]
    v['phase_codes'] = joined(codes)
    # the by-stage sentence and the one small cell
    e1s = {r.stage: r for r in F['results'].itertuples() if r.reading == 'drugged_only' and r.stage != 'all'}
    st3s = {r.stage: r for r in F['results'].itertuples() if r.reading == 'dr_step3' and r.stage != 'all'}
    shown = [c for c in codes if e1s[c].shown]
    small = [c for c in codes if not e1s[c].shown]
    assert len(small) == 1 and len(shown) == len(codes) - 1, 'the stage cells that show counts only are no longer a single one'
    assert not st3s[small[0]].shown and all(st3s[c].shown for c in shown), 'the stage cells of step 3 differ from those of the drugged-only model'
    c = small[0]
    v['stage_text'] = joined(f'{f3(e1s[x].auroc)} at code {x}' for x in shown)
    v['small_text'] = (f'code {c} holds {num(e1s[c].genes)} genes with {num(e1s[c].events)} flagged in the drugged-only model '
                       f'and {num(st3s[c].genes)} genes with {num(st3s[c].events)} flagged in step 3')
    return v


def sentences(v, F, K, src):
    """The twelve statements of the section "At a glance", each made from the figures and each one line of one sentence."""
    k = src.facts
    s = {}
    stage_n = F['stage'][(F['stage']['universe'] == 'Open Targets file') & F['stage']['phase'].isin([str(c) for c in PHASES])].set_index('phase')['genes']
    s['s_summary'] = (f"The safety model asks whether Open Targets lists a known safety event for a gene, learning from its Open Targets columns, loss of function genetics, constraint and cell type expression, "
                      f"with a drugged-only classifier and a DR+PU pipeline that scores all {v['dr_n']} genes.")
    s['s_unit'] = (f"One row is one gene, an Ensembl id of the Open Targets file (never a drug, a trial or a patient), and every file read here holds one row per id with no duplicate id and no duplicate symbol.")
    s['s_label'] = (f"A positive is a gene whose `hasSafetyEvent` is `-1` ({v['ot_flag']} of the {v['ot_n']} genes of the file, {v['ot_flag_drugged']} of them drugged), an assumed negative is a drugged gene without that flag ({v['ot_safe']}), "
                    f"and a gene with neither a phase nor a flag is unlabeled ({v['ot_unl']}) and carries no label, so no gene is confirmed free of a safety event.")
    s['s_population'] = (f"The drugged-only model reads the {v['e1_n']} genes that have a maximum phase ({v['e1_events']} flagged), and the DR+PU pipeline reads {v['dr_n']} genes "
                         f"({v['dr_drugged']} drugged, {v['dr_events']} of them flagged) and fits its outcome model on the {v['dr_drugged']} drugged ones.")
    s['s_test'] = (f"Stratified {v['folds']}-fold cross-validation with shuffling and seed {v['seed']} scores each gene once out of fold (XGBoost, early stopping on the scored fold), "
                   f"with 95% percentile intervals from {v['draws']} bootstrap draws over genes at seed {v['seed']}, no split in time and no test set held back.")
    s['s_metric'] = (f"AUROC and average precision (AUPRC), each with a 95% bootstrap interval over genes, beside the number of genes, the number of flagged genes and the base rate, which is the AUPRC of a random ranking.")
    s['s_baseline'] = (f"The repo compares each model only with chance (AUROC 0.5, AUPRC equal to the base rate) and fits no baseline model, and for context this page reads each of the {v['n_alone']} Open Targets columns alone, "
                       f"the best being `{v['best_col']}` at AUROC {v['best_auc']} in its best direction, against {v['e1_auc_point']} for the drugged-only model.")
    s['s_not_comparable'] = (f"These are out-of-fold AUROCs against assumed negatives at one snapshot with random folds, so they cannot be compared with the forecast of a later event, with a score on confirmed outcomes, "
                             f"with a deployed score that was fitted to its own labels, or with a test that groups paralogs or splits in time.")
    s['s_how_many'] = (f"The drugged-only model scores {v['e1_n']} genes ({v['e1_events']} flagged, base rate {v['e1_base']}) and the DR+PU pipeline {v['dr_n']} genes ({v['dr_drugged']} drugged, {v['dr_events']} of them flagged), "
                       f"and the files name no drug and no trial, so drug and trial counts are not in this repo's files.")
    s['s_cuts'] = (f"Every file holds one row per gene, cancer only is not run because no file holds an indication ({v['drivers_drugged']} of the {v['e1_n']} drugged genes carry a cancer driver flag), "
                   f"combination does not apply to genes, and the stage cut by maximum phase holds {num(stage_n['0.25'])}, {num(stage_n['0.5'])}, {num(stage_n['0.75'])} and {num(stage_n['1.0'])} drugged genes at the codes {v['phase_codes']}.")
    top = F['shap'].head(3)
    lead = joined(f"{r.feature} {f4(r.mean_abs_shap)}" for r in top.itertuples())
    s['s_features'] = (f"In the drugged-only model the {v['n_ot']} Open Targets columns carry {v['ot_share']} of the mean absolute SHAP ({lead} lead), probably because the other {v['n_nonot']} columns are filled for few genes, "
                       f"and the DR+PU importance is saved only as a picture.")
    s['s_llm'] = (f"Not run: the labels are Open Targets {v['ot_release']} safety flags, which carry no date and were all published before the {v['release']} release of {v['llm']}, "
                  f"and well known safety targets are public ({v['textbook_flagged']} of {v['textbook_n']} listed ones carry the flag), so a language model would be graded on recall and not on a forecast.")
    for name, text in s.items():
        assert '\n' not in text and len(text) < LIMIT, (name, len(text))
    return s


def findings(v):
    """The statements of the section "What the page finds"."""
    return {
        'f1': f"The drugged-only model reads AUROC {v['e1_auc']} and AUPRC {v['e1_ap']} out of fold on {v['e1_n']} drugged genes ({v['e1_events']} flagged, base rate {v['e1_base']}).",
        'f2': (f"The outcome model of the DR+PU pipeline reads AUROC {v['dr3_auc']} and AUPRC {v['dr3_ap']} on {v['dr3_n']} drugged genes, and its AUROC minus that of the drugged-only model on the {v['shared_n']} genes "
               f"that both score is {v['paired_auc']}, an interval that holds 0."),
        'f3': f"By stage the drugged-only AUROC reads {v['stage_text']}, and the cell of {v['small_text']}, so it shows counts only.",
        'f4': (f"The deployed DR+PU score cannot be evaluated on drugged genes: every flagged drugged gene scores above every unflagged one in sample (AUROC {v['insample_auc']}), because the pseudo-outcome "
               f"it is fitted to holds the label."),
        'f5': (f"On the {v['und_n']} undrugged genes, whose labels no step read, the deployed score reads AUROC {v['und_dep_auc']} and the outcome model {v['und_out_auc']} for the {v['und_pos']} flagged genes "
               f"(base rate {v['und_base']}), against unlabeled genes that are unknown and not safe."),
        'f6': f"{v['n_differ']} of the {v['n_claims']} numbers that the written docs state about the outcome model, the prior, the pseudo-outcome and the labels differ from the saved files and the code (section 8).",
        'f7': f"Cancer only is not run, combination does not apply to genes, the files name no drug and no trial, and no {v['llm']} arm exists (sections 3 and 6).",
    }


DIRECT = ('ot_release', 'draws', 'seed', 'e1_stop_line', 'dr_stop_line', 'e1_trees', 'e1_depth', 'e1_lr', 'e1_sub', 'e1_col', 'e1_stop', 'e2_depth', 'dr_reg_trees', 'folds', 'dr_und_flag', 'dr_flag',
          'n_cols', 'n_features', 'gene_list_n', 'meth_universe', 'dr_n', 'phase_line', 'phase_codes', 'e1_n', 'drivers_all', 'drivers_undrugged', 'min_genes', 'min_class', 'dr3_n', 'shared_n', 'und_n',
          'und_pos', 'pi_line', 'und_neg', 'und_zero', 'und_flag_zero', 'min_flag_score', 'max_safe_score', 'small_text', 'n_saved', 'shap_rows', 'shap_cols', 'ot_share', 'n_nonot', 'nonot_share',
          'genes_list_n108', 'genes_line', 'symbol_line', 'n_sym', 'n_sym_drugged', 'n_nosym', 'n_expr', 'n_expr_old', 'shap_insample', 'rho', 'rho_shuffled', 'best_col', 'n_alone', 'best_auc', 'e1_auc',
          'dr_picture', 'readme_feat_line', 'dr_feat_total', 'dr_feat_groups', 'llm', 'release', 'n_claims', 'n_differ', 'dr3_auc_point', 'dr3_ap_point', 'pi', 'pi_doc_line', 'pi_doc',
          'clip_range', 'clip_lo', 'clip_line', 'shift_comment_lines', 'weight_line', 'app_shift_line', 'shift_from', 'shift_to')


def parts(F, K, src):
    """Every figure, sentence and table that the template names, and no other."""
    v = values(F, K, src)
    k = src.facts
    out = {name: v[name] for name in DIRECT}
    out.update(sentences(v, F, K, src))
    out.update(findings(v))
    out.update(table_experiments=experiments_table(F, K, k), table_words=words_table(), table_labels=labels_table(F), table_setup=setup_table(F), table_dedup=dedup_table(F),
               table_cut_words=cut_words_table(K), table_stage_ot=stage_table(F, 'Open Targets file'), table_stage_dr=stage_table(F, 'DR+PU file'), table_driver=driver_table(F),
               table_results=results_table(F, K), table_stage_results=stage_results_table(F), table_groups=groups_table(F), table_shap_top=shap_top_table(F), table_columns=columns_table(F),
               table_llm=llm_table(K), table_claims=claims_table(F), table_not_run=not_run_table(K, k))
    return out


# ---------------------------------------------------------------------------------------------------------------------------------------------
# build, write and check
# ---------------------------------------------------------------------------------------------------------------------------------------------
KEYS = {'labels': ['universe'], 'setup': ['reading'], 'dedup': ['file'], 'stage': ['universe', 'phase'], 'driver': ['group'], 'results': ['reading', 'stage'], 'paired': [], 'shap': ['feature'],
        'groups': ['group'], 'columns': ['column'], 'claims': ['file', 'line', 'item'], 'facts_llm': ['item']}


def numbers_csv(F, K):
    """Every table of the page, and the facts behind its prose, in one long CSV: table, row, field, value."""
    def fmt(x):
        if isinstance(x, (bool, np.bool_)):
            return str(bool(x))
        if isinstance(x, (int, np.integer)):
            return str(int(x))
        if isinstance(x, (float, np.floating)):
            return '' if np.isnan(x) else f'{x:.6g}'
        return '' if x is None else str(x)
    rows = []
    for name, df in F.items():
        keys = KEYS[name]
        for i, r in enumerate(df.to_dict('records')):
            row = '|'.join(str(r[c]) for c in keys) if keys else str(i)
            for field in df.columns:
                if field not in keys:
                    rows.append((name, row, field, fmt(r[field])))
    for field, value in sorted(K.items()):
        if isinstance(value, (int, float, str, bool, np.integer, np.floating)):
            rows.append(('facts', '', field, fmt(value)))
    out = io.StringIO()
    w = csv.writer(out, lineterminator='\n')
    w.writerow(('table', 'row', 'field', 'value'))
    w.writerows(rows)
    return out.getvalue()


def render(F, K, src, template=None):
    template = read_text(TEMPLATE) if template is None else template
    return template.format(**parts(F, K, src))


def build(src=None):
    src = Sources() if src is None else src
    F, K = compute(src)
    page = render(F, K, src)
    return dict(F=F, K=K, src=src, page=page, numbers=numbers_csv(F, K))


def main(argv=None):
    ap = argparse.ArgumentParser(description=(__doc__ or '').splitlines()[0])
    ap.add_argument('--check', action='store_true', help='fail when the committed page or CSV is not what the files give')
    ap.add_argument('--doc', default=DOC)
    ap.add_argument('--numbers', default=NUMBERS)
    a = ap.parse_args(argv)
    made = build()
    if a.check:
        stale = []
        for path, text in ((a.doc, made['page']), (a.numbers, made['numbers'])):
            if not os.path.exists(path) or read_text(path) != text:
                stale.append(os.path.relpath(path, ROOT) if path.startswith(ROOT) else path)
        if stale:
            sys.exit(f"{' and '.join(stale)} is stale: run `uv run --with pandas --with numpy --with pyarrow --with scikit-learn python scripts/grid_doc.py`")
        return 0
    for path, text in ((a.doc, made['page']), (a.numbers, made['numbers'])):
        with open(path, 'w', encoding='utf-8', newline='\n') as f:
            f.write(text)
    return 0


if __name__ == '__main__':
    sys.exit(main())
