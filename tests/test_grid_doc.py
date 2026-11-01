"""Tests of the grid page (docs/GRID.md) and of what makes it.

    uv run --with pandas --with numpy --with pyarrow --with scikit-learn --with pytest python -m pytest tests/test_grid_doc.py -q

They check that the committed page and CSV are what the files give, that every figure of the template is made and every one made is used, that the prose types no figure,
that a claim of the text that the files do not hold stops the maker, that the counts add up, that the bootstrap is the AUROC and the average precision of scikit-learn on every draw,
that no data file or column of the repo names a drug or a trial, that a file is in the repo when git tracks it, that the page says what it does not hold,
and that the page keeps its rules: one sentence a line, no dash, twelve short statements at a glance, and no other project named.
Each check has a planted bug, so that a check that sees nothing cannot pass.
"""
import copy
import csv
import io
import os
import re
import shutil
import string
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import average_precision_score, roc_auc_score

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import grid_doc as G  # noqa: E402

SENTENCES = re.compile(r'[.!?] +[A-Z`]')
MARKER = re.compile(r'^\s*\d+\. ')
DASHES = (chr(0x2014), chr(0x2013))                                     # the em dash and the en dash
ALLOWED_NUMBERS = {'0', '95%'}                                          # zero (the null) and the level of the intervals: fixed by the method, not by a file
GLANCE = ['Summary', 'Unit', 'Label', 'Population', 'Test', 'Metric', 'Baseline', 'Not comparable', 'How many', 'Cuts', 'Features that matter', 'Language model only']
OTHER_PROJECTS = re.compile('|'.join(['hyth' + 'erx', 'atl' + 'as', 'vau' + 'lt', 'asset' + 'select', 'asset-' + 'forecast', 'trial-' + 'forecast', 'early-' + 'asset', 'cancer-' + 'portal']), re.I)
PROTOCOL_ID = re.compile(r'\bP\d{3}\b')


def read(path):
    with open(path, encoding='utf-8') as f:
        return f.read()


@pytest.fixture(scope='module')
def made():
    return G.build()


def named(template):
    return {name for _, name, _, _ in string.Formatter().parse(template) if name}


def sentence_problems(page):
    """The lines of prose that hold two sentences: outside the tables, the headings and the code, a line holds at most one."""
    bad, in_code = [], False
    for line in page.splitlines():
        if line.startswith('```'):
            in_code = not in_code
        elif not in_code and not line.startswith(('|', '#')) and SENTENCES.search(MARKER.sub('', line)):
            bad.append(line[:90])
    return bad


def dash_problems(text):
    return [d for d in DASHES if d in text]


def typed_numbers(template):
    """The numbers that the prose of the template types by hand: those that are not a figure (a placeholder), code (a code span), a heading, a list number, a step or a section."""
    out, in_code = [], False
    for line in template.splitlines():
        if line.startswith('```'):
            in_code = not in_code
            continue
        if in_code or line.startswith(('#', '|')):
            continue
        line = MARKER.sub('', line)
        line = re.sub(r'`[^`]*`|\{[a-z0-9_]+\}', ' ', line)
        line = re.sub(r'\bsteps? [0-9]+\b|\bsections? [0-9]+\b', ' ', line)
        out += re.findall(r'\d+(?:,\d{3})*(?:\.\d+)?%?', line)
    return [n for n in out if n not in ALLOWED_NUMBERS]


def section(page, heading):
    """The text of a section, from its heading to the next heading of the same or a higher level."""
    level = len(heading) - len(heading.lstrip('#'))
    start = page.index('\n' + heading + '\n') + 1
    ends = [m.start() for m in re.finditer(r'^#{1,%d} ' % level, page[start + len(heading):], flags=re.M)]
    return page[start:start + len(heading) + ends[0]] if ends else page[start:]


def glance_problems(page):
    """The twelve statements of the section "At a glance": each is one line, one sentence, shorter than the limit, under its own heading, in the order of GLANCE."""
    parts = re.split(r'^### (.+)$', section(page, '## At a glance'), flags=re.M)
    names, bodies = parts[1::2], parts[2::2]
    problems = [] if names == GLANCE else [f'the headings are {names}']
    for name, body in zip(names, bodies):
        lines = [x for x in body.splitlines() if x.strip()]
        if len(lines) != 1:
            problems.append(f'{name}: {len(lines)} lines')
            continue
        line = lines[0]
        if len(line) >= G.LIMIT:
            problems.append(f'{name}: {len(line)} characters')
        if not line.endswith('.') or SENTENCES.search(line):
            problems.append(f'{name}: not one sentence')
    return problems


CAVEATS = {                                                              # what the page must say that it does not hold or cannot claim, by statement
    'Label': ['no gene is confirmed free of a safety event'],
    'Test': ['out of fold', 'no split in time', 'no test set held back'],
    'Baseline': ['fits no baseline model'],
    'Not comparable': ['cannot be compared'],
    'How many': ["not in this repo's files"],
    'Cuts': ['cancer only is not run', 'combination does not apply'],
    'Language model only': ['Not run:'],
}
NOT_RUN = ['Result of the all-labeled model', 'Cancer only and not cancer', 'Single agent, combination, mixed', 'Counts of distinct drugs and trials', 'Folds that group paralogs',
           'A split in time or a held-back test', 'Out-of-fold evaluation of the deployed score on drugged genes', 'The DeepSeek V4 Flash only arm']


def caveat_problems(page):
    """The caveats of CAVEATS in their statements, the not-run cells in the table of section 7 with a reason and a need, and the unrun sections under their headings."""
    parts = re.split(r'^### (.+)$', section(page, '## At a glance'), flags=re.M)
    said = dict(zip(parts[1::2], parts[2::2]))
    problems = [f'{name}: "{phrase}" is missing' for name, phrases in CAVEATS.items() for phrase in phrases if phrase not in said.get(name, '')]
    rows = [[c.strip() for c in line.strip('|').split('|')] for line in section(page, '## 7. The cells that were not run').splitlines() if line.startswith('| ')][1:]
    problems += [f'section 7 has no row {cell}' for cell in NOT_RUN if not any(r[0].startswith(cell) for r in rows)]
    problems += [f'section 7: a row has an empty cell: {r[0]}' for r in rows if len(r) != 3 or not all(r)]
    for heading, phrase in (('### Cancer only', 'Cancer only is not run.'), ('### Combination therapy', 'does not apply'), ('## 6. The language model only arm: not run', 'No such arm exists')):
        if phrase not in section(page, heading):
            problems.append(f'{heading}: "{phrase}" is missing')
    return problems


# ---------------------------------------------------------------------------------------------------------------------------------------------
# the page is what the files give
# ---------------------------------------------------------------------------------------------------------------------------------------------
def test_the_page_and_the_numbers_are_what_the_files_give(made):
    assert read(G.DOC) == made['page'], 'docs/GRID.md is stale: run `uv run --with pandas --with numpy --with pyarrow --with scikit-learn python scripts/grid_doc.py`'
    assert read(G.NUMBERS) == made['numbers'], 'docs/grid/numbers.csv is stale: run `uv run --with pandas --with numpy --with pyarrow --with scikit-learn python scripts/grid_doc.py`'
    # the planted bug: a result that moved by a thousandth changes the page and the numbers
    frames = copy.deepcopy(made['F'])
    hit = (frames['results']['reading'] == 'drugged_only') & (frames['results']['stage'] == 'all')
    frames['results'].loc[hit, 'auroc'] += 0.001
    assert G.render(frames, made['K'], made['src']) != made['page'], 'a moved result leaves the page as it was'
    assert G.numbers_csv(frames, made['K']) != made['numbers'], 'a moved result leaves the numbers as they were'


def test_check_mode_fails_on_a_stale_page_or_numbers(made, monkeypatch, tmp_path):
    monkeypatch.setattr(G, 'build', lambda: made)
    assert G.main(['--check']) == 0
    page, numbers = tmp_path / 'GRID.md', tmp_path / 'numbers.csv'
    page.write_text(made['page'], encoding='utf-8')
    numbers.write_text(made['numbers'], encoding='utf-8')
    assert G.main(['--check', '--doc', str(page), '--numbers', str(numbers)]) == 0
    page.write_text(made['page'] + 'x', encoding='utf-8')                                      # the planted bugs: a stale page, a stale CSV, a missing page
    with pytest.raises(SystemExit) as e:
        G.main(['--check', '--doc', str(page), '--numbers', str(numbers)])
    assert 'stale' in str(e.value)
    page.write_text(made['page'], encoding='utf-8')
    changed = made['numbers'].replace(str(made['K']['e1_n']), str(made['K']['e1_n'] + 1), 1)
    assert changed != made['numbers']
    numbers.write_text(changed, encoding='utf-8')
    with pytest.raises(SystemExit):
        G.main(['--check', '--doc', str(page), '--numbers', str(numbers)])
    with pytest.raises(SystemExit):
        G.main(['--check', '--doc', str(tmp_path / 'missing.md'), '--numbers', str(numbers)])


def test_the_command_line_check_passes():
    done = subprocess.run([sys.executable, os.path.join(ROOT, 'scripts', 'grid_doc.py'), '--check'], capture_output=True, text=True, cwd=ROOT)
    assert done.returncode == 0, done.stderr


def test_every_figure_of_the_template_is_made_and_every_one_made_is_used(made):
    template, parts = read(G.TEMPLATE), set(G.parts(made['F'], made['K'], made['src']))
    assert named(template) == parts
    assert named(template + ' {nobody_makes_this}') != parts, 'a figure that nobody makes goes unseen'
    assert named(template.replace('{e1_auc}', '')) != parts, 'a figure that nobody uses goes unseen'


def test_the_template_types_no_figure():
    template = read(G.TEMPLATE)
    assert typed_numbers(template) == []
    assert typed_numbers(template + '\nThe model reads 0.717 on 1,564 genes.\n') == ['0.717', '1,564'], 'a typed number goes unseen'
    assert typed_numbers('A cut of 95% and a null of 0, in `code 12`, after step 3 and in section 4, with {n_claims} claims.\n') == [], 'a fixed number is read as a typed one'


# ---------------------------------------------------------------------------------------------------------------------------------------------
# the rules of the page
# ---------------------------------------------------------------------------------------------------------------------------------------------
def test_prose_holds_one_sentence_a_line(made):
    assert sentence_problems(made['page']) == []
    assert sentence_problems(read(G.TEMPLATE)) == []
    assert sentence_problems(made['page'] + '\nOne sentence here. Another one here.\n') != [], 'a line of two sentences goes unseen'


def test_no_dash_anywhere(made):
    for name, text in (('page', made['page']), ('numbers', made['numbers']), ('template', read(G.TEMPLATE)), ('script', read(os.path.join(ROOT, 'scripts', 'grid_doc.py'))),
                       ('test', read(os.path.join(HERE, 'test_grid_doc.py')))):
        assert dash_problems(text) == [], name
    assert dash_problems('a ' + chr(0x2014) + ' b') and dash_problems('1' + chr(0x2013) + '2'), 'a dash goes unseen'


def test_the_twelve_statements_at_a_glance(made):
    assert glance_problems(made['page']) == []
    page = made['page']
    long_one = re.sub(r'(### Summary\n\n)(.*)', lambda m: m[1] + m[2] + ' ' + 'x' * G.LIMIT, page, count=1)
    assert glance_problems(long_one) != [], 'a line over the limit goes unseen'
    two = re.sub(r'(### Unit\n\n)(.*)', lambda m: m[1] + m[2] + ' Another sentence.', page, count=1)
    assert glance_problems(two) != [], 'two sentences go unseen'
    assert glance_problems(page.replace('### Unit\n\n', '### Units\n\n')) != [], 'a renamed heading goes unseen'


def test_the_page_says_what_it_does_not_hold(made):
    page = made['page']
    assert caveat_problems(page) == []
    assert section(page, '## 6. The language model only arm: not run').startswith('## 6.') and section(page, '### Combination therapy').count('\n#') == 0
    # the planted bugs: a caveat that is dropped, a not-run row that is dropped, a row without a reason, a section that no longer says it
    assert caveat_problems(page.replace("not in this repo's files", 'in this repo')) != [], 'a dropped caveat goes unseen'
    assert caveat_problems(page.replace('| Counts of distinct drugs and trials |', '| Counts |')) != [], 'a dropped not-run row goes unseen'
    assert caveat_problems(page.replace('| a field that names the disease of the trial, or an area for each target |', '|  |')) != [], 'a row without a need goes unseen'
    assert caveat_problems(page.replace('Cancer only is not run.', 'Cancer only is run.')) != [], 'a section that drops its words goes unseen'


def test_the_page_is_about_humanproof_only(made):
    for name, text in (('page', made['page']), ('numbers', made['numbers']), ('template', read(G.TEMPLATE)), ('script', read(os.path.join(ROOT, 'scripts', 'grid_doc.py')))):
        assert not OTHER_PROJECTS.search(text), f'{name} names another project: {OTHER_PROJECTS.search(text)[0]}'
        assert not PROTOCOL_ID.search(text), f'{name} holds a protocol id: {PROTOCOL_ID.search(text)[0]}'
    assert made['page'].startswith('# The grid of experiments, results and features of HumanProof')
    assert OTHER_PROJECTS.search('a sibling repo of ' + 'hyth' + 'erx') and PROTOCOL_ID.search('protocol ' + 'P' + '011'), 'a mention goes unseen'


# ---------------------------------------------------------------------------------------------------------------------------------------------
# the counts
# ---------------------------------------------------------------------------------------------------------------------------------------------
def test_the_inputs_agree(made):
    G.verify(made['src'])
    bad = copy.copy(made['src'])
    bad.oof = made['src'].oof.iloc[1:]
    with pytest.raises(AssertionError):
        G.verify(bad)
    bad = copy.copy(made['src'])
    bad.dr = made['src'].dr.assign(safety_label='unlabeled')
    with pytest.raises(AssertionError):
        G.verify(bad)


def test_the_files_name_no_drug_and_no_trial(made):
    src = made['src']
    assert src.data_files == list(G.DATA_EXPECTED)
    assert G.NAMING.search('drug_name') and G.NAMING.search('nct_id') and G.NAMING.search('indication') and G.NAMING.search('Disease') and not G.NAMING.search('mouseKOScore'), 'a name goes unseen'
    named = {c for df in (src.ot, src.oof, src.dr) for c in df.columns if G.NAMING.search(c)}
    assert named == set(G.EXPLAINED), 'the columns that the page explains are not those that match'
    # the planted bugs: another data file in the repo, a column that names a drug, a list that is not a list of symbols
    bad = copy.copy(src)
    bad.data_files = src.data_files + ['data/drugs.csv']
    with pytest.raises(AssertionError, match='data file'):
        G.verify(bad)
    bad = copy.copy(src)
    bad.dr = src.dr.assign(drug_name='x')
    with pytest.raises(AssertionError, match='names a drug'):
        G.verify(bad)
    bad = copy.copy(src)
    bad.gene_list = src.gene_list + [{'drug': 'x'}]
    with pytest.raises(AssertionError, match='symbols'):
        G.verify(bad)


def test_a_file_is_in_the_repo_when_git_tracks_it(made, tmp_path):
    src, K = made['src'], made['K']
    done = subprocess.run(['git', 'ls-files', '-z'], capture_output=True, text=True, cwd=ROOT, check=True)
    assert src.tracked == {p for p in done.stdout.split('\0') if p}
    assert all(K['present'][p] == (p in src.tracked) for p in K['present'])
    assert all(p in src.tracked for p in G.DATA_EXPECTED), 'the page reads a file that git does not track'
    assert not any(p in src.tracked for p in G.ABSENT), 'a file that the page calls absent is tracked'
    assert set(src.ignored) == set(G.ABSENT), 'a rule of .gitignore does not exclude an absent file'
    # the planted bug: a tracked file is present, and a rule of .gitignore that is gone is seen
    assert G.tracked_files(ROOT) == src.tracked and G.ignored_among(ROOT, []) == set()
    assert G.ignored_among(ROOT, ['data/derived.csv.gz', 'data/safety_model_output/dr/predictions.csv']) == {'data/derived.csv.gz'}
    with pytest.raises(RuntimeError, match='needs a git checkout'):
        G.tracked_files(str(tmp_path))                                                              # an empty folder is no checkout


def test_the_page_is_the_same_in_a_fresh_checkout_and_with_ignored_files(made, tmp_path):
    """Only what git tracks counts: a new repository of the tracked files gives the same page, and files that are untracked or ignored (a database, a model output) change nothing."""
    src = made['src']
    fresh = tmp_path / 'fresh'
    keep = [p for p in sorted(src.tracked) if not p.startswith(G.APART) or p == G.APP_PAGE]
    for p in keep:
        (fresh / p).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(os.path.join(ROOT, p), fresh / p)
    subprocess.run(['git', 'init', '-q'], cwd=fresh, check=True, capture_output=True)
    subprocess.run(['git', 'add', '--'] + keep, cwd=fresh, check=True, capture_output=True)
    extras = [G.MODEL_A_FILES[0], G.INPUT_FILES[0], 'data/untracked_and_not_ignored.csv']
    for p in extras:
        (fresh / p).parent.mkdir(parents=True, exist_ok=True)
        (fresh / p).write_bytes(b'x')
    again = G.Sources(root=str(fresh))
    assert again.tracked == set(keep) and not any(again.present[p] for p in G.ABSENT) and again.data_files == src.data_files
    frames, facts = G.compute(again)
    assert G.render(frames, facts, again) == made['page'] and G.numbers_csv(frames, facts) == made['numbers']
    # the planted bug: the same file, once git tracks it, is in the repo, and the page cannot be made until it is described
    subprocess.run(['git', 'add', '-f', '--', G.INPUT_FILES[0]], cwd=fresh, check=True, capture_output=True)
    tracked = G.Sources(root=str(fresh))
    assert tracked.present[G.INPUT_FILES[0]] and G.INPUT_FILES[0] in tracked.data_files
    with pytest.raises(AssertionError, match='data file'):
        G.compute(tracked)


def test_the_counts_add_up(made):
    F, K = made['F'], made['K']
    labels = F['labels'].set_index('universe')
    for name, r in labels.iterrows():
        assert r['positive'] + r['drugged_safe'] + r['unlabeled'] == r['genes'], name
        assert r['positive_drugged'] + r['positive_no_phase'] == r['positive'], name
    assert labels.loc['Open Targets file', 'genes'] == K['ot_n'] and labels.loc['Drugged-only rows', 'genes'] == K['e1_n'] and labels.loc['DR+PU file', 'genes'] == K['dr_n']
    assert labels.loc['Drugged-only rows', 'positive'] + labels.loc['Drugged-only rows', 'drugged_safe'] == K['e1_n']
    assert labels.loc['DR+PU file', 'positive_drugged'] + labels.loc['DR+PU file', 'drugged_safe'] == K['dr_drugged']
    for universe, sub in F['stage'].groupby('universe'):
        coded, total = sub[~sub['phase'].isin(['all'])], sub[sub['phase'] == 'all'].iloc[0]
        assert coded['genes'].sum() == total['genes'] and coded['positive'].sum() == total['positive'], universe
        assert coded['drugged_safe'].sum() == total['drugged_safe'] and coded['unlabeled'].sum() == total['unlabeled'], universe
        assert (coded['positive'] + coded['drugged_safe'] + coded['unlabeled'] == coded['genes']).all(), universe
    drugged_rows = F['stage'][(F['stage']['universe'] == 'Open Targets file') & F['stage']['phase'].isin([str(c) for c in G.PHASES])]
    assert drugged_rows['genes'].sum() == K['e1_n'], 'the four coded stages are not the genes of the drugged-only model'
    assert F['driver'].iloc[0]['genes'] + F['driver'].iloc[1]['genes'] == F['driver'].iloc[2]['genes'] == K['e1_n']
    assert (F['driver']['flagged'] + F['driver']['not_flagged'] == F['driver']['genes']).all()
    for r in F['setup'].dropna(subset=['fold_train_min']).itertuples():
        assert r.fold_train_min + r.fold_score_max == r.trained and r.fold_train_max + r.fold_score_min == r.trained, r.reading
        assert r.scored >= r.trained and r.fold_score_max - r.fold_score_min <= 2 and abs(r.fold_score_min - r.trained / r.folds) < 2, r.reading
        assert abs(r.base_rate - r.events / r.trained) < 1e-12, r.reading
    dedup = F['dedup']
    assert (dedup['duplicate_ids'].dropna() == 0).all() and (dedup['duplicate_symbols'].dropna() == 0).all()
    ids = dedup.dropna(subset=['unique_ids'])
    assert len(ids) == 3 and (ids['rows'] == ids['unique_ids']).all()
    sym = dedup.dropna(subset=['unique_symbols'])
    assert len(sym) == 3 and (sym['rows_with_symbol'] == sym['unique_symbols']).all()
    assert F['groups']['columns'].sum() == K['n_saved'] == len(F['shap']) and abs(F['groups']['share'].sum() - 1) < 1e-9
    assert len(F['claims']) == K['n_claims'] and (F['claims']['verdict'] == 'differs').sum() == K['n_differ']
    for r in F['results'].itertuples():
        assert 0 < r.events < r.genes and abs(r.base_rate - r.events / r.genes) < 1e-12, (r.reading, r.stage)
        assert r.shown == (r.genes >= G.MIN_GENES and r.events >= G.MIN_CLASS and r.genes - r.events >= G.MIN_CLASS) or r.reading == 'dr_deployed_drugged', (r.reading, r.stage)
        if r.shown:
            assert r.auroc_lo <= r.auroc <= r.auroc_hi and r.auprc_lo <= r.auprc <= r.auprc_hi, (r.reading, r.stage)
            assert r.draws_used >= G.DRAWS - 5, (r.reading, r.stage)
    p = F['paired'].iloc[0]
    assert p['auroc_diff_lo'] <= p['auroc_diff'] <= p['auroc_diff_hi'] and p['auprc_diff_lo'] <= p['auprc_diff'] <= p['auprc_diff_hi']
    assert abs(p['auroc_diff'] - (p['auroc_2'] - p['auroc_1'])) < 1e-12 and p['genes'] == K['shared_n']


def test_the_results_are_those_of_scikit_learn(made):
    src, res = made['src'], made['F']['results'].set_index(['reading', 'stage'])
    oof = src.oof
    assert abs(res.loc[('drugged_only', 'all'), 'auroc'] - roc_auc_score(oof['true_label'], oof['safety_score_oof'])) < 1e-12
    assert abs(res.loc[('drugged_only', 'all'), 'auprc'] - average_precision_score(oof['true_label'], oof['safety_score_oof'])) < 1e-12
    dr = src.dr
    d = dr[dr['is_drugged'] == 1]
    y = (d['hasSafetyEvent'] == -1.0).astype(int)
    assert abs(res.loc[('dr_step3', 'all'), 'auroc'] - roc_auc_score(y, d['m_hat'])) < 1e-12
    assert abs(res.loc[('dr_step3', 'all'), 'auprc'] - average_precision_score(y, d['m_hat'])) < 1e-12
    assert abs(res.loc[('dr_step1', 'all'), 'auroc'] - roc_auc_score(dr['is_drugged'], dr['pi_hat'])) < 1e-12
    u = dr[dr['is_drugged'] == 0]
    yu = (u['hasSafetyEvent'] == -1.0).astype(int)
    assert abs(res.loc[('dr_und_deployed', 'all'), 'auroc'] - roc_auc_score(yu, u['safety_score_dr'])) < 1e-12
    assert abs(res.loc[('dr_und_outcome', 'all'), 'auprc'] - average_precision_score(yu, u['m_hat'])) < 1e-12
    stage = oof[oof['maxClinicalTrialPhase'] == 1.0]
    assert abs(res.loc[('drugged_only', '1.0'), 'auroc'] - roc_auc_score(stage['true_label'], stage['safety_score_oof'])) < 1e-12


def test_the_bootstrap_is_scikit_learn_on_every_draw():
    rng = np.random.default_rng(3)
    n = 400
    y = (rng.random(n) < 0.3).astype(int)
    s = np.round(rng.normal(size=n) + 0.8 * y, 1)                              # rounded, so the scores tie
    assert len(np.unique(s)) < n
    a, p = G.point(y, s)
    assert abs(a - roc_auc_score(y, s)) < 1e-12 and abs(p - average_precision_score(y, s)) < 1e-12
    order, ys, starts = G._prepare(y, s)
    for counts in G._draws(n, 20, 5):
        auc, ap = G._metrics(ys, starts, counts[:, order].astype(float))
        for row in range(3):
            idx = np.repeat(np.arange(n), counts[row])
            assert abs(auc[row] - roc_auc_score(y[idx], s[idx])) < 1e-12
            assert abs(ap[row] - average_precision_score(y[idx], s[idx])) < 1e-12
        assert not np.isnan(auc).any() and not np.isnan(ap).any()
    # the planted bug: a score that is not the one of the draw gives another number
    assert abs(G.point(y, -s)[0] - roc_auc_score(y, s)) > 0.1
    # a draw with no flagged gene gives nan and not a number
    auc, ap = G._metrics(ys, starts, np.where(ys == 1, 0.0, 1.0)[None, :])
    assert np.isnan(auc).all() and np.isnan(ap).all()


def test_the_bootstrap_is_deterministic_and_holds_its_point(made):
    oof = made['src'].oof
    y, s = oof['true_label'].values, oof['safety_score_oof'].values
    one, two = G.bootstrap(y, s, draws=200), G.bootstrap(y, s, draws=200)
    assert one == two, 'the same seed gives two intervals'
    assert G.bootstrap(y, s, draws=200, seed=G.SEED + 1) != one, 'another seed gives the same interval'
    assert one['auroc_lo'] < one['auroc'] < one['auroc_hi'] and one['auprc_lo'] < one['auprc'] < one['auprc_hi']
    pair = G.bootstrap_pair(y, s, s, draws=200)
    assert pair['auroc_diff'] == 0 and pair['auroc_diff_lo'] == 0 and pair['auroc_diff_hi'] == 0, 'a model against itself differs'


def test_fold_sizes_follow_the_class_counts():
    # by hand: 10 flagged genes give 2 to each of 5 folds, and 93 others give 18 to two folds and 19 to three, so a fold scores 20 or 21 genes and trains on 82 or 83
    assert G.fold_sizes(10, 93, 5, 42) == (82, 83, 20, 21)
    assert G.fold_sizes(10, 93, 5, 1) == (82, 83, 20, 21), 'the seed changes a fold size'
    assert G.fold_sizes(10, 90, 5, 42) == (80, 80, 20, 20)
    assert G.fold_sizes(10, 90, 2, 42) == (50, 50, 50, 50)


# ---------------------------------------------------------------------------------------------------------------------------------------------
# the claims of the text
# ---------------------------------------------------------------------------------------------------------------------------------------------
def test_a_claim_the_files_do_not_hold_stops_the_maker(made):
    G.values(made['F'], made['K'], made['src'])
    frames, facts = copy.deepcopy(made['F']), dict(made['K'])
    facts['min_flag_score'] = facts['max_safe_score'] - 0.1
    with pytest.raises(AssertionError):
        G.values(frames, facts, made['src'])
    frames = copy.deepcopy(made['F'])
    frames['paired'].loc[0, 'auroc_diff_lo'] = 0.01
    with pytest.raises(AssertionError):
        G.values(frames, dict(made['K']), made['src'])
    facts = dict(made['K'])
    facts['shap_rho'] = 0.0
    with pytest.raises(AssertionError):
        G.values(copy.deepcopy(made['F']), facts, made['src'])
    # the app page gives for the fixed prior an AUROC that is not the saved one
    K, src = made['K'], copy.copy(made['src'])
    page = [x.replace(f"from {K['shift_from']} to {K['shift_to']}", f"from {K['shift_from']} to {K['shift_from']}") for x in made['src'].lines[G.APP_PAGE]]
    assert page != made['src'].lines[G.APP_PAGE]
    src.lines = {**made['src'].lines, G.APP_PAGE: page}
    with pytest.raises(AssertionError, match='AUROC that the app page gives'):
        G.compute_claims(src, copy.deepcopy(made['F']), dict(K))
    # the script no longer clips the deployed score, and the comment on the prior is gone
    texts = {n: made['src'].scripts[f'data/{n}'] for n in ('safety_model.py', 'safety_model_dr.py')}
    G.code_facts(texts)
    with pytest.raises(AssertionError, match='clipped'):
        G.code_facts({**texts, 'safety_model_dr.py': texts['safety_model_dr.py'].replace('np.clip(final_model.predict(', 'np.abs(final_model.predict(')})
    with pytest.raises(AssertionError, match='monotone shift'):
        G.code_facts({**texts, 'safety_model_dr.py': texts['safety_model_dr.py'].replace('monotone logit-shift', 'shift')})


def test_the_claims_of_the_summary_hold(made):
    F, K, src = made['F'], made['K'], made['src']
    res = F['results'].set_index(['reading', 'stage'])
    assert res.loc[('dr_deployed_drugged', 'all'), 'auroc'] == 1.0 and K['min_flag_score'] > K['max_safe_score']
    assert K['ot_share'] > 0.9 and K['shap_rho'] > 0.5 > abs(K['shap_rho_shuffled'])
    assert K['n_date_cols'] == 0 and K['textbook_n'] == len(G.TEXTBOOK) and 0 <= K['textbook_flagged'] <= K['textbook_n']
    assert K['pi'] == src.facts['pi_override'] and K['n_differ'] > 0
    assert abs(float(K['shift_to']) - res.loc[('dr_step3', 'all'), 'auroc']) < 0.0015 and float(K['shift_from']) < float(K['shift_to'])
    assert K['splitters'] == ['StratifiedKFold'] and K['learners'] == ['XGBClassifier', 'XGBRegressor']
    assert not any(K['present'][p] for p in G.ABSENT) and all(K['present'][p] for p in G.DR_PICTURES)
    lo, hi = src.facts['clip_range']
    assert lo < hi and src.dr['safety_score_dr'].min() == lo and src.dr['safety_score_dr'].max() <= hi
    und = src.dr[src.dr['is_drugged'] == 0]
    assert (und['Y_tilde'] == und['m_hat']).all() and K['und_zero'] == int((und['safety_score_dr'] == lo).sum()) > 0
    cv = src.facts['e1_cv']
    assert cv['shuffle'] is True and cv['random_state'] == G.SEED and cv['n_splits'] == K['folds']
    release_end = pd.Timestamp(K['ot_month']) + pd.offsets.MonthEnd(0)
    assert release_end < pd.Timestamp(G.RELEASE), 'the Open Targets release ends after the release date of the language model'


def test_the_numbers_csv_holds_the_tables(made):
    rows = list(csv.reader(io.StringIO(made['numbers'])))
    assert rows[0] == ['table', 'row', 'field', 'value'] and all(len(r) == 4 for r in rows)
    tables = {r[0] for r in rows[1:]}
    assert set(made['F']) | {'facts'} == tables
    got = {(r[0], r[1], r[2]): r[3] for r in rows[1:]}
    assert len(got) == len(rows) - 1, 'two values share one key'
    oof, K = made['src'].oof, made['K']
    assert abs(float(got[('results', 'drugged_only|all', 'auroc')]) - roc_auc_score(oof['true_label'], oof['safety_score_oof'])) < 1e-6
    assert got[('results', 'drugged_only|all', 'events')] == str(K['e1_events']) and got[('facts', '', 'dr_drugged')] == str(K['dr_drugged'])
    assert got[('labels', 'Open Targets file', 'genes')] == str(K['ot_n'])
    assert not any(r[3] == 'nan' for r in rows[1:]), 'a nan reaches the CSV'
    # the planted bug: a table that changes, changes the CSV
    frames = copy.deepcopy(made['F'])
    frames['labels'].loc[0, 'genes'] += 1
    assert G.numbers_csv(frames, K) != made['numbers']
