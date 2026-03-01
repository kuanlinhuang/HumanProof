"""Load HumanProof database with REAL data from all three sources:

  1. Expression:  data/cellxgene/ -- CZ CellxGene census aggregated per cell type
                  (log1p mean expression + fraction expressing)
  2. Dosage:      data/LOEUF_scores.csv.gz -- LOEUF constraint scores
                  (https://grr.iossifovlab.com/gene_properties/gene_scores/LOEUF/)
  3. pLoF:        data/genebass_pLoF_filtered.pkl -- Genebass burden-test results
                  (same real data source as before)

Run from the backend/ directory:
    python load_real_data.py
"""

import json
import math
import os
import pickle
import sqlite3
import sys
from pathlib import Path

import numpy as np
import pandas as pd

def _resolve_data_dir() -> Path:
    env_data_dir = os.getenv("HUMANPROOF_DATA_DIR")
    if env_data_dir:
        return Path(env_data_dir)

    script_dir = Path(__file__).resolve().parent
    candidates = [
        script_dir.parent / "data",  # local: <repo>/backend/load_real_data.py -> <repo>/data
        script_dir / "data",         # Railway: /app/load_real_data.py -> /app/data
        Path("/app/data"),
    ]
    for path in candidates:
        if path.exists():
            return path

    return script_dir.parent / "data"


DATA_DIR = _resolve_data_dir()
CELLXGENE_DIR = DATA_DIR / "cellxgene"
PLOF_PICKLE = DATA_DIR / "genebass_pLoF_filtered.pkl"
LOEUF_FILE = DATA_DIR / "LOEUF_scores.csv.gz"

# ────────────────────────────────────────────────────────────────────────────
# Gene list: kept only for reference (used by HumanProof ML model / SHAP export)
# The DB is now built for ALL ~19K LOEUF protein-coding genes.
# ────────────────────────────────────────────────────────────────────────────
GENES = {
    "TP53": "ENSG00000141510", "BRCA1": "ENSG00000012048", "BRCA2": "ENSG00000139618",
    "EGFR": "ENSG00000146648", "ERBB2": "ENSG00000141736", "KRAS": "ENSG00000133703",
    "BRAF": "ENSG00000157764", "PIK3CA": "ENSG00000121879", "PTEN": "ENSG00000171862",
    "MYC": "ENSG00000136997", "RB1": "ENSG00000139687", "APC": "ENSG00000134982",
    "VHL": "ENSG00000134086", "WT1": "ENSG00000184937", "NF1": "ENSG00000196712",
    "NF2": "ENSG00000186575", "RET": "ENSG00000165731", "MET": "ENSG00000105976",
    "ALK": "ENSG00000171094", "FGFR1": "ENSG00000077782", "FGFR2": "ENSG00000066468",
    "FGFR3": "ENSG00000068078", "PDGFRA": "ENSG00000134853", "KIT": "ENSG00000157404",
    "ABL1": "ENSG00000097007", "JAK2": "ENSG00000096968", "CDK4": "ENSG00000135446",
    "CDK6": "ENSG00000105810", "CCND1": "ENSG00000110092", "MDM2": "ENSG00000135679",
    "BCL2": "ENSG00000171791", "BAX": "ENSG00000087088", "CASP3": "ENSG00000164305",
    "TNF": "ENSG00000232810", "IL6": "ENSG00000136244", "IL1B": "ENSG00000125538",
    "IL2": "ENSG00000109471", "IL10": "ENSG00000136634", "IFNG": "ENSG00000111537",
    "TGFB1": "ENSG00000105329", "VEGFA": "ENSG00000112715", "FLT1": "ENSG00000102755",
    "KDR": "ENSG00000128052", "PDCD1": "ENSG00000188389", "CD274": "ENSG00000120217",
    "CTLA4": "ENSG00000163599", "CD28": "ENSG00000178562", "FOXP3": "ENSG00000049768",
    "CD4": "ENSG00000010610", "CD8A": "ENSG00000153563", "HLA-A": "ENSG00000206503",
    "B2M": "ENSG00000166710", "PCNA": "ENSG00000132646", "TOP2A": "ENSG00000131747",
    "ESR1": "ENSG00000091831", "AR": "ENSG00000169083", "ERBB3": "ENSG00000065361",
    "ERBB4": "ENSG00000178568", "IGF1R": "ENSG00000140443", "INSR": "ENSG00000171105",
    "SRC": "ENSG00000197122", "RAF1": "ENSG00000132155", "MAP2K1": "ENSG00000169032",
    "MAPK1": "ENSG00000100030", "MAPK3": "ENSG00000102882", "AKT1": "ENSG00000142208",
    "MTOR": "ENSG00000198793", "NOTCH1": "ENSG00000148400", "WNT1": "ENSG00000125084",
    "SHH": "ENSG00000164690", "PTCH1": "ENSG00000185920", "SMO": "ENSG00000128602",
    "STAT3": "ENSG00000168610", "NFKB1": "ENSG00000109320", "RELA": "ENSG00000173039",
    "CREBBP": "ENSG00000005339", "EP300": "ENSG00000100393", "HDAC1": "ENSG00000116478",
    "DNMT1": "ENSG00000130816", "TET2": "ENSG00000168769", "IDH1": "ENSG00000138413",
    "IDH2": "ENSG00000182054", "ARID1A": "ENSG00000117713", "SMAD4": "ENSG00000141646",
    "CDKN2A": "ENSG00000147889", "CDKN1A": "ENSG00000124762", "GAPDH": "ENSG00000111640",
    "ACTB": "ENSG00000075624", "ALB": "ENSG00000163631", "INS": "ENSG00000254647",
    "GCG": "ENSG00000115263", "EPO": "ENSG00000130427", "THPO": "ENSG00000090534",
    "CSF3": "ENSG00000108342", "APOE": "ENSG00000130203", "APP": "ENSG00000142192",
    "MAPT": "ENSG00000186868", "SNCA": "ENSG00000145335", "HTT": "ENSG00000197386",
    "SOD1": "ENSG00000142168", "CFTR": "ENSG00000001626", "DMD": "ENSG00000198947",
    "F8": "ENSG00000185010", "HBB": "ENSG00000244734", "HBA1": "ENSG00000206172",
    "FGA": "ENSG00000171560", "CYP3A4": "ENSG00000160868", "CYP2D6": "ENSG00000100197",
}


# ────────────────────────────────────────────────────────────────────────────
# Gene universe builder
# ────────────────────────────────────────────────────────────────────────────

def build_gene_universe() -> dict[str, tuple[str, str]]:
    """Build {gene_symbol: (ensembl_id, col_name_in_matrix)} for all LOEUF
    protein-coding genes that have CellxGene expression data.

    Uses LOEUF_scores.csv.gz as the master protein-coding gene list (~19,704 genes)
    and gene_metadata.csv to map column names → Ensembl IDs.
    """
    print("  Loading LOEUF gene list (master protein-coding universe)...")
    loeuf_df = pd.read_csv(LOEUF_FILE).dropna(subset=["LOEUF"])
    loeuf_genes: set[str] = set(loeuf_df["gene"])  # duplicates resolved by set
    print(f"    {len(loeuf_genes):,} protein-coding genes in LOEUF (NaN LOEUF excluded)")

    print("  Loading CellxGene gene_metadata (col_name → Ensembl ID)...")
    gm = pd.read_csv(CELLXGENE_DIR / "gene_metadata.csv")
    # gene_metadata columns: gene_id (Ensembl), gene_symbol (matrix column name)
    # For plain genes: gene_symbol == gene_symbol (e.g. "EGFR")
    # For multi-ENSG genes: gene_symbol == "GENE_ENSGXXXXXXXX" (e.g. "F8_ENSG00000185010")

    gene_universe: dict[str, tuple[str, str]] = {}
    for col_name, ensembl_id in zip(gm["gene_symbol"], gm["gene_id"]):
        # Parse gene symbol from column name
        if "_ENSG" in col_name:
            gene_sym = col_name.split("_ENSG")[0]
        else:
            gene_sym = col_name

        # Only include LOEUF protein-coding genes
        if gene_sym not in loeuf_genes:
            continue

        # First occurrence wins (avoids duplicates for multi-ENSG genes)
        if gene_sym not in gene_universe:
            gene_universe[gene_sym] = (ensembl_id, col_name)

    print(f"    {len(gene_universe):,} genes have CellxGene expression data "
          f"(out of {len(loeuf_genes):,} LOEUF genes)")
    return gene_universe

# ────────────────────────────────────────────────────────────────────────────
# CellxGene cell type → (display_label, tissue, organ) mapping
#
# Keys are exact cell_type_label values from celltype_metadata.csv.
# Selected for biological relevance, data quality (total_cells ≥ 5k),
# and diversity across tissue types.
# ────────────────────────────────────────────────────────────────────────────
CELLXGENE_TISSUE_MAP: dict[str, tuple[str, str, str]] = {
    # ── Brain ─────────────────────────────────────────────────────────────
    "neuron":                                   ("Neurons",                 "Brain",       "Brain"),
    "glutamatergic neuron":                     ("Glutamatergic neurons",   "Brain",       "Brain"),
    "GABAergic neuron":                         ("GABAergic neurons",       "Brain",       "Brain"),
    "interneuron":                              ("Interneurons",            "Brain",       "Brain"),
    "oligodendrocyte":                          ("Oligodendrocytes",        "Brain",       "Brain"),
    "astrocyte":                                ("Astrocytes",              "Brain",       "Brain"),
    "microglial cell":                          ("Microglia",               "Brain",       "Brain"),
    # ── Heart ─────────────────────────────────────────────────────────────
    "cardiac muscle cell":                      ("Cardiomyocytes",          "Heart",       "Heart"),
    "cardiac endothelial cell":                 ("Cardiac endothelial",     "Heart",       "Heart"),
    "fibroblast of cardiac tissue":             ("Cardiac fibroblasts",     "Heart",       "Heart"),
    # ── Liver ─────────────────────────────────────────────────────────────
    "hepatocyte":                               ("Hepatocytes",             "Liver",       "Liver"),
    "Kupffer cell":                             ("Kupffer cells",           "Liver",       "Liver"),
    "hepatic stellate cell":                    ("Hepatic stellate cells",  "Liver",       "Liver"),
    "cholangiocyte":                            ("Cholangiocytes",          "Liver",       "Liver"),
    # ── Lung ──────────────────────────────────────────────────────────────
    "type I pneumocyte":                        ("Type I pneumocytes",      "Lung",        "Lung"),
    "type II pneumocyte":                       ("Type II pneumocytes",     "Lung",        "Lung"),
    "alveolar macrophage":                      ("Alveolar macrophages",    "Lung",        "Lung"),
    "ciliated epithelial cell":                 ("Ciliated epithelium",     "Lung",        "Lung"),
    "club cell":                                ("Club cells",              "Lung",        "Lung"),
    "lung macrophage":                          ("Lung macrophages",        "Lung",        "Lung"),
    # ── Kidney ────────────────────────────────────────────────────────────
    "epithelial cell of proximal tubule":       ("Proximal tubule",         "Kidney",      "Kidney"),
    "kidney tubule cell":                       ("Tubular cells",           "Kidney",      "Kidney"),
    "podocyte":                                 ("Podocytes",               "Kidney",      "Kidney"),
    "kidney collecting duct cell":              ("Collecting duct",         "Kidney",      "Kidney"),
    "epithelial cell of distal tubule":         ("Distal tubule",           "Kidney",      "Kidney"),
    # ── Intestine ─────────────────────────────────────────────────────────
    "enterocyte":                               ("Enterocytes",             "Intestine",   "Intestine"),
    "goblet cell":                              ("Goblet cells",            "Intestine",   "Intestine"),
    "enteroendocrine cell":                     ("Enteroendocrine cells",   "Intestine",   "Intestine"),
    "intestinal crypt stem cell":               ("Intestinal stem cells",   "Intestine",   "Intestine"),
    # ── Bone Marrow ───────────────────────────────────────────────────────
    "hematopoietic stem cell":                  ("HSCs",                    "Bone marrow", "Bone marrow"),
    "erythroblast":                             ("Erythroblasts",           "Bone marrow", "Bone marrow"),
    "erythrocyte":                              ("Erythrocytes",            "Bone marrow", "Bone marrow"),
    "megakaryocyte":                            ("Megakaryocytes",          "Bone marrow", "Bone marrow"),
    "neutrophil":                               ("Neutrophils",             "Bone marrow", "Bone marrow"),
    "plasma cell":                              ("Plasma cells",            "Bone marrow", "Bone marrow"),
    # ── Lymph Node / Immune ───────────────────────────────────────────────
    "T cell":                                   ("T cells",                 "Lymph node",  "Lymph node"),
    "CD4-positive, alpha-beta T cell":          ("CD4+ T cells",            "Lymph node",  "Lymph node"),
    "CD8-positive, alpha-beta T cell":          ("CD8+ T cells",            "Lymph node",  "Lymph node"),
    "regulatory T cell":                        ("Regulatory T cells",      "Lymph node",  "Lymph node"),
    "B cell":                                   ("B cells",                 "Lymph node",  "Lymph node"),
    "natural killer cell":                      ("NK cells",                "Lymph node",  "Lymph node"),
    "dendritic cell":                           ("Dendritic cells",         "Lymph node",  "Lymph node"),
    # ── Spleen ────────────────────────────────────────────────────────────
    "macrophage":                               ("Macrophages",             "Spleen",      "Spleen"),
    "monocyte":                                 ("Monocytes",               "Spleen",      "Spleen"),
    # ── Skin ──────────────────────────────────────────────────────────────
    "keratinocyte":                             ("Keratinocytes",           "Skin",        "Skin"),
    "melanocyte":                               ("Melanocytes",             "Skin",        "Skin"),
    "skin fibroblast":                          ("Skin fibroblasts",        "Skin",        "Skin"),
    "fibroblast of dermis":                     ("Dermal fibroblasts",      "Skin",        "Skin"),
    # ── Pancreas ──────────────────────────────────────────────────────────
    "type B pancreatic cell":                   ("Beta cells",              "Pancreas",    "Pancreas"),
    "pancreatic A cell":                        ("Alpha cells",             "Pancreas",    "Pancreas"),
    "acinar cell":                              ("Acinar cells",            "Pancreas",    "Pancreas"),
    "pancreatic ductal cell":                   ("Ductal cells",            "Pancreas",    "Pancreas"),
    # ── Muscle ────────────────────────────────────────────────────────────
    "skeletal muscle fiber":                    ("Skeletal muscle fibers",  "Muscle",      "Muscle"),
    "cell of skeletal muscle":                  ("Skeletal muscle cells",   "Muscle",      "Muscle"),
    "skeletal muscle satellite cell":           ("Satellite cells",         "Muscle",      "Muscle"),
    # ── Adrenal ───────────────────────────────────────────────────────────
    "cortical cell of adrenal gland":           ("Adrenal cortical cells",  "Adrenal",     "Adrenal"),
    "chromaffin cell":                          ("Chromaffin cells",        "Adrenal",     "Adrenal"),
    # ── Breast ────────────────────────────────────────────────────────────
    "mammary gland epithelial cell":            ("Mammary epithelial",      "Breast",      "Breast"),
    "luminal epithelial cell of mammary gland": ("Luminal epithelial",      "Breast",      "Breast"),
    "fibroblast of mammary gland":              ("Mammary fibroblasts",     "Breast",      "Breast"),
}


# ────────────────────────────────────────────────────────────────────────────
# Phenotype categorizer (unchanged from original)
# ────────────────────────────────────────────────────────────────────────────

PHENOTYPE_RULES = [
    (["cancer", "carcinoma", "melanoma", "lymphoma", "leukemia", "leukaemia",
      "neoplasm", "tumor", "tumour", "malignant", "myeloma", "sarcoma",
      "mesothelioma", "glioma"], "cancer", "Cancer"),
    (["heart", "cardiac", "coronary", "myocardial", "atrial", "ventricular",
      "angina", "arrhythmia", "aortic", "hypertension", "blood pressure",
      "pulse rate", "heart rate", "systolic", "diastolic", "vascular",
      "peripheral vascular", "arterial", "stroke", "cerebrovascular",
      "thrombosis", "embolism", "aneurysm", "ischaemic heart",
      "cardiomyopathy", "valve", "ecg", "electrocardiogr"], "cardiovascular", "Cardiovascular"),
    (["haemoglobin", "hemoglobin", "platelet", "erythrocyte", "leukocyte",
      "lymphocyte", "monocyte", "neutrophil", "eosinophil", "basophil",
      "reticulocyte", "haematocrit", "hematocrit", "red blood cell",
      "white blood cell", "mean corpuscular", "mcv", "mch", "mchc",
      "red cell", "blood count", "anaemia", "anemia", "bleeding",
      "coagulation", "clotting", "fibrinogen", "von willebrand",
      "sickle cell", "thalassaemia", "thalassemia", "haemophilia",
      "hemophilia", "iron deficiency"], "hematologic", "Hematologic"),
    (["alzheimer", "parkinson", "dementia", "epilepsy", "seizure",
      "multiple sclerosis", "neuropathy", "migraine", "headache",
      "neuralgia", "brain", "cerebral", "meningitis", "encephalitis",
      "depression", "anxiety", "schizophrenia", "bipolar", "psychosis",
      "psychiatric", "insomnia", "sleep", "neurodegenerative", "motor neuron",
      "als", "huntington", "tremor", "cognitive", "intelligence",
      "reaction time", "prospective memory", "fluid intelligence",
      "nerve", "mental health"], "neurological", "Neurological"),
    (["diabetes", "glucose", "insulin", "hba1c", "glycated", "cholesterol",
      "triglyceride", "lipid", "hdl", "ldl", "apolipoprotein", "bmi",
      "body mass index", "obesity", "waist", "hip circumference", "weight",
      "body fat", "impedance", "basal metabolic", "metabolic", "gout",
      "urate", "uric acid", "adiposity", "fat mass", "lean mass",
      "thyroid", "tsh", "thyroxine", "testosterone", "oestradiol",
      "estradiol", "shbg", "igf-1", "growth hormone", "cortisol",
      "vitamin d", "calcium", "phosphate", "parathyroid",
      "glycoprotein acetyls", "fatty acid"], "metabolic", "Metabolic"),
    (["liver", "hepat", "bilirubin", "albumin", "alt", "ast", "ggt",
      "gamma-glutamyl", "alkaline phosphatase", "cirrhosis", "fatty liver",
      "gallstone", "gallbladder", "cholecyst", "jaundice", "bile",
      "cholelithiasis"], "hepatic", "Hepatic"),
    (["kidney", "renal", "creatinine", "cystatin", "glomerular",
      "nephr", "urea", "urinary", "urine", "bladder", "dialysis",
      "microalbumin", "proteinuria", "kidney stone"], "renal", "Renal"),
    (["lung", "pulmonary", "respiratory", "asthma", "copd", "pneumonia",
      "bronch", "fev1", "fvc", "peak expiratory", "spirometry",
      "emphysema", "fibrosis", "pleural", "sleep apn", "oxygen sat",
      "breathless"], "respiratory", "Respiratory"),
    (["bone", "fracture", "osteopor", "osteoarth", "arthritis", "joint",
      "muscle", "muscular", "dystrophy", "back pain", "spine", "spinal",
      "scoliosis", "grip strength", "hand grip", "sarcopenia",
      "rheumatoid", "gout", "height", "standing height", "sitting height",
      "heel bone", "bone density", "mineral density"], "musculoskeletal", "Musculoskeletal"),
    (["immune", "autoimmune", "lupus", "psoriasis", "crohn", "colitis",
      "inflammatory bowel", "celiac", "coeliac", "sjogren", "vasculitis",
      "sarcoidosis", "allergy", "allergic", "atopic", "eczema", "dermatitis",
      "urticaria", "hay fever", "c-reactive protein", "crp",
      "immunoglobulin", "igg", "ige", "complement"], "immunologic", "Immunologic"),
    (["eye", "ocular", "retinal", "glaucoma", "cataract", "macular",
      "myopia", "hypermetropia", "astigmatism", "intraocular pressure",
      "corneal", "visual acuity", "refractive error", "optic",
      "blindness", "spectacle", "glasses"], "ophthalmologic", "Ophthalmologic"),
    (["skin", "dermat", "acne", "psoriasis", "eczema", "melanoma",
      "hair", "alopecia", "nail", "pigment", "sunburn",
      "skin colour", "skin color"], "dermatologic", "Dermatologic"),
    (["prostate", "ovarian", "uterine", "cervical", "breast",
      "endometri", "menstrua", "menopause", "fertility", "pregnancy",
      "birth weight", "gestational", "miscarriage", "contraceptive",
      "erectile", "impotence", "testicular"], "reproductive", "Reproductive"),
    (["gastric", "stomach", "intestin", "bowel", "colon", "rectal",
      "colorectal", "oesophag", "esophag", "reflux", "hernia",
      "appendic", "diverticul", "irritable bowel", "ibs",
      "constipation", "diarrhoea", "diarrhea", "peptic", "ulcer",
      "pancreat", "abdominal"], "gastrointestinal", "Gastrointestinal"),
    (["hearing", "tinnitus", "ear", "deafness", "otitis",
      "noise exposure", "cochlear"], "audiologic", "Audiologic"),
    (["arm span", "trunk", "leg", "body size", "birth weight",
      "comparative body size", "handedness"], "anthropometric", "Anthropometric"),
]

_CATEGORY_CACHE: dict[str, tuple[str, str]] = {}


def classify_phenotype(description: str) -> tuple[str, str]:
    if description in _CATEGORY_CACHE:
        return _CATEGORY_CACHE[description]
    desc_lower = description.lower()
    for keywords, category, organ in PHENOTYPE_RULES:
        for kw in keywords:
            if kw in desc_lower:
                _CATEGORY_CACHE[description] = (category, organ)
                return category, organ
    if any(kw in desc_lower for kw in ["medication", "treatment", "prescribed", "tablet"]):
        _CATEGORY_CACHE[description] = ("medication", "Medication")
        return "medication", "Medication"
    if any(kw in desc_lower for kw in ["operation", "surgery", "procedure", "operative"]):
        _CATEGORY_CACHE[description] = ("procedural", "Procedural")
        return "procedural", "Procedural"
    _CATEGORY_CACHE[description] = ("other", "Other")
    return "other", "Other"


# ────────────────────────────────────────────────────────────────────────────
# Real data loaders
# ────────────────────────────────────────────────────────────────────────────

def load_real_expression_data(gene_universe: dict[str, tuple[str, str]]) -> list[dict]:
    """Load real cell-type expression from CZ CellxGene census aggregations.

    Source files (all in data/cellxgene/):
      - celltype_log1p_mean_expression.csv.gz  -- log1p(mean raw count) per cell type/gene
      - celltype_fraction_expressing.csv.gz    -- fraction of cells with count > 0
      - celltype_metadata.csv                  -- total_cells per cell type

    gene_universe: {gene_symbol: (ensembl_id, col_name_in_matrix)}

    Returns long-format records ready for SQLite insertion.
    """
    print("  Loading CellxGene cell type metadata...")
    meta = pd.read_csv(CELLXGENE_DIR / "celltype_metadata.csv")
    n_cells_map = dict(zip(meta["cell_type_label"], meta["total_cells"].astype(int)))

    log1p_path = CELLXGENE_DIR / "celltype_log1p_mean_expression.csv"
    frac_path  = CELLXGENE_DIR / "celltype_fraction_expressing.csv"

    # Build col_name → gene_symbol reverse map; filter to genes in universe
    col_to_gene: dict[str, str] = {col: gene for gene, (_, col) in gene_universe.items()}
    cols_to_load = ["cell_type_label"] + list(col_to_gene.keys())
    selected_ct = list(CELLXGENE_TISSUE_MAP.keys())

    print(f"  Loading log1p mean expression for {len(gene_universe):,} genes "
          f"× {len(selected_ct)} cell types...")
    log1p_df = pd.read_csv(
        log1p_path,
        compression="gzip",
        usecols=cols_to_load,
        index_col="cell_type_label",
    )
    log1p_df = log1p_df.loc[log1p_df.index.isin(selected_ct)]

    print("  Loading fraction expressing...")
    frac_df = pd.read_csv(
        frac_path,
        compression="gzip",
        usecols=cols_to_load,
        index_col="cell_type_label",
    )
    frac_df = frac_df.loc[frac_df.index.isin(selected_ct)]

    # Rename matrix column names → clean gene symbols
    log1p_df.rename(columns=col_to_gene, inplace=True)
    frac_df.rename(columns=col_to_gene, inplace=True)

    print("  Converting to long format...")
    records = []
    for ct_label, (display_label, tissue, organ) in CELLXGENE_TISSUE_MAP.items():
        if ct_label not in log1p_df.index:
            print(f"    WARNING: '{ct_label}' not found in expression matrix, skipping")
            continue

        n_cells = n_cells_map.get(ct_label, 0)
        for gene, (ensembl_id, _) in gene_universe.items():
            if gene not in log1p_df.columns:
                continue
            mean_expr = float(log1p_df.at[ct_label, gene])
            pct_expr  = float(frac_df.at[ct_label, gene]) if ct_label in frac_df.index else 0.0

            # Sanitize NaN/Inf
            if not math.isfinite(mean_expr):
                mean_expr = 0.0
            if not math.isfinite(pct_expr):
                pct_expr = 0.0
            pct_expr = max(0.0, min(1.0, pct_expr))

            records.append({
                "gene_symbol":     gene,
                "ensembl_id":      ensembl_id,
                "cell_type":       display_label,
                "tissue":          tissue,
                "organ":           organ,
                "mean_expression": round(mean_expr, 5),
                "pct_expressed":   round(pct_expr, 5),
                "n_cells":         n_cells,
            })

    print(f"  Prepared {len(records):,} real expression records "
          f"({len(gene_universe):,} genes × {len(CELLXGENE_TISSUE_MAP)} cell types)")
    return records


def load_real_plof_data(top_n_per_gene: int = 200) -> list[dict]:
    """Load real Genebass pLoF burden-test data for ALL genes.

    Keeps the top_n_per_gene associations per gene by p-value to keep DB size
    manageable (~18K genes × 200 = ~3.6M rows max).
    """
    print(f"  Loading {PLOF_PICKLE} ...")
    df = pd.read_pickle(PLOF_PICKLE)
    print(f"  Full dataset: {len(df):,} rows, {df['gene'].nunique():,} genes")

    df = df.dropna(subset=["Pvalue"]).copy()

    # Keep top N associations per gene by ascending p-value
    print(f"  Keeping top {top_n_per_gene} associations per gene by p-value...")
    df = (
        df.sort_values("Pvalue")
          .groupby("gene", sort=False)
          .head(top_n_per_gene)
          .reset_index(drop=True)
    )
    print(f"  After top-{top_n_per_gene} cap: {len(df):,} rows, {df['gene'].nunique():,} genes")

    for col in ["BETA_Burden", "SE_Burden", "Pvalue_Burden", "Pvalue_SKAT"]:
        if col in df.columns:
            df[col] = df[col].replace([np.inf, -np.inf], np.nan)

    print("  Classifying phenotypes...")
    categories = df["pheno_description"].apply(classify_phenotype)
    df["phenotype_category"] = [c[0] for c in categories]
    df["organ_system"] = [c[1] for c in categories]

    records = []
    for _, row in df.iterrows():
        beta = row["BETA_Burden"]
        se   = row["SE_Burden"]
        if pd.isna(beta):
            beta = 0.0
        if pd.isna(se):
            se = 0.0
        p_burden = row.get("Pvalue_Burden")
        p_skat   = row.get("Pvalue_SKAT")
        if pd.isna(p_burden):
            p_burden = None
        if pd.isna(p_skat):
            p_skat = None
        records.append({
            "gene_symbol":        row["gene"],
            "phenotype":          row["pheno_description"],
            "phenotype_category": row["phenotype_category"],
            "organ_system":       row["organ_system"],
            "p_value":            float(row["Pvalue"]),
            "p_value_burden":     float(p_burden) if p_burden is not None else None,
            "p_value_skat":       float(p_skat) if p_skat is not None else None,
            "beta":               round(float(beta), 6),
            "se":                 round(float(se), 6),
            "n_carriers":         None,
            "direction":          "loss" if beta < 0 else "gain",
        })
    print(f"  Prepared {len(records):,} pLoF association records "
          f"({df['gene'].nunique():,} genes)")
    return records


def load_real_loeuf_data(gene_universe: dict[str, tuple[str, str]]) -> list[dict]:
    """Load LOEUF gene constraint scores from grr.iossifovlab.com.

    Source: https://grr.iossifovlab.com/gene_properties/gene_scores/LOEUF/
    File:   data/LOEUF_scores.csv.gz

    Columns in source: gene, LOEUF, LOEUF_rank
      LOEUF: 90% upper CI for observed/expected pLoF ratio (0.03–2.0)
             Lower = more constrained (intolerant to LOF)
      LOEUF_rank: rank among 19,200 genes (1 = most constrained)

    pLI is not available from this source; mis_z is set to 0.0 (placeholder).
    risk_class is derived from LOEUF thresholds following gnomAD conventions:
      critical  : LOEUF < 0.35  (equivalent to pLI ~ 0.9+)
      high      : LOEUF < 0.70
      moderate  : LOEUF < 1.00
      low       : LOEUF >= 1.00
    """
    print(f"  Loading {LOEUF_FILE} ...")
    df = pd.read_csv(LOEUF_FILE)
    df = df.dropna(subset=["LOEUF"])   # drop genes without LOEUF (non-coding / sparse)
    # Keep minimum LOEUF per gene (most constrained) for the ~42 duplicated gene symbols
    df = df.groupby("gene", as_index=False)["LOEUF"].min()
    print(f"  LOEUF dataset: {len(df):,} unique genes (after dedup + NaN drop)")

    records = []
    for _, row in df.iterrows():
        gene: str = row["gene"]
        loeuf: float = float(row["LOEUF"])

        # Get Ensembl ID from gene_universe (CellxGene) if available
        ensembl_id = gene_universe.get(gene, ("", ""))[0]

        # Derive risk class from LOEUF (gnomAD-aligned thresholds)
        if loeuf < 0.35:
            risk_class = "critical"
        elif loeuf < 0.70:
            risk_class = "high"
        elif loeuf < 1.00:
            risk_class = "moderate"
        else:
            risk_class = "low"

        records.append({
            "gene_symbol":  gene,
            "ensembl_id":   ensembl_id,
            "pli_score":    0.0,        # not available from this source
            "loeuf_score":  round(loeuf, 4),
            "mis_z_score":  0.0,        # not available from this source
            "risk_class":   risk_class,
        })

    print(f"  Prepared {len(records):,} LOEUF records")
    return records


# ────────────────────────────────────────────────────────────────────────────
# Database creation
# ────────────────────────────────────────────────────────────────────────────

def create_database(db_path: str) -> None:
    conn   = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute("PRAGMA journal_mode=WAL")
    cursor.execute("PRAGMA synchronous=NORMAL")

    # ── Expression table ────────────────────────────────────────
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS expression_summary (
            id              INTEGER PRIMARY KEY AUTOINCREMENT,
            gene_symbol     TEXT    NOT NULL,
            ensembl_id      TEXT    NOT NULL,
            cell_type       TEXT    NOT NULL,
            tissue          TEXT    NOT NULL,
            organ           TEXT    NOT NULL,
            mean_expression REAL    NOT NULL,
            pct_expressed   REAL    NOT NULL,
            n_cells         INTEGER NOT NULL
        )
    """)
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_expr_gene ON expression_summary(gene_symbol)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_expr_gene_cell ON expression_summary(gene_symbol, cell_type)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_expr_tissue ON expression_summary(tissue)")

    # ── pLOF table ───────────────────────────────────────────────
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS plof_associations (
            id                  INTEGER PRIMARY KEY AUTOINCREMENT,
            gene_symbol         TEXT    NOT NULL,
            phenotype           TEXT    NOT NULL,
            phenotype_category  TEXT    NOT NULL,
            organ_system        TEXT    NOT NULL,
            p_value             REAL    NOT NULL,
            p_value_burden      REAL,
            p_value_skat        REAL,
            beta                REAL    NOT NULL,
            se                  REAL    NOT NULL,
            n_carriers          INTEGER,
            direction           TEXT    NOT NULL
        )
    """)
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_plof_gene ON plof_associations(gene_symbol)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_plof_gene_pheno ON plof_associations(gene_symbol, phenotype)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_plof_pvalue ON plof_associations(p_value)")

    # ── Dosage sensitivity table ─────────────────────────────────
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS gene_dosage_sensitivity (
            id           INTEGER PRIMARY KEY AUTOINCREMENT,
            gene_symbol  TEXT    NOT NULL UNIQUE,
            ensembl_id   TEXT    NOT NULL,
            pli_score    REAL    NOT NULL,
            loeuf_score  REAL    NOT NULL,
            mis_z_score  REAL    NOT NULL,
            risk_class   TEXT    NOT NULL
        )
    """)
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_dosage_gene ON gene_dosage_sensitivity(gene_symbol)")

    # ── Prediction jobs table ────────────────────────────────────
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS prediction_jobs (
            id              TEXT    PRIMARY KEY,
            status          TEXT    NOT NULL DEFAULT 'pending',
            sequence_type   TEXT    NOT NULL,
            sequence_name   TEXT    DEFAULT '',
            sequence        TEXT    NOT NULL,
            heavy_chain     TEXT,
            light_chain     TEXT,
            created_at      DATETIME DEFAULT CURRENT_TIMESTAMP,
            completed_at    DATETIME,
            n_targets_found INTEGER DEFAULT 0,
            error_message   TEXT,
            predictor_used  TEXT    DEFAULT 'mock'
        )
    """)
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_job_status ON prediction_jobs(status)")

    # ── Binding predictions table ────────────────────────────────
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS binding_predictions (
            id              INTEGER PRIMARY KEY AUTOINCREMENT,
            job_id          TEXT    NOT NULL,
            gene_symbol     TEXT    NOT NULL,
            ensembl_id      TEXT    NOT NULL,
            binding_score   REAL    NOT NULL,
            confidence      REAL    NOT NULL,
            binding_site    TEXT,
            interaction_type TEXT   NOT NULL,
            delta_g         REAL    NOT NULL,
            kd_nm           REAL    NOT NULL,
            rank            INTEGER NOT NULL
        )
    """)
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_binding_job  ON binding_predictions(job_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_binding_gene ON binding_predictions(gene_symbol)")

    # ── Build gene universe (LOEUF × CellxGene intersection) ────
    print("\nBuilding gene universe...")
    gene_universe = build_gene_universe()

    # ── Insert real expression data ──────────────────────────────
    print("\nLoading REAL CellxGene expression data...")
    expression_data = load_real_expression_data(gene_universe)
    cursor.executemany(
        "INSERT INTO expression_summary "
        "(gene_symbol, ensembl_id, cell_type, tissue, organ, mean_expression, pct_expressed, n_cells) "
        "VALUES (:gene_symbol, :ensembl_id, :cell_type, :tissue, :organ, :mean_expression, :pct_expressed, :n_cells)",
        expression_data,
    )
    print(f"  Inserted {len(expression_data):,} REAL expression records")

    # ── Insert real pLoF data ────────────────────────────────────
    print("\nLoading REAL Genebass pLoF data...")
    plof_data = load_real_plof_data()
    cursor.executemany(
        "INSERT INTO plof_associations "
        "(gene_symbol, phenotype, phenotype_category, organ_system, p_value, p_value_burden, p_value_skat, beta, se, n_carriers, direction) "
        "VALUES (:gene_symbol, :phenotype, :phenotype_category, :organ_system, :p_value, :p_value_burden, :p_value_skat, :beta, :se, :n_carriers, :direction)",
        plof_data,
    )
    print(f"  Inserted {len(plof_data):,} REAL pLoF association records")

    # ── Insert real LOEUF dosage data ────────────────────────────
    print("\nLoading REAL LOEUF dosage sensitivity data...")
    dosage_data = load_real_loeuf_data(gene_universe)
    cursor.executemany(
        "INSERT INTO gene_dosage_sensitivity "
        "(gene_symbol, ensembl_id, pli_score, loeuf_score, mis_z_score, risk_class) "
        "VALUES (:gene_symbol, :ensembl_id, :pli_score, :loeuf_score, :mis_z_score, :risk_class)",
        dosage_data,
    )
    print(f"  Inserted {len(dosage_data):,} REAL LOEUF records")

    conn.commit()

    # ── Summary ──────────────────────────────────────────────────
    cursor.execute("SELECT COUNT(*) FROM plof_associations WHERE p_value < 5e-8")
    n_sig = cursor.fetchone()[0]
    cursor.execute("SELECT COUNT(DISTINCT gene_symbol) FROM gene_dosage_sensitivity")
    n_genes_total = cursor.fetchone()[0]
    cursor.execute("SELECT COUNT(DISTINCT gene_symbol) FROM plof_associations")
    n_genes_plof = cursor.fetchone()[0]
    cursor.execute("SELECT COUNT(DISTINCT phenotype) FROM plof_associations")
    n_phenos = cursor.fetchone()[0]
    cursor.execute("SELECT COUNT(DISTINCT tissue) FROM expression_summary")
    n_tissues = cursor.fetchone()[0]
    cursor.execute("SELECT COUNT(DISTINCT cell_type) FROM expression_summary")
    n_ct = cursor.fetchone()[0]
    cursor.execute("SELECT risk_class, COUNT(*) FROM gene_dosage_sensitivity GROUP BY risk_class")
    rc_counts = dict(cursor.fetchall())

    print(f"\n{'='*60}")
    print(f"Database: {db_path}")
    print(f"  Expression (CellxGene):  {len(expression_data):,} records")
    print(f"                           {n_tissues} tissues, {n_ct} cell types")
    print(f"  pLoF (Genebass):         {len(plof_data):,} records")
    print(f"                           {n_genes_plof} genes, {n_phenos:,} phenotypes")
    print(f"                           {n_sig:,} genome-wide significant (p<5e-8)")
    print(f"  Total genes in DB:       {n_genes_total:,}")
    print(f"  Dosage (LOEUF):          {len(dosage_data):,} records")
    for rc, cnt in sorted(rc_counts.items()):
        print(f"    {rc}: {cnt}")
    print(f"{'='*60}")

    # ── Save JSON sidecar for frontend ───────────────────────────
    demo_dir = DATA_DIR / "demo"
    demo_dir.mkdir(parents=True, exist_ok=True)

    # Use all loaded LOEUF genes (from dosage_data) as the gene list
    all_gene_symbols = sorted({r["gene_symbol"] for r in dosage_data})
    with open(demo_dir / "genes.json", "w") as f:
        json.dump(all_gene_symbols, f, indent=2)

    tissues_dict: dict[str, dict] = {}
    for ct_label, (display_label, tissue, organ) in CELLXGENE_TISSUE_MAP.items():
        if tissue not in tissues_dict:
            tissues_dict[tissue] = {"organ": organ, "cell_types": []}
        if display_label not in tissues_dict[tissue]["cell_types"]:
            tissues_dict[tissue]["cell_types"].append(display_label)
    with open(demo_dir / "tissues.json", "w") as f:
        json.dump(tissues_dict, f, indent=2)

    conn.close()


if __name__ == "__main__":
    db_path = Path(__file__).parent / "humanproof.db"
    if db_path.exists():
        db_path.unlink()
    for wal in [db_path.with_suffix(".db-wal"), db_path.with_suffix(".db-shm")]:
        if wal.exists():
            wal.unlink()
    create_database(str(db_path))
