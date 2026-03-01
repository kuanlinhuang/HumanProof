"""Generate HumanProof database with REAL Genebass pLoF data + synthetic expression.

Loads real pLoF burden-test associations from:
    data/genebass_pLoF_filtered.pkl

Expression data remains synthetic for the 108-gene set.
Dosage sensitivity remains synthetic (real gnomAD scores for Phase 2).
"""

import json
import math
import pickle
import random
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd

random.seed(42)
np.random.seed(42)

DATA_DIR = Path(__file__).parent.parent / "data"
PLOF_PICKLE = DATA_DIR / "genebass_pLoF_filtered.pkl"

# ────────────────────────────────────────────────────────────────────────────
# Gene list: 108 well-known human genes
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
# Phenotype categorizer for UK Biobank / Genebass phenotype descriptions
# ────────────────────────────────────────────────────────────────────────────

# Map of keyword patterns → (category, organ_system)
# Order matters: first match wins, so more specific patterns go first
PHENOTYPE_RULES = [
    # Cancer / Neoplasms
    (["cancer", "carcinoma", "melanoma", "lymphoma", "leukemia", "leukaemia",
      "neoplasm", "tumor", "tumour", "malignant", "myeloma", "sarcoma",
      "mesothelioma", "glioma"], "cancer", "Cancer"),
    # Cardiovascular
    (["heart", "cardiac", "coronary", "myocardial", "atrial", "ventricular",
      "angina", "arrhythmia", "aortic", "hypertension", "blood pressure",
      "pulse rate", "heart rate", "systolic", "diastolic", "vascular",
      "peripheral vascular", "arterial", "stroke", "cerebrovascular",
      "thrombosis", "embolism", "aneurysm", "ischaemic heart",
      "cardiomyopathy", "valve", "ecg", "electrocardiogr"], "cardiovascular", "Cardiovascular"),
    # Hematologic / Blood
    (["haemoglobin", "hemoglobin", "platelet", "erythrocyte", "leukocyte",
      "lymphocyte", "monocyte", "neutrophil", "eosinophil", "basophil",
      "reticulocyte", "haematocrit", "hematocrit", "red blood cell",
      "white blood cell", "mean corpuscular", "mcv", "mch", "mchc",
      "red cell", "blood count", "anaemia", "anemia", "bleeding",
      "coagulation", "clotting", "fibrinogen", "von willebrand",
      "sickle cell", "thalassaemia", "thalassemia", "haemophilia",
      "hemophilia", "iron deficiency"], "hematologic", "Hematologic"),
    # Neurological / Psychiatric
    (["alzheimer", "parkinson", "dementia", "epilepsy", "seizure",
      "multiple sclerosis", "neuropathy", "migraine", "headache",
      "neuralgia", "brain", "cerebral", "meningitis", "encephalitis",
      "depression", "anxiety", "schizophrenia", "bipolar", "psychosis",
      "psychiatric", "insomnia", "sleep", "neurodegenerative", "motor neuron",
      "als", "huntington", "tremor", "cognitive", "intelligence",
      "reaction time", "prospective memory", "fluid intelligence",
      "nerve", "mental health"], "neurological", "Neurological"),
    # Metabolic / Endocrine
    (["diabetes", "glucose", "insulin", "hba1c", "glycated", "cholesterol",
      "triglyceride", "lipid", "hdl", "ldl", "apolipoprotein", "bmi",
      "body mass index", "obesity", "waist", "hip circumference", "weight",
      "body fat", "impedance", "basal metabolic", "metabolic", "gout",
      "urate", "uric acid", "adiposity", "fat mass", "lean mass",
      "thyroid", "tsh", "thyroxine", "testosterone", "oestradiol",
      "estradiol", "shbg", "igf-1", "growth hormone", "cortisol",
      "vitamin d", "calcium", "phosphate", "parathyroid",
      "glycoprotein acetyls", "fatty acid"], "metabolic", "Metabolic"),
    # Hepatic / Liver
    (["liver", "hepat", "bilirubin", "albumin", "alt", "ast", "ggt",
      "gamma-glutamyl", "alkaline phosphatase", "cirrhosis", "fatty liver",
      "gallstone", "gallbladder", "cholecyst", "jaundice", "bile",
      "cholelithiasis"], "hepatic", "Hepatic"),
    # Renal / Kidney
    (["kidney", "renal", "creatinine", "cystatin", "glomerular",
      "nephr", "urea", "urinary", "urine", "bladder", "dialysis",
      "microalbumin", "proteinuria", "kidney stone"], "renal", "Renal"),
    # Respiratory / Pulmonary
    (["lung", "pulmonary", "respiratory", "asthma", "copd", "pneumonia",
      "bronch", "fev1", "fvc", "peak expiratory", "spirometry",
      "emphysema", "fibrosis", "pleural", "sleep apn", "oxygen sat",
      "breathless"], "respiratory", "Respiratory"),
    # Musculoskeletal
    (["bone", "fracture", "osteopor", "osteoarth", "arthritis", "joint",
      "muscle", "muscular", "dystrophy", "back pain", "spine", "spinal",
      "scoliosis", "grip strength", "hand grip", "sarcopenia",
      "rheumatoid", "gout", "height", "standing height", "sitting height",
      "heel bone", "bone density", "mineral density"], "musculoskeletal", "Musculoskeletal"),
    # Immunologic / Autoimmune
    (["immune", "autoimmune", "lupus", "psoriasis", "crohn", "colitis",
      "inflammatory bowel", "celiac", "coeliac", "sjogren", "vasculitis",
      "sarcoidosis", "allergy", "allergic", "atopic", "eczema", "dermatitis",
      "urticaria", "hay fever", "c-reactive protein", "crp",
      "immunoglobulin", "igg", "ige", "complement"], "immunologic", "Immunologic"),
    # Ophthalmologic / Eye
    (["eye", "ocular", "retinal", "glaucoma", "cataract", "macular",
      "myopia", "hypermetropia", "astigmatism", "intraocular pressure",
      "corneal", "visual acuity", "refractive error", "optic",
      "blindness", "spectacle", "glasses"], "ophthalmologic", "Ophthalmologic"),
    # Dermatologic / Skin
    (["skin", "dermat", "acne", "psoriasis", "eczema", "melanoma",
      "hair", "alopecia", "nail", "pigment", "sunburn",
      "skin colour", "skin color"], "dermatologic", "Dermatologic"),
    # Reproductive / Urogenital
    (["prostate", "ovarian", "uterine", "cervical", "breast",
      "endometri", "menstrua", "menopause", "fertility", "pregnancy",
      "birth weight", "gestational", "miscarriage", "contraceptive",
      "erectile", "impotence", "testicular"], "reproductive", "Reproductive"),
    # Gastrointestinal
    (["gastric", "stomach", "intestin", "bowel", "colon", "rectal",
      "colorectal", "oesophag", "esophag", "reflux", "hernia",
      "appendic", "diverticul", "irritable bowel", "ibs",
      "constipation", "diarrhoea", "diarrhea", "peptic", "ulcer",
      "pancreat", "abdominal"], "gastrointestinal", "Gastrointestinal"),
    # Audiologic / ENT
    (["hearing", "tinnitus", "ear", "deafness", "otitis",
      "noise exposure", "cochlear"], "audiologic", "Audiologic"),
    # Anthropometric (general body measurements not fitting above)
    (["arm span", "trunk", "leg", "body size", "birth weight",
      "comparative body size", "handedness"], "anthropometric", "Anthropometric"),
]

# Build a compiled lookup for speed
_CATEGORY_CACHE: dict[str, tuple[str, str]] = {}


def classify_phenotype(description: str) -> tuple[str, str]:
    """Classify a phenotype description into (category, organ_system)."""
    if description in _CATEGORY_CACHE:
        return _CATEGORY_CACHE[description]

    desc_lower = description.lower()
    for keywords, category, organ in PHENOTYPE_RULES:
        for kw in keywords:
            if kw in desc_lower:
                _CATEGORY_CACHE[description] = (category, organ)
                return category, organ

    # Fallback: if it contains common medication patterns
    if any(kw in desc_lower for kw in ["medication", "treatment", "prescribed", "tablet"]):
        _CATEGORY_CACHE[description] = ("medication", "Medication")
        return "medication", "Medication"

    # Fallback: procedural
    if any(kw in desc_lower for kw in ["operation", "surgery", "procedure", "operative"]):
        _CATEGORY_CACHE[description] = ("procedural", "Procedural")
        return "procedural", "Procedural"

    _CATEGORY_CACHE[description] = ("other", "Other")
    return "other", "Other"


# ────────────────────────────────────────────────────────────────────────────
# Tissues, organs, and cell types (synthetic expression)
# ────────────────────────────────────────────────────────────────────────────
TISSUE_CELL_TYPES = {
    "Brain": {
        "organ": "Brain",
        "cell_types": ["Neurons", "Astrocytes", "Microglia", "Oligodendrocytes", "Endothelial cells"],
    },
    "Heart": {
        "organ": "Heart",
        "cell_types": ["Cardiomyocytes", "Cardiac fibroblasts", "Endothelial cells", "Smooth muscle cells"],
    },
    "Liver": {
        "organ": "Liver",
        "cell_types": ["Hepatocytes", "Kupffer cells", "Stellate cells", "Cholangiocytes", "Liver sinusoidal endothelial cells"],
    },
    "Lung": {
        "organ": "Lung",
        "cell_types": ["Alveolar type I", "Alveolar type II", "Alveolar macrophages", "Ciliated epithelial", "Club cells"],
    },
    "Kidney": {
        "organ": "Kidney",
        "cell_types": ["Proximal tubular", "Distal tubular", "Podocytes", "Collecting duct", "Mesangial cells"],
    },
    "Intestine": {
        "organ": "Intestine",
        "cell_types": ["Enterocytes", "Goblet cells", "Paneth cells", "Enteroendocrine", "Intestinal stem cells"],
    },
    "Bone marrow": {
        "organ": "Bone marrow",
        "cell_types": ["HSCs", "Erythroblasts", "Myeloid progenitors", "Megakaryocytes", "Plasma cells"],
    },
    "Lymph node": {
        "organ": "Lymph node",
        "cell_types": ["B cells", "T cells", "Dendritic cells", "Macrophages", "NK cells"],
    },
    "Spleen": {
        "organ": "Spleen",
        "cell_types": ["Red pulp macrophages", "B cells", "T cells", "Dendritic cells"],
    },
    "Skin": {
        "organ": "Skin",
        "cell_types": ["Keratinocytes", "Melanocytes", "Fibroblasts", "Langerhans cells"],
    },
    "Pancreas": {
        "organ": "Pancreas",
        "cell_types": ["Beta cells", "Alpha cells", "Delta cells", "Acinar cells", "Ductal cells"],
    },
    "Muscle": {
        "organ": "Muscle",
        "cell_types": ["Skeletal myocytes", "Satellite cells", "Fibroblasts"],
    },
    "Thyroid": {
        "organ": "Thyroid",
        "cell_types": ["Follicular cells", "Parafollicular C cells"],
    },
    "Adrenal": {
        "organ": "Adrenal",
        "cell_types": ["Cortical zona glomerulosa", "Cortical zona fasciculata", "Medullary chromaffin"],
    },
    "Breast": {
        "organ": "Breast",
        "cell_types": ["Luminal epithelial", "Myoepithelial", "Adipocytes", "Fibroblasts"],
    },
}

GENE_TISSUE_BIAS = {
    "ALB": {"Liver": 5.0}, "INS": {"Pancreas": 5.0},
    "GCG": {"Pancreas": 4.0, "Intestine": 2.0}, "EPO": {"Kidney": 4.0},
    "HBB": {"Bone marrow": 5.0}, "HBA1": {"Bone marrow": 5.0},
    "CYP3A4": {"Liver": 5.0, "Intestine": 2.0}, "CYP2D6": {"Liver": 4.0},
    "CFTR": {"Lung": 3.0, "Intestine": 3.0, "Pancreas": 2.0},
    "DMD": {"Muscle": 5.0, "Heart": 3.0}, "ERBB2": {"Breast": 3.0},
    "EGFR": {"Lung": 3.0, "Skin": 2.5, "Intestine": 2.0},
    "PDCD1": {"Lymph node": 4.0, "Spleen": 3.0},
    "CD274": {"Lymph node": 3.0, "Spleen": 2.0},
    "FOXP3": {"Lymph node": 4.0},
    "CD4": {"Lymph node": 4.0, "Spleen": 3.0},
    "CD8A": {"Lymph node": 4.0, "Spleen": 3.0},
    "APP": {"Brain": 4.0}, "MAPT": {"Brain": 5.0}, "SNCA": {"Brain": 4.0},
    "HTT": {"Brain": 3.0}, "SOD1": {"Brain": 2.0, "Muscle": 2.0},
    "TNF": {"Lymph node": 3.0, "Bone marrow": 2.0},
    "IL6": {"Bone marrow": 2.0, "Lymph node": 2.0},
    "VEGFA": {"Kidney": 2.0, "Lung": 2.0},
    "APOE": {"Liver": 3.0, "Brain": 3.0},
    "NOTCH1": {"Bone marrow": 3.0, "Intestine": 2.0},
    "F8": {"Liver": 4.0}, "FGA": {"Liver": 5.0},
    "THPO": {"Liver": 4.0}, "CSF3": {"Bone marrow": 3.0},
}


# ────────────────────────────────────────────────────────────────────────────
# Data generators
# ────────────────────────────────────────────────────────────────────────────

def generate_expression_data():
    """Generate synthetic expression data for all genes across tissues/cell types."""
    records = []
    for gene, ensembl_id in GENES.items():
        biases = GENE_TISSUE_BIAS.get(gene, {})
        for tissue, info in TISSUE_CELL_TYPES.items():
            tissue_bias = biases.get(tissue, 1.0)
            base_expr = random.uniform(0.1, 2.0) * tissue_bias
            for cell_type in info["cell_types"]:
                ct_multiplier = random.uniform(0.2, 2.0)
                mean_expr = max(0.0, base_expr * ct_multiplier + random.gauss(0, 0.3))
                pct_base = min(1.0, mean_expr / 5.0)
                pct_expressed = max(0.0, min(1.0, pct_base + random.gauss(0, 0.1)))
                n_cells = random.randint(50, 5000)
                records.append({
                    "gene_symbol": gene, "ensembl_id": ensembl_id,
                    "cell_type": cell_type, "tissue": tissue, "organ": info["organ"],
                    "mean_expression": round(mean_expr, 4),
                    "pct_expressed": round(pct_expressed, 4), "n_cells": n_cells,
                })
    return records


def load_real_plof_data():
    """Load real Genebass pLoF data for our 108 genes.

    Returns list of dicts ready for SQLite insertion.
    """
    print(f"  Loading {PLOF_PICKLE} ...")
    df = pd.read_pickle(PLOF_PICKLE)
    print(f"  Full dataset: {len(df):,} rows, {df['gene'].nunique()} genes")

    # Filter to our 108 genes
    gene_set = set(GENES.keys())
    df = df[df["gene"].isin(gene_set)].copy()
    print(f"  After filtering to {len(gene_set)} target genes: {len(df):,} rows")

    # Drop rows where combined Pvalue is NaN
    df = df.dropna(subset=["Pvalue"])

    # Clean infinite and NaN values in numeric columns
    for col in ["BETA_Burden", "SE_Burden", "Pvalue_Burden", "Pvalue_SKAT"]:
        if col in df.columns:
            df[col] = df[col].replace([np.inf, -np.inf], np.nan)

    # Classify phenotypes
    print("  Classifying phenotypes...")
    categories = df["pheno_description"].apply(classify_phenotype)
    df["phenotype_category"] = [c[0] for c in categories]
    df["organ_system"] = [c[1] for c in categories]

    # Report category distribution
    cat_counts = df["phenotype_category"].value_counts()
    print(f"  Phenotype categories ({len(cat_counts)}):")
    for cat, count in cat_counts.items():
        print(f"    {cat}: {count:,}")

    # Build records
    records = []
    for _, row in df.iterrows():
        beta = row["BETA_Burden"]
        se = row["SE_Burden"]

        # Handle NaN beta/se: default to 0
        if pd.isna(beta):
            beta = 0.0
        if pd.isna(se):
            se = 0.0

        p_burden = row.get("Pvalue_Burden")
        p_skat = row.get("Pvalue_SKAT")
        if pd.isna(p_burden):
            p_burden = None
        if pd.isna(p_skat):
            p_skat = None

        records.append({
            "gene_symbol": row["gene"],
            "phenotype": row["pheno_description"],
            "phenotype_category": row["phenotype_category"],
            "organ_system": row["organ_system"],
            "p_value": float(row["Pvalue"]),
            "p_value_burden": float(p_burden) if p_burden is not None else None,
            "p_value_skat": float(p_skat) if p_skat is not None else None,
            "beta": round(float(beta), 6),
            "se": round(float(se), 6),
            "n_carriers": None,  # Not available in Genebass data
            "direction": "loss" if beta < 0 else "gain",
        })

    print(f"  Prepared {len(records):,} pLoF association records")
    return records


def generate_dosage_data():
    """Generate synthetic gene dosage sensitivity scores."""
    records = []
    for gene, ensembl_id in GENES.items():
        pli = random.betavariate(2, 5)
        loeuf = random.uniform(0.05, 2.0)
        mis_z = random.gauss(0, 2)
        if pli > 0.9 and loeuf < 0.2:
            risk_class = "critical"
        elif pli > 0.5 or loeuf < 0.5:
            risk_class = "high"
        elif pli > 0.2 or loeuf < 1.0:
            risk_class = "moderate"
        else:
            risk_class = "low"
        if gene in ("TP53", "RB1", "BRCA1", "BRCA2", "PTEN", "APC"):
            pli = random.uniform(0.95, 1.0)
            loeuf = random.uniform(0.05, 0.15)
            risk_class = "critical"
        elif gene in ("GAPDH", "ACTB"):
            pli = random.uniform(0.8, 0.95)
            loeuf = random.uniform(0.1, 0.3)
            risk_class = "high"
        records.append({
            "gene_symbol": gene, "ensembl_id": ensembl_id,
            "pli_score": round(pli, 4), "loeuf_score": round(loeuf, 4),
            "mis_z_score": round(mis_z, 4), "risk_class": risk_class,
        })
    return records


# ────────────────────────────────────────────────────────────────────────────
# Database creation
# ────────────────────────────────────────────────────────────────────────────

def create_database(db_path: str):
    """Create SQLite database and populate with real pLoF + synthetic expression."""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Pragmas for faster bulk insert
    cursor.execute("PRAGMA journal_mode=WAL")
    cursor.execute("PRAGMA synchronous=NORMAL")

    # ── Expression table ────────────────────────────────────────
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS expression_summary (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            gene_symbol TEXT NOT NULL,
            ensembl_id TEXT NOT NULL,
            cell_type TEXT NOT NULL,
            tissue TEXT NOT NULL,
            organ TEXT NOT NULL,
            mean_expression REAL NOT NULL,
            pct_expressed REAL NOT NULL,
            n_cells INTEGER NOT NULL
        )
    """)
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_expr_gene ON expression_summary(gene_symbol)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_expr_gene_cell ON expression_summary(gene_symbol, cell_type)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_expr_tissue ON expression_summary(tissue)")

    # ── pLOF table (updated schema) ─────────────────────────────
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS plof_associations (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            gene_symbol TEXT NOT NULL,
            phenotype TEXT NOT NULL,
            phenotype_category TEXT NOT NULL,
            organ_system TEXT NOT NULL,
            p_value REAL NOT NULL,
            p_value_burden REAL,
            p_value_skat REAL,
            beta REAL NOT NULL,
            se REAL NOT NULL,
            n_carriers INTEGER,
            direction TEXT NOT NULL
        )
    """)
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_plof_gene ON plof_associations(gene_symbol)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_plof_gene_pheno ON plof_associations(gene_symbol, phenotype)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_plof_pvalue ON plof_associations(p_value)")

    # ── Dosage sensitivity table ────────────────────────────────
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS gene_dosage_sensitivity (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            gene_symbol TEXT NOT NULL UNIQUE,
            ensembl_id TEXT NOT NULL,
            pli_score REAL NOT NULL,
            loeuf_score REAL NOT NULL,
            mis_z_score REAL NOT NULL,
            risk_class TEXT NOT NULL
        )
    """)
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_dosage_gene ON gene_dosage_sensitivity(gene_symbol)")

    # ── Prediction jobs table ───────────────────────────────────
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS prediction_jobs (
            id TEXT PRIMARY KEY,
            status TEXT NOT NULL DEFAULT 'pending',
            sequence_type TEXT NOT NULL,
            sequence_name TEXT DEFAULT '',
            sequence TEXT NOT NULL,
            heavy_chain TEXT,
            light_chain TEXT,
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
            completed_at DATETIME,
            n_targets_found INTEGER DEFAULT 0,
            error_message TEXT,
            predictor_used TEXT DEFAULT 'mock'
        )
    """)
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_job_status ON prediction_jobs(status)")

    # ── Binding predictions table ───────────────────────────────
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS binding_predictions (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            job_id TEXT NOT NULL,
            gene_symbol TEXT NOT NULL,
            ensembl_id TEXT NOT NULL,
            binding_score REAL NOT NULL,
            confidence REAL NOT NULL,
            binding_site TEXT,
            interaction_type TEXT NOT NULL,
            delta_g REAL NOT NULL,
            kd_nm REAL NOT NULL,
            rank INTEGER NOT NULL
        )
    """)
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_binding_job ON binding_predictions(job_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_binding_gene ON binding_predictions(gene_symbol)")

    # ── Insert expression data ──────────────────────────────────
    print("Generating synthetic expression data...")
    expression_data = generate_expression_data()
    cursor.executemany(
        "INSERT INTO expression_summary (gene_symbol, ensembl_id, cell_type, tissue, organ, mean_expression, pct_expressed, n_cells) VALUES (:gene_symbol, :ensembl_id, :cell_type, :tissue, :organ, :mean_expression, :pct_expressed, :n_cells)",
        expression_data,
    )
    print(f"  Inserted {len(expression_data):,} expression records")

    # ── Insert REAL pLOF data ───────────────────────────────────
    print("\nLoading REAL Genebass pLoF data...")
    plof_data = load_real_plof_data()
    cursor.executemany(
        "INSERT INTO plof_associations (gene_symbol, phenotype, phenotype_category, organ_system, p_value, p_value_burden, p_value_skat, beta, se, n_carriers, direction) VALUES (:gene_symbol, :phenotype, :phenotype_category, :organ_system, :p_value, :p_value_burden, :p_value_skat, :beta, :se, :n_carriers, :direction)",
        plof_data,
    )
    print(f"  Inserted {len(plof_data):,} REAL pLoF association records")

    # ── Insert dosage data ──────────────────────────────────────
    print("\nGenerating synthetic dosage sensitivity data...")
    dosage_data = generate_dosage_data()
    cursor.executemany(
        "INSERT INTO gene_dosage_sensitivity (gene_symbol, ensembl_id, pli_score, loeuf_score, mis_z_score, risk_class) VALUES (:gene_symbol, :ensembl_id, :pli_score, :loeuf_score, :mis_z_score, :risk_class)",
        dosage_data,
    )
    print(f"  Inserted {len(dosage_data):,} dosage records")

    conn.commit()

    # ── Summary stats ───────────────────────────────────────────
    cursor.execute("SELECT COUNT(*) FROM plof_associations WHERE p_value < 5e-8")
    n_sig = cursor.fetchone()[0]
    cursor.execute("SELECT COUNT(DISTINCT gene_symbol) FROM plof_associations")
    n_genes_plof = cursor.fetchone()[0]
    cursor.execute("SELECT COUNT(DISTINCT phenotype) FROM plof_associations")
    n_phenos = cursor.fetchone()[0]

    print(f"\n{'='*60}")
    print(f"Database created at: {db_path}")
    print(f"  Expression: {len(expression_data):,} records (synthetic)")
    print(f"  pLoF:       {len(plof_data):,} records (REAL Genebass)")
    print(f"              {n_genes_plof} genes, {n_phenos:,} phenotypes")
    print(f"              {n_sig:,} genome-wide significant (p < 5e-8)")
    print(f"  Dosage:     {len(dosage_data):,} records (synthetic)")
    print(f"{'='*60}")

    # Save JSON for frontend
    demo_dir = Path(__file__).parent.parent / "data" / "demo"
    demo_dir.mkdir(parents=True, exist_ok=True)
    with open(demo_dir / "genes.json", "w") as f:
        json.dump(list(GENES.keys()), f, indent=2)
    with open(demo_dir / "tissues.json", "w") as f:
        tissues_list = {
            tissue: {"organ": info["organ"], "cell_types": info["cell_types"]}
            for tissue, info in TISSUE_CELL_TYPES.items()
        }
        json.dump(tissues_list, f, indent=2)

    conn.close()


if __name__ == "__main__":
    db_path = Path(__file__).parent / "humanproof.db"
    if db_path.exists():
        db_path.unlink()
    create_database(str(db_path))
