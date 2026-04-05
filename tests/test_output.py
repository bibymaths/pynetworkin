import pandas as pd
import tempfile
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from output import write_tsv, write_cytoscape, STANDARD_COLUMNS


SAMPLE_PREDICTIONS = [
    {
        "Name": "SRC",
        "Position": "S17",
        "Tree": "KIN",
        "Motif Group": "PRKACA",
        "Kinase/Phosphatase/Phospho-binding domain": "PRKACA",
        "NetworKIN score": 0.78,
        "Motif probability": 0.82,
        "STRING score": 0.74,
        "Target STRING ID": "9606.ENSP00000418787",
        "Kinase STRING ID": "9606.ENSP00000356969",
        "Target Name": "SRC",
        "Kinase Name": "PRKACA",
        "Target description": "Proto-oncogene tyrosine-protein kinase Src",
        "Kinase description": "cAMP-dependent protein kinase catalytic subunit alpha",
        "Peptide sequence window": "RLIEDNEytARQGAK",
        "Intermediate nodes": "notdef",
        "recovered": False,
        "recovery_method": "",
    },
    {
        "Name": "RB1",
        "Position": "S807",
        "Tree": "KIN",
        "Motif Group": "CDK2",
        "Kinase/Phosphatase/Phospho-binding domain": "CDK2",
        "NetworKIN score": 0.65,
        "Motif probability": -1.0,
        "STRING score": 0.65,
        "Target STRING ID": "9606.ENSP00000267163",
        "Kinase STRING ID": "9606.ENSP00000266970",
        "Target Name": "RB1",
        "Kinase Name": "CDK2",
        "Target description": "Retinoblastoma-associated protein",
        "Kinase description": "Cyclin-dependent kinase 2",
        "Peptide sequence window": "ISPPPKKsPPKKKLP",
        "Intermediate nodes": "notdef",
        "recovered": True,
        "recovery_method": "context_proximity",
    },
]


def test_tsv_has_all_standard_columns():
    with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as f:
        path = f.name
    write_tsv(SAMPLE_PREDICTIONS, path)
    df = pd.read_csv(path, sep="\t")
    for col in STANDARD_COLUMNS:
        assert col in df.columns, f"Missing column: {col}"
    os.unlink(path)


def test_tsv_column_order_matches_standard():
    with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as f:
        path = f.name
    write_tsv(SAMPLE_PREDICTIONS, path)
    df = pd.read_csv(path, sep="\t")
    assert list(df.columns[:len(STANDARD_COLUMNS)]) == STANDARD_COLUMNS
    os.unlink(path)


def test_recovered_column_is_boolean():
    with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as f:
        path = f.name
    write_tsv(SAMPLE_PREDICTIONS, path)
    df = pd.read_csv(path, sep="\t")
    assert set(df["recovered"].unique()).issubset({True, False, "True", "False"})
    os.unlink(path)


def test_networkin_score_in_range():
    with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as f:
        path = f.name
    write_tsv(SAMPLE_PREDICTIONS, path)
    df = pd.read_csv(path, sep="\t")
    scored = df[~df["recovered"].astype(str).isin(["True"])]
    assert (scored["NetworKIN score"].between(0.0, 1.0)).all(), \
        "Non-recovered NetworKIN score must be in [0, 1]"
    os.unlink(path)


def test_no_null_kinase_or_substrate():
    with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as f:
        path = f.name
    write_tsv(SAMPLE_PREDICTIONS, path)
    df = pd.read_csv(path, sep="\t")
    assert df["Kinase/Phosphatase/Phospho-binding domain"].notna().all()
    assert df["Name"].notna().all()
    os.unlink(path)


def test_cytoscape_has_required_columns():
    with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as f:
        path = f.name
    write_cytoscape(SAMPLE_PREDICTIONS, path)
    df = pd.read_csv(path, sep="\t")
    assert set(df.columns) == {"source", "interaction", "target", "weight"}
    os.unlink(path)


def test_cytoscape_no_null_edges():
    with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as f:
        path = f.name
    write_cytoscape(SAMPLE_PREDICTIONS, path)
    df = pd.read_csv(path, sep="\t")
    assert df["source"].notna().all()
    assert df["target"].notna().all()
    os.unlink(path)
