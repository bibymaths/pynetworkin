from __future__ import annotations

import glob
import gzip
import os
import platform
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from pynetworkin.inputs.phosphosites import fetch_phosphosite
from pynetworkin.inputs.string_network import fetch_string_network
from pynetworkin.likelihood import ConvertScore2L, ReadConversionTableBin, ReadConversionTableFromMemory
from pynetworkin.logger import logger
from pynetworkin.graph_scoring import filter_and_rank_predictions
from pynetworkin.motif_scoring import score_sequences
from pynetworkin.output import write_output
from pynetworkin.recovery import recover_false_negatives
from pynetworkin.resources import conversion_parquet_path

ALPHAS = {"9606": 0.85, "4932": 0.65}
SPECIES_NAME = {"9606": "human", "4932": "yeast"}
PENALTY = {
    "9606": {"hub penalty": 100, "length penalty": 800},
    "4932": {"hub penalty": 170, "length penalty": 1000},
}

NETWORKIN_SITE_FILE = 1
PROTEOME_DISCOVERER_SITE_FILE = 2
MAX_QUANT_DIRECT_OUTPUT_FILE = 3
LEGACY_SITE_FILE = 4
MS_MCMC_FILE = 5


@dataclass(slots=True)
class AppConfig:
    organism: str
    fasta_path: str
    sites_path: str | None
    datadir: str
    blast_dir: str
    threads: int = 1
    active_threads: int = 2
    mode: str | bool = False
    path_mode: str = "direct"
    verbose: bool = False
    fast: bool = False
    leave_intermediates: bool = False
    string_for_uncovered: bool = False
    refresh: bool = False
    result_dir: str = "results"
    temp_dir: str = "tmp"

    @property
    def species_name(self) -> str:
        return SPECIES_NAME[self.organism]

    @property
    def fasta_stem(self) -> str:
        return Path(self.fasta_path).name

    @property
    def blast_output_path(self) -> str:
        return f"{self.fasta_path}.{self.organism}.blast.out"


@dataclass(slots=True)
class RunArtifacts:
    id_seq: dict[str, str]
    id_pos_res: dict[str, dict[int, str]]
    incoming2string: dict[str, dict[str, bool]]
    string2incoming: dict[str, dict[str, bool]]
    map_group_to_domain: dict[str, dict[str, list[str]]]
    string_alias: dict[str, str]
    string_desc: dict[str, str]
    name_hash: dict[str, list[str]]
    tree_pred_string_data: dict[str, dict[str, dict[str, Any]]]
    id_pos_tree_pred: dict[str, Any]


@dataclass(slots=True)
class PredictionStats:
    proteins_seen: int = 0
    proteins_mapped: int = 0
    proteins_unmapped: int = 0
    string_targets_seen: int = 0
    string_scores_seen: int = 0
    uncovered_branch_hits: int = 0
    covered_branch_hits: int = 0
    result_rows: int = 0
    unmatched_names: int = 0

    def as_dict(self) -> dict[str, int]:
        return {
            "proteins_seen": self.proteins_seen,
            "proteins_mapped": self.proteins_mapped,
            "proteins_unmapped": self.proteins_unmapped,
            "string_targets_seen": self.string_targets_seen,
            "string_scores_seen": self.string_scores_seen,
            "uncovered_branch_hits": self.uncovered_branch_hits,
            "covered_branch_hits": self.covered_branch_hits,
            "result_rows": self.result_rows,
            "unmatched_names": self.unmatched_names,
        }


class NetworkinError(RuntimeError):
    pass


class CommandError(NetworkinError):
    pass


def ensure_dirs(*paths: str) -> None:
    for path in paths:
        Path(path).mkdir(parents=True, exist_ok=True)
        logger.info("Ensured directory exists: {}", path)


def run_command(cmd: str) -> str:
    # Use bash explicitly with pipefail so downstream pipe tools (like sort) don't swallow blastp's failure exit code
    cmd_with_pipefail = f"set -o pipefail; {cmd}"
    if platform.system() == "Windows":
        command = f'wsl bash -c "{cmd_with_pipefail}"'
    else:
        command = f"bash -c '{cmd_with_pipefail}'"

    logger.muted("Executing command: {}", command)

    try:
        completed = subprocess.run(
            command,
            shell=True,
            check=True,
            capture_output=True,
            text=True,
        )
        return completed.stdout
    except subprocess.CalledProcessError as exc:
        stderr = (exc.stderr or "").strip()
        logger.error("Command failed with exit code {}: {}", exc.returncode, exc.cmd)
        if stderr:
            logger.error("stderr: {}", stderr)
        raise CommandError(stderr or str(exc)) from exc


def read_fasta_file(path: str) -> dict[str, str]:
    aminoacids = re.compile(r"[^ACDEFGHIKLMNPQRSTVWYXB]")
    id_seq: dict[str, str] = {}
    current_id: str | None = None
    chunks: list[str] = []

    def _parse_fasta_id(header: str) -> str:
        header = header.strip()
        token = header.split(None, 1)[0]
        if token.startswith(("sp|", "tr|")):
            parts = token.split("|")
            if len(parts) >= 2 and parts[1]:
                return parts[1]
        return token

    with open(path, encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith(";"):
                continue
            if line.startswith(">"):
                if current_id is not None:
                    id_seq[current_id] = "".join(chunks)
                current_id = _parse_fasta_id(line[1:])
                chunks = []
                continue
            chunks.append(aminoacids.sub("", line))

    if current_id is not None:
        id_seq[current_id] = "".join(chunks)

    logger.success("Loaded {} FASTA sequences", len(id_seq))
    return id_seq


def detect_site_file_type(path: str) -> int:
    with open(path, encoding="utf-8") as handle:
        line = handle.readline().strip()

    if not line:
        raise NetworkinError(f"Empty sites file: {path}")

    if "\t" in line:
        tokens = line.split("\t")
    elif "," in line:
        tokens = line.split(",")
    else:
        tokens = line.split()

    if len(tokens) == 3:
        return NETWORKIN_SITE_FILE
    if len(tokens) == 2:
        return PROTEOME_DISCOVERER_SITE_FILE
    if len(tokens) > 4 and tokens[0] == "Proteins" and (
        tokens[4] == "Leading" or (len(tokens) > 2 and tokens[2].startswith("Leading"))
    ):
        return MAX_QUANT_DIRECT_OUTPUT_FILE
    if len(tokens) > 1 and tokens[1] == "phospho":
        return LEGACY_SITE_FILE
    if Path(path).name.startswith("MS"):
        return MS_MCMC_FILE

    raise NetworkinError(f"Unknown site file format: {path}")


def read_networkin_sites(path: str) -> dict[str, dict[int, str]]:
    id_pos_res: dict[str, dict[int, str]] = {}
    with open(path, encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            tokens = line.rstrip("\n").split("\t")
            if len(tokens) < 2:
                logger.warning("Skipping malformed line {} in {}", line_no, path)
                continue
            protein_id = tokens[0]
            try:
                pos = int(tokens[1])
            except ValueError as exc:
                raise NetworkinError(
                    f"Invalid position at line {line_no}: {line.rstrip()}"
                ) from exc
            residue = tokens[2].strip() if len(tokens) > 2 else ""
            id_pos_res.setdefault(protein_id, {})[pos] = residue

    logger.success("Loaded phosphosite assignments for {} proteins", len(id_pos_res))
    return id_pos_res


def read_proteome_discoverer_sites(path: str) -> dict[str, dict[int, str]]:
    """Parse a Proteome Discoverer 2-column phosphosite file.

    Each line contains a protein ID and either:

    * a residue+position string such as ``"S15"`` or ``"T103"``,
    * a plain integer position (residue defaults to ``""``), or
    * a peptide sequence where lowercase letters mark phosphorylated
      S/T/Y residues (position is 1-based within the peptide string).

    Returns
    -------
    dict[str, dict[int, str]]
        Mapping of protein ID → {1-based position: residue character}.
    """
    _pos_with_residue = re.compile(r"^([A-Za-z])(\d+)$")
    id_pos_res: dict[str, dict[int, str]] = {}
    with open(path, encoding="utf-8") as handle:
        for line_no, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                logger.warning("Skipping malformed line {} in {}", line_no, path)
                continue
            protein_id = parts[0].strip()
            raw_pos = parts[1].strip()
            if not protein_id or not raw_pos:
                continue

            m = _pos_with_residue.match(raw_pos)
            if m:
                # "S15" or "T103" style
                residue = m.group(1).upper()
                pos = int(m.group(2))
                id_pos_res.setdefault(protein_id, {})[pos] = residue
            elif raw_pos.isdigit():
                # Plain integer position
                id_pos_res.setdefault(protein_id, {})[int(raw_pos)] = ""
            else:
                # Peptide sequence – lowercase S/T/Y indicate phosphosites
                for i, ch in enumerate(raw_pos):
                    if ch in "sty":
                        id_pos_res.setdefault(protein_id, {})[i + 1] = ch.upper()

    logger.success("Loaded phosphosite assignments for {} proteins", len(id_pos_res))
    return id_pos_res


def read_max_quant_sites(path: str) -> dict[str, dict[int, str]]:
    """Parse a MaxQuant phosphosite direct-output file.

    The first line is a tab-separated header.  Required columns:

    * ``"Proteins"``             – semicolon-separated protein IDs
    * ``"Positions within proteins"`` – corresponding semicolon-separated positions

    Optional column (residue defaults to ``""`` when absent):

    * ``"Amino acid"`` – single-character residue (same for all proteins in the row)

    Returns
    -------
    dict[str, dict[int, str]]
        Mapping of protein ID → {1-based position: residue character}.
    """
    id_pos_res: dict[str, dict[int, str]] = {}
    with open(path, encoding="utf-8") as handle:
        header_line = handle.readline().rstrip("\n")
        headers = header_line.split("\t")

        try:
            prot_idx = headers.index("Proteins")
            pos_idx = headers.index("Positions within proteins")
        except ValueError as exc:
            raise NetworkinError(
                f"MaxQuant file is missing a required column: {exc}"
            ) from exc

        try:
            res_idx: int | None = headers.index("Amino acid")
        except ValueError:
            res_idx = None

        for line_no, raw_line in enumerate(handle, start=2):
            raw_line = raw_line.rstrip("\n")
            if not raw_line:
                continue
            tokens = raw_line.split("\t")
            if len(tokens) <= max(prot_idx, pos_idx):
                logger.warning("Skipping malformed line {} in {}", line_no, path)
                continue

            proteins = tokens[prot_idx].split(";")
            positions = tokens[pos_idx].split(";")
            residue = (
                tokens[res_idx].strip().upper()
                if res_idx is not None and res_idx < len(tokens)
                else ""
            )

            for prot, pos_str in zip(proteins, positions):
                prot = prot.strip()
                pos_str = pos_str.strip()
                if not prot or not pos_str:
                    continue
                try:
                    pos = int(pos_str)
                except ValueError:
                    logger.warning(
                        "Invalid position '{}' at line {} in {}", pos_str, line_no, path
                    )
                    continue
                id_pos_res.setdefault(prot, {})[pos] = residue

    logger.success("Loaded phosphosite assignments for {} proteins", len(id_pos_res))
    return id_pos_res


def read_sheet(path: str, offset: int = 0) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with open(path, encoding="utf-8") as handle:
        for _ in range(offset):
            next(handle)
        columns = handle.readline().strip().split("\t")
        for line_no, line in enumerate(handle, start=offset + 2):
            fields = line.strip().split("\t")
            if len(fields) > len(columns):
                raise NetworkinError(f"Too many columns at line {line_no} in {path}")
            row = {col: fields[i] if i < len(fields) else "" for i, col in enumerate(columns)}
            rows.append(row)
    logger.info("Loaded {} tabular rows from {}", len(rows), path)
    return rows


def insert_multilevel_list(mapping: dict[str, Any], keys: list[str], value: Any) -> None:
    cursor = mapping
    for key in keys[:-1]:
        cursor = cursor.setdefault(key, {})
    cursor.setdefault(keys[-1], []).append(value)


def set_multilevel_value(mapping: dict[str, Any], keys: list[str], value: Any) -> None:
    cursor = mapping
    for key in keys[:-1]:
        cursor = cursor.setdefault(key, {})
    previous = cursor.get(keys[-1])
    if previous is not None and previous != value:
        logger.warning("Replacing value at nested key path: {}", " -> ".join(keys))
    cursor[keys[-1]] = value


def read_group_to_domain_map(default_path: str) -> dict[str, dict[str, list[str]]]:
    mapping: dict[str, dict[str, list[str]]] = {}
    curated_path = Path(default_path).parent / "group_human_protein_name_map_curated.tsv"
    source = curated_path if curated_path.exists() else Path(default_path)

    with open(source, encoding="utf-8") as handle:
        for line in handle:
            tokens = line.split()
            if len(tokens) < 3:
                continue
            name = tokens[3] if source == curated_path and len(tokens) > 3 else tokens[2]
            insert_multilevel_list(mapping, tokens[:2], name)

    logger.success("Loaded group-to-domain map from {}", source)
    return mapping


def read_alias_files(
    organism: str, datadir: str, map_group_to_domain: dict[str, dict[str, list[str]]]
) -> tuple[dict[str, str], dict[str, str], dict[str, list[str]]]:
    alias_hash: dict[str, str] = {}
    desc_hash: dict[str, str] = {}
    name_hash: dict[str, list[str]] = {}

    allowed_names = {
        name
        for tree_map in map_group_to_domain.values()
        for pred_names in tree_map.values()
        for name in pred_names
    }

    desc_file = os.path.join(datadir, f"{organism}.text_best.v9.0.tsv.gz")
    alias_file = os.path.join(datadir, f"{organism}.protein.aliases.v12.0.txt.gz")

    try:
        for line in run_command(f"gzip -cd {desc_file}").splitlines():
            if not line.strip():
                continue
            tax_id, seq_id, desc = line.split("\t")[:3]
            desc_hash[seq_id] = desc
    except Exception as exc:
        logger.warning("Description aliases unavailable for {}: {}", organism, exc)

    try:
        for line in run_command(f"gzip -cd {alias_file}").splitlines():
            if not line.strip():
                continue
            seq_id, alias, source = line.split("\t")
            key = seq_id[5:]
            name_hash.setdefault(key, [])
            if alias not in name_hash[key]:
                name_hash[key].append(alias)
    except Exception as exc:
        logger.warning("Protein aliases unavailable for {}: {}", organism, exc)

    for seq_id, aliases in name_hash.items():
        alias_hash[seq_id] = next((a for a in aliases if a in allowed_names), seq_id)

    logger.success(
        "Loaded {} preferred aliases and {} descriptions", len(alias_hash), len(desc_hash)
    )
    return alias_hash, desc_hash, name_hash


def map_peptides_to_string(
    config: AppConfig, id_pos_res: dict[str, dict[int, str]], id_seq: dict[str, str]
) -> tuple[dict[str, dict[str, bool]], dict[str, dict[str, bool]]]:
    logger.header("Mapping query proteins to STRING IDs via BLAST")
    incoming2string: dict[str, dict[str, bool]] = {}
    string2incoming: dict[str, dict[str, bool]] = {}

    with tempfile.NamedTemporaryFile(mode="w", delete=False) as blast_tmpfile:
        iterable = id_pos_res.keys() if id_pos_res else id_seq.keys()
        for protein_id in iterable:
            seq = id_seq.get(protein_id)
            if not seq:
                logger.warning("No sequence available for {}", protein_id)
                continue
            blast_tmpfile.write(f">{protein_id}\n{seq}\n")
        query_path = blast_tmpfile.name

    blast_db = os.path.join(config.datadir, f"{config.organism}.protein.sequences.v12.0.fa")

    # --- BLAST PATH RESOLUTION ---
    resolved_blast_dir = config.blast_dir or os.environ.get("BLAST_PATH", "")
    blastp_exe = "blastp"
    makeblastdb_exe = "makeblastdb"

    if resolved_blast_dir:
        test_blastp = os.path.join(resolved_blast_dir, "blastp")
        # Check if the file actually exists and is executable
        if os.path.isfile(test_blastp) and os.access(test_blastp, os.X_OK):
            blastp_exe = test_blastp
            makeblastdb_exe = os.path.join(resolved_blast_dir, "makeblastdb")
        else:
            logger.warning(
                "blastp not found in '{}'. Falling back to system PATH.", resolved_blast_dir
            )

    # Final safety check: does blastp exist in the system at all?
    if not shutil.which(blastp_exe):
        raise CommandError(
            "blastp is not installed or not found in system PATH. Please install NCBI BLAST+ (e.g., sudo apt-get install ncbi-blast+)"
        )

    if not os.path.isfile(blast_db + ".pin"):
        command = f'{makeblastdb_exe} -in "{blast_db}" -out "{blast_db}" -parse_seqids -dbtype prot'
        logger.warning("BLAST DB missing. Initializing database.")
        run_command(command)

    command = (
        f"{blastp_exe} -num_threads {config.threads} -evalue 1e-10 "
        f'-db "{blast_db}" -query "{query_path}" -outfmt 6 | sort -k12gr'
    )

    if config.fast and os.path.isfile(config.blast_output_path):
        logger.info("Using cached BLAST output: {}", config.blast_output_path)
        with open(config.blast_output_path, encoding="utf-8") as handle:
            blast_lines = handle.read().splitlines()
    else:
        blast_lines = run_command(command).splitlines()
        if config.leave_intermediates:
            with open(config.blast_output_path, "w", encoding="utf-8") as handle:
                handle.write("\n".join(blast_lines))

    for line in blast_lines:
        if not line.strip():
            continue
        tokens = line.split("\t")
        incoming = tokens[0]
        if incoming in incoming2string:
            continue

        string_id = tokens[1].replace(f"{config.organism}.", "")
        identity = float(tokens[2])
        evalue = float(tokens[-2])

        if identity < 90:
            logger.warning("Best hit for {} has low identity: {:.2f}%", incoming, identity)
        if evalue > 1e-40:
            logger.warning("Best hit for {} has weak E-value: {:.2e}", incoming, evalue)
        if string_id in string2incoming:
            logger.warning("Best hit for {} is not reciprocal", incoming)

        incoming2string.setdefault(incoming, {})[string_id] = True
        string2incoming.setdefault(string_id, {})[incoming] = True

    logger.success("Mapped {} input proteins to STRING IDs", len(incoming2string))
    return incoming2string, string2incoming


def load_string_data(
    config: AppConfig, string2incoming: dict[str, dict[str, bool]], alias_hash: dict[str, str]
) -> dict[str, dict[str, dict[str, Any]]]:
    fn_bestpath = os.path.join(
        config.datadir, "string_data", f"{config.organism}.bestpath_0340_0950.v9.tsv.gz"
    )
    if not os.path.isfile(fn_bestpath):
        raise NetworkinError(f"Best path file missing: {fn_bestpath}")

    tree_pred_string_data: dict[str, dict[str, dict[str, Any]]] = {}
    edge_count = 0

    with gzip.open(fn_bestpath, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            edge_count += 1
            tokens = line.rstrip("\n").split("\t")

            parsed = _parse_string_line(tokens, config.organism, alias_hash)
            if parsed is None:
                continue

            name, string1, string2, direct_score, indirect_score, path = parsed
            if string2 not in string2incoming:
                continue

            tree_pred_string_data.setdefault(string2, {}).setdefault(string1, {"_name": name})
            score = direct_score if config.path_mode == "direct" else indirect_score
            tree_pred_string_data[string2][string1]["_score"] = score
            tree_pred_string_data[string2][string1]["_path"] = path

    logger.success("Loaded {} STRING best-path edges", edge_count)
    return tree_pred_string_data


def _parse_string_line(
    tokens: list[str], organism: str, alias_hash: dict[str, str]
) -> tuple[str, str, str, float, float, str] | None:
    if len(tokens) == 8:
        tree, group, name, string1, string2, s_direct, s_indirect, path = tokens
        return name, string1, string2, float(s_direct), float(s_indirect), path
    if len(tokens) == 7:
        tree, group, name, string1, string2, s_direct, s_indirect = tokens
        return name, string1, string2, float(s_direct), float(s_indirect), "notdef"
    if len(tokens) == 6:
        tree, group, name, string1, string2, s_direct = tokens
        string1 = string1.replace(f"{organism}.", "")
        string2 = string2.replace(f"{organism}.", "")
        score = float(s_direct) / 1000.0
        return name, string1, string2, score, score, "notdef"
    if len(tokens) == 3:
        string1, string2, s_direct = tokens
        string1 = string1.replace(f"{organism}.", "")
        string2 = string2.replace(f"{organism}.", "")
        score = float(s_direct) / 1000.0
        return alias_hash.get(string1, string1), string1, string2, score, score, "notdef"
    return None


def load_conversion_tables(parquet_path: str, species: str, verbose: bool = False) -> dict[str, Any]:
    """Load likelihood conversion tables from a Parquet file.

    Args:
        parquet_path: Path to the Parquet file produced by migrate_to_parquet.py.
        species: Species name to filter on (e.g. ``"human"`` or ``"yeast"``).
        verbose: Emit a log line for each table loaded.

    Returns:
        Nested dictionary ``tables[species][tree][player_name][score_kind]``
        where each leaf value is a ``list[CConvEntry]`` as returned by
        :func:`ReadConversionTableFromMemory`.
    """
    tables: dict[str, Any] = {}
    df = pd.read_parquet(parquet_path)
    species_df = df[df["species"] == species]

    for row in species_df.itertuples(index=False):
        conversion_tbl = ReadConversionTableFromMemory(row.raw_data)
        set_multilevel_value(
            tables, [row.species, row.tree, row.player_name, row.score_kind], conversion_tbl
        )
        if verbose:
            logger.muted(
                "Loaded conversion table: {} {} {} {}", row.species, row.tree, row.player_name, row.score_kind
            )

    logger.success("Loaded likelihood conversion tables for species {}", species)
    return tables


def _select_conversion_tables(
    dlr: dict[str, Any], species: str, tree: str, name: str
) -> tuple[Any, Any]:
    if species == "human" and tree in ["1433", "BRCT", "WW", "PTB", "WD40", "FHA"]:
        return dlr[species]["SH2"]["general"]["motif"], dlr[species]["SH2"]["general"]["string"]

    tree_tables = dlr[species][tree]
    if name in tree_tables:
        return tree_tables[name]["motif"], tree_tables[name]["string"]
    return tree_tables["general"]["motif"], tree_tables["general"]["string"]


def build_prediction_row(
    target_id: str,
    pos: int,
    residue: str,
    tree: str,
    pred: str,
    name: str,
    peptide: str,
    motif_score: float,
    string1: str,
    string2: str,
    string_score: float,
    path: str,
    networkin_score: float,
    string_alias: dict[str, str],
    string_desc: dict[str, str],
    recovered: bool = False,
    recovery_method: str = "",
) -> dict[str, Any]:
    return {
        "Name": target_id,
        "Position": f"{residue}{pos}" if pos else "",
        "Tree": tree,
        "Motif Group": pred,
        "Kinase/Phosphatase/Phospho-binding domain": name,
        "NetworKIN score": round(float(networkin_score), 4),
        "Motif probability": round(float(motif_score), 4),
        "STRING score": round(float(string_score), 4),
        "Target STRING ID": string1,
        "Kinase STRING ID": string2,
        "Target Name": string_alias.get(string1, string1),
        "Kinase Name": string_alias.get(string2, string2),
        "Target description": string_desc.get(string1, "notdef"),
        "Kinase description": string_desc.get(string2, "notdef"),
        "Peptide sequence window": peptide,
        "Intermediate nodes": path,
        "recovered": recovered,
        "recovery_method": recovery_method,
    }


def compile_predictions(
    config: AppConfig,
    id_pos_tree_pred: dict[str, Any],
    tree_pred_string_data: dict[str, dict[str, dict[str, Any]]],
    incoming2string: dict[str, dict[str, bool]],
    string_alias: dict[str, str],
    string_desc: dict[str, str],
    map_group_to_domain: dict[str, dict[str, list[str]]],
    likelihood_path: str,
    fasta_path: str,
) -> tuple[list[dict[str, Any]], PredictionStats]:
    stats = PredictionStats()
    predictions: list[dict[str, Any]] = []
    unmatched_names: set[str] = set()
    conversion_tables = load_conversion_tables(likelihood_path, config.species_name, config.verbose)

    for protein_id, pos_data in id_pos_tree_pred.items():
        stats.proteins_seen += 1
        mapped_strings = incoming2string.get(protein_id)
        if not mapped_strings:
            stats.proteins_unmapped += 1
            continue

        stats.proteins_mapped += 1
        for pos, tree_data in pos_data.items():
            for tree, pred_data in tree_data.items():
                for pred, (residue, peptide, motif_score) in pred_data.items():
                    allowed_names = set(map_group_to_domain.get(tree, {}).get(pred, []))
                    for string1 in mapped_strings:
                        neighbors = tree_pred_string_data.get(string1, {})
                        if not neighbors:
                            continue
                        stats.string_targets_seen += 1
                        for string2, edge_data in neighbors.items():
                            stats.string_scores_seen += 1
                            name = edge_data["_name"]
                            string_score = edge_data["_score"]
                            path = edge_data.get("_path", "notdef")

                            if name not in allowed_names:
                                unmatched_names.add(name)
                                if not config.string_for_uncovered:
                                    continue
                                motif_likelihood = 1.0
                                if config.species_name == "human" and tree in [
                                    "1433",
                                    "BRCT",
                                    "WW",
                                    "PTB",
                                    "WD40",
                                    "FHA",
                                ]:
                                    string_tbl = conversion_tables[config.species_name]["SH2"][
                                        "general"
                                    ]["string"]
                                else:
                                    string_tbl = conversion_tables[config.species_name][tree][
                                        "general"
                                    ]["string"]
                                string_likelihood = ConvertScore2L(string_score, string_tbl)
                                networkin_score = motif_likelihood * string_likelihood
                                stats.uncovered_branch_hits += 1
                            else:
                                motif_tbl, string_tbl = _select_conversion_tables(
                                    conversion_tables, config.species_name, tree, name
                                )
                                motif_likelihood = ConvertScore2L(motif_score, motif_tbl)
                                string_likelihood = ConvertScore2L(string_score, string_tbl)
                                networkin_score = motif_likelihood * string_likelihood
                                stats.covered_branch_hits += 1

                            if networkin_score < 0.02:
                                continue

                            row = build_prediction_row(
                                target_id=protein_id,
                                pos=pos,
                                residue=residue,
                                tree=tree,
                                pred=pred,
                                name=name,
                                peptide=peptide,
                                motif_score=motif_score,
                                string1=string1,
                                string2=string2,
                                string_score=string_score,
                                path=path,
                                networkin_score=networkin_score,
                                string_alias=string_alias,
                                string_desc=string_desc,
                            )
                            predictions.append(row)
                            stats.result_rows += 1

    stats.unmatched_names = len(unmatched_names)
    if unmatched_names:
        logger.warning("Unmatched kinase/domain names: {}", len(unmatched_names))
    logger.success("Compiled {} raw prediction rows", len(predictions))
    return predictions, stats


def recover_predictions(
    id_pos_tree_pred: dict[str, Any],
    incoming2string: dict[str, dict[str, bool]],
    tree_pred_string_data: dict[str, dict[str, dict[str, Any]]],
    string_alias: dict[str, str],
    string_desc: dict[str, str],
) -> list[dict[str, Any]]:
    motif_score_dict: dict[tuple[str, str], float] = {}
    for prot_id, pos_data in id_pos_tree_pred.items():
        mapped_strings = incoming2string.get(prot_id)
        if not mapped_strings:
            continue
        for string1 in mapped_strings:
            if string1 not in tree_pred_string_data:
                continue
            for _pos, tree_data in pos_data.items():
                for _tree, pred_data in tree_data.items():
                    for _pred_name, (_res, _peptide, motif_score) in pred_data.items():
                        for string2 in tree_pred_string_data[string1]:
                            key = (string2, string1)
                            motif_score_dict[key] = max(
                                motif_score, motif_score_dict.get(key, -1.0)
                            )

    all_pairs: list[tuple[str, str]] = []
    protein_ids: set[str] = set()
    for sub, kin_map in tree_pred_string_data.items():
        protein_ids.add(sub)
        for kin in kin_map:
            protein_ids.add(kin)
            all_pairs.append((kin, sub))

    protein_list = sorted(protein_ids)
    node_index = {pid: i for i, pid in enumerate(protein_list)}
    n = len(protein_list)
    dist_matrix = np.full((n, n), np.inf, dtype=np.float32)
    np.fill_diagonal(dist_matrix, 0.0)

    for sub, kin_map in tree_pred_string_data.items():
        si = node_index[sub]
        for kin, data in kin_map.items():
            ki = node_index[kin]
            score = data.get("_score", 0.0)
            if score > 0:
                dist_matrix[ki, si] = min(dist_matrix[ki, si], float(1.0 / score - 1.0))

    recovered = recover_false_negatives(
        candidates=all_pairs,
        dist_matrix=dist_matrix,
        node_index=node_index,
        motif_scores=motif_score_dict,
    )

    recovered_rows = [
        {
            "Name": string_alias.get(r["substrate_uniprot"], r["substrate_uniprot"]),
            "Position": "",
            "Tree": "",
            "Motif Group": "",
            "Kinase/Phosphatase/Phospho-binding domain": string_alias.get(
                r["kinase_id"], r["kinase_id"]
            ),
            "NetworKIN score": round(float(r["networkin_score"]), 4),
            "Motif probability": round(float(r["motif_score"]), 4),
            "STRING score": round(float(r["context_score"]), 4),
            "Target STRING ID": r["substrate_uniprot"],
            "Kinase STRING ID": r["kinase_id"],
            "Target Name": string_alias.get(r["substrate_uniprot"], r["substrate_uniprot"]),
            "Kinase Name": string_alias.get(r["kinase_id"], r["kinase_id"]),
            "Target description": string_desc.get(r["substrate_uniprot"], "notdef"),
            "Kinase description": string_desc.get(r["kinase_id"], "notdef"),
            "Peptide sequence window": "",
            "Intermediate nodes": "notdef",
            "recovered": r["recovered"],
            "recovery_method": r["recovery_method"],
        }
        for r in recovered
    ]

    logger.success("Recovered {} false-negative predictions", len(recovered_rows))
    return recovered_rows


def run_pipeline(config: AppConfig) -> dict[str, Any]:
    ensure_dirs(config.result_dir, config.temp_dir)
    logger.header("Starting NetworKIN pipeline")
    logger.info("Organism: {}", config.organism)
    logger.info("FASTA: {}", config.fasta_path)
    logger.info("Sites: {}", config.sites_path or "<all S/T/Y residues>")

    id_seq = read_fasta_file(config.fasta_path)

    if config.sites_path:
        file_type = detect_site_file_type(config.sites_path)
        if file_type == NETWORKIN_SITE_FILE:
            id_pos_res = read_networkin_sites(config.sites_path)
        elif file_type == PROTEOME_DISCOVERER_SITE_FILE:
            id_pos_res = read_proteome_discoverer_sites(config.sites_path)
        elif file_type == MAX_QUANT_DIRECT_OUTPUT_FILE:
            id_pos_res = read_max_quant_sites(config.sites_path)
        else:
            # handle LEGACY_SITE_FILE and MS_MCMC later...
            raise NetworkinError(
                "Only the basic Networkin site reader was migrated in this core refactor. Add the other adapters next."
            )
    else:
        id_pos_res = {}
        logger.info(
            "No sites file supplied; downstream scorer must enumerate all candidate S/T/Y sites"
        )

    logger.info("Fetching live phosphorylation reference data")
    fetch_phosphosite(refresh=config.refresh)

    group_map_path = (
        os.path.join(config.datadir, "group_human_protein_name_map.tsv")
        if config.organism == "9606"
        else os.path.join(config.datadir, "group_yeast_KIN.tsv")
    )
    map_group_to_domain = read_group_to_domain_map(group_map_path)
    string_alias, string_desc, name_hash = read_alias_files(
        config.organism, config.datadir, map_group_to_domain
    )
    incoming2string, string2incoming = map_peptides_to_string(config, id_pos_res, id_seq)

    logger.info("Fetching live STRING associations for reference")
    fetch_string_network(proteins=list(id_seq.keys()), refresh=config.refresh)

    tree_pred_string_data = load_string_data(config, string2incoming, string_alias)
    logger.info("Running motif scorer")
    id_pos_tree_pred = score_sequences(id_seq, id_pos_res)

    with conversion_parquet_path(config.path_mode) as likelihood_path:
        predictions, stats = compile_predictions(
            config=config,
            id_pos_tree_pred=id_pos_tree_pred,
            tree_pred_string_data=tree_pred_string_data,
            incoming2string=incoming2string,
            string_alias=string_alias,
            string_desc=string_desc,
            map_group_to_domain=map_group_to_domain,
            likelihood_path=str(likelihood_path),
            fasta_path=config.fasta_path,
        )

    predictions.extend(
        recover_predictions(
            id_pos_tree_pred=id_pos_tree_pred,
            incoming2string=incoming2string,
            tree_pred_string_data=tree_pred_string_data,
            string_alias=string_alias,
            string_desc=string_desc,
        )
    )

    final_predictions = filter_and_rank_predictions(predictions)
    output_path = os.path.join(config.result_dir, f"{Path(config.fasta_path).name}.result.tsv")
    write_output(final_predictions, output_path)

    logger.success("Wrote {} final predictions to {}", len(final_predictions), output_path)
    logger.header("Pipeline complete")
    for key, value in stats.as_dict().items():
        logger.score("{} = {}", key, value)

    return {
        "output_path": output_path,
        "stats": stats.as_dict(),
        "prediction_count": len(final_predictions),
    }