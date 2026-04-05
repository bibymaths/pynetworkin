# inputs/maxquant_processor.py
#
# Processes MaxQuant Site Table (Phospho (STY)Sites.txt) output for
# PyNetworKIN integration.
#
# Steps:
#   1. Parse the MaxQuant site table (tab-separated).
#   2. Clean & classify every protein ID in the "Proteins" column.
#   3. Auto-detect ID type (UniProt, RefSeq, Ensembl, other).
#   4. Download FASTA sequences:
#      - Ensembl protein IDs (ENS*P...) are fetched DIRECTLY from the
#        Ensembl REST sequence endpoint – no UniProt mapping required.
#      - UniProt accessions are fetched directly from UniProt.
#      - RefSeq IDs are still mapped via UniProt ID-mapping before fetching.
#   5. Write a clean FASTA file, a cleaned site table, an ID-mapping CSV,
#      and a JSON processing report to the requested output directory.
#
# Ensembl REST endpoint (protein FASTA):
#   https://rest.ensembl.org/sequence/id/{ID}?type=protein
#
# UniProt REST endpoints:
#   - Fetch by accession: https://rest.uniprot.org/uniprotkb/{ID}.fasta
#   - ID-mapping service:  https://rest.uniprot.org/idmapping/
#
# Rate limiting: honour rate_limit_delay (seconds) between individual
# requests; batch up to UNIPROT_BATCH_SIZE accessions into one search
# request when fetching UniProt IDs directly.

from __future__ import annotations

import asyncio
import json
import logging
import re
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
from urllib.parse import quote

import httpx
import pandas as pd

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_IDMAPPING_RUN_URL = "https://rest.uniprot.org/idmapping/run"
UNIPROT_IDMAPPING_STATUS_URL = "https://rest.uniprot.org/idmapping/status/{job_id}"
UNIPROT_IDMAPPING_DETAILS_URL = "https://rest.uniprot.org/idmapping/details/{job_id}"

# Ensembl protein sequence endpoint – used directly for ENS*P... IDs so
# that UniProt ID-mapping is no longer needed for Ensembl accessions.
ENSEMBL_SEQUENCE_URL = "https://rest.ensembl.org/sequence/id/{id}?type=protein"

UNIPROT_BATCH_SIZE = 50   # max accessions per search request
MAX_RETRIES = 3
RETRY_DELAY = 2.0          # base delay in seconds (doubles each retry)
IDMAPPING_POLL_INTERVAL = 2.0
IDMAPPING_POLL_MAX = 30    # seconds before giving up on a job


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------


@dataclass
class ProteinIDInfo:
    """Parsed and classified protein identifier from MaxQuant output."""

    original_id: str
    cleaned_id: str
    id_type: str          # 'uniprot', 'refseq', 'ensembl', 'other'
    is_contaminant: bool
    is_reverse: bool
    isoform: str | None = None


@dataclass
class ProcessingReport:
    """Statistics produced by :meth:`MaxQuantProcessor.process_site_table`."""

    input_file: str = ""
    total_rows: int = 0
    total_protein_entries: int = 0
    unique_original_ids: int = 0
    unique_cleaned_ids: int = 0
    contaminants_removed: int = 0
    reverse_removed: int = 0
    uniprot_ids: int = 0
    refseq_ids: int = 0
    ensembl_ids: int = 0
    other_ids: int = 0
    sequences_downloaded: int = 0
    sequences_missing: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "input_file": self.input_file,
            "total_rows": self.total_rows,
            "total_protein_entries": self.total_protein_entries,
            "unique_original_ids": self.unique_original_ids,
            "unique_cleaned_ids": self.unique_cleaned_ids,
            "contaminants_removed": self.contaminants_removed,
            "reverse_removed": self.reverse_removed,
            "id_types": {
                "uniprot": self.uniprot_ids,
                "refseq": self.refseq_ids,
                "ensembl": self.ensembl_ids,
                "other": self.other_ids,
            },
            "sequences_downloaded": self.sequences_downloaded,
            "sequences_missing": self.sequences_missing,
            "errors": self.errors,
        }


# ---------------------------------------------------------------------------
# Core processor
# ---------------------------------------------------------------------------


class MaxQuantProcessor:
    """Process MaxQuant Site Table files for PyNetworKIN compatibility.

    Parameters
    ----------
    rate_limit_delay:
        Minimum pause (seconds) between UniProt HTTP requests.
    """

    # Pre-compiled ID patterns
    _RE_SWISSPROT = re.compile(r"sp\|([A-Z0-9]+)\|")
    _RE_TREMBL = re.compile(r"tr\|([A-Z0-9]+)\|")
    # Bare UniProt accession: [A-N,R-Z][0-9][A-Z][A-Z0-9]{2}[0-9] (6-char) or
    # [O,P,Q][0-9][A-Z0-9]{3}[0-9] (6-char) or 10-char secondary accessions
    _RE_UNIPROT_BARE = re.compile(
        r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"
    )
    _RE_REFSEQ = re.compile(r"^(NP_|XP_|WP_)\d+")
    _RE_ENSEMBL = re.compile(r"^ENS[A-Z]*P\d{11}")
    _RE_CONTAMINANT = re.compile(r"^CON__")
    _RE_REVERSE = re.compile(r"^REV__")
    _RE_ISOFORM = re.compile(r"-(\d+)$")
    _RE_REFSEQ_VERSION = re.compile(r"\.\d+$")

    def __init__(self, rate_limit_delay: float = 0.1) -> None:
        self.rate_limit_delay = rate_limit_delay

    # ------------------------------------------------------------------
    # ID cleaning helpers
    # ------------------------------------------------------------------

    def _strip_prefix(self, raw: str) -> tuple[str, bool, bool]:
        """Return (stripped_id, is_contaminant, is_reverse)."""
        is_contaminant = bool(self._RE_CONTAMINANT.match(raw))
        is_reverse = bool(self._RE_REVERSE.match(raw))
        stripped = self._RE_CONTAMINANT.sub("", raw)
        stripped = self._RE_REVERSE.sub("", stripped)
        # Remove sp| / tr| wrappers: sp|ACC|ENTRY → ACC
        m = self._RE_SWISSPROT.match(stripped) or self._RE_TREMBL.match(stripped)
        if m:
            stripped = m.group(1)
        return stripped, is_contaminant, is_reverse

    def _classify(self, acc: str) -> tuple[str, str, str | None]:
        """Classify and normalise a bare accession string.

        Strips isoform suffixes and RefSeq version numbers, then detects the
        ID type.

        Parameters
        ----------
        acc:
            Bare accession (after prefix removal by :meth:`_strip_prefix`).

        Returns
        -------
        tuple of (id_type, cleaned_accession, isoform_or_None)
            ``id_type`` is one of ``'uniprot'``, ``'refseq'``, ``'ensembl'``,
            or ``'other'``.  ``cleaned_accession`` has isoform and version
            suffixes removed.
        """
        isoform: str | None = None
        m_iso = self._RE_ISOFORM.search(acc)
        if m_iso:
            isoform = m_iso.group(1)
            acc = acc[: m_iso.start()]

        # Remove RefSeq version suffix (NP_001234.1 → NP_001234)
        if self._RE_REFSEQ.match(acc):
            acc = self._RE_REFSEQ_VERSION.sub("", acc)
            return "refseq", acc, isoform

        if self._RE_ENSEMBL.match(acc):
            return "ensembl", acc, isoform

        if self._RE_UNIPROT_BARE.match(acc):
            return "uniprot", acc, isoform

        return "other", acc, isoform

    def clean_protein_ids(self, proteins_str: str) -> list[ProteinIDInfo]:
        """Clean and classify protein IDs from a MaxQuant *Proteins* cell.

        Parameters
        ----------
        proteins_str:
            Semicolon-separated protein IDs as they appear in MaxQuant output.

        Returns
        -------
        list[ProteinIDInfo]
            One entry per ID token (including contaminants and reversed hits so
            callers can decide what to keep).
        """
        results: list[ProteinIDInfo] = []
        for raw in proteins_str.split(";"):
            raw = raw.strip()
            if not raw:
                continue

            stripped, is_contaminant, is_reverse = self._strip_prefix(raw)
            id_type, cleaned_id, isoform = self._classify(stripped)

            results.append(
                ProteinIDInfo(
                    original_id=raw,
                    cleaned_id=cleaned_id,
                    id_type=id_type,
                    is_contaminant=is_contaminant,
                    is_reverse=is_reverse,
                    isoform=isoform,
                )
            )

        return results

    def _clean_proteins_column(self, proteins_str: str) -> str:
        """Return a normalised semicolon-joined string of cleaned IDs.

        Contaminants and reversed hits are excluded; duplicates within a
        cell are collapsed.
        """
        infos = self.clean_protein_ids(proteins_str)
        seen: set[str] = set()
        cleaned: list[str] = []
        for info in infos:
            if info.is_contaminant or info.is_reverse:
                continue
            if info.cleaned_id not in seen:
                seen.add(info.cleaned_id)
                cleaned.append(info.cleaned_id)
        return ";".join(cleaned)

    # ------------------------------------------------------------------
    # UniProt download
    # ------------------------------------------------------------------

    async def _fetch_one_fasta(
        self,
        client: httpx.AsyncClient,
        accession: str,
    ) -> str | None:
        """Fetch a single UniProt FASTA; return raw FASTA text or None."""
        url = UNIPROT_FASTA_URL.format(accession=quote(accession, safe=""))
        for attempt in range(MAX_RETRIES):
            try:
                resp = await client.get(url, timeout=30)
                if resp.status_code == 200:
                    text = resp.text.strip()
                    return text if text.startswith(">") else None
                if resp.status_code == 429:
                    # Rate limited – back off
                    await asyncio.sleep(RETRY_DELAY * (attempt + 1))
                    continue
                if resp.status_code == 404:
                    return None
                resp.raise_for_status()
            except httpx.TimeoutException:
                logger.warning("Timeout fetching %s (attempt %d)", accession, attempt + 1)
                await asyncio.sleep(RETRY_DELAY * (attempt + 1))
            except httpx.HTTPStatusError as exc:
                logger.warning("HTTP error fetching %s: %s", accession, exc)
                await asyncio.sleep(RETRY_DELAY * (attempt + 1))
        return None

    async def _fetch_ensembl_fasta(
        self,
        client: httpx.AsyncClient,
        ensembl_id: str,
    ) -> str | None:
        """Fetch a single Ensembl protein FASTA directly from Ensembl REST API.

        Ensembl protein IDs (``ENS*P...``) are resolved here without going
        through UniProt ID-mapping, which avoids stalling on that service.

        Parameters
        ----------
        client:
            Shared async HTTP client.
        ensembl_id:
            Cleaned Ensembl protein accession (e.g. ``ENSP00000449404``).

        Returns
        -------
        Raw FASTA text (header + sequence lines) or ``None`` if not found.
        """
        url = ENSEMBL_SEQUENCE_URL.format(id=quote(ensembl_id, safe=""))
        for attempt in range(MAX_RETRIES):
            try:
                resp = await client.get(
                    url,
                    headers={"Accept": "text/x-fasta"},
                    timeout=30,
                )
                if resp.status_code == 200:
                    text = resp.text.strip()
                    return text if text.startswith(">") else None
                if resp.status_code == 429:
                    # Rate limited – back off
                    await asyncio.sleep(RETRY_DELAY * (attempt + 1))
                    continue
                if resp.status_code == 404:
                    return None
                resp.raise_for_status()
            except httpx.TimeoutException:
                logger.warning(
                    "Timeout fetching Ensembl %s (attempt %d)", ensembl_id, attempt + 1
                )
                await asyncio.sleep(RETRY_DELAY * (attempt + 1))
            except httpx.HTTPStatusError as exc:
                logger.warning("HTTP error fetching Ensembl %s: %s", ensembl_id, exc)
                await asyncio.sleep(RETRY_DELAY * (attempt + 1))
        return None

    async def _submit_idmapping_job(
        self,
        client: httpx.AsyncClient,
        ids: list[str],
        from_db: str,
        to_db: str = "UniProtKB",
    ) -> str | None:
        """Submit an ID-mapping job; return job ID or None on failure."""
        try:
            resp = await client.post(
                UNIPROT_IDMAPPING_RUN_URL,
                data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
                timeout=30,
            )
            resp.raise_for_status()
            data = resp.json()
            return data.get("jobId")
        except Exception as exc:
            logger.warning("ID-mapping submission failed: %s", exc)
            return None

    async def _get_idmapping_results_url(
        self,
        client: httpx.AsyncClient,
        job_id: str,
    ) -> str | None:
        """Return the UniProt redirect URL for a finished ID-mapping job."""
        resp = await client.get(
            UNIPROT_IDMAPPING_DETAILS_URL.format(job_id=job_id),
            timeout=30,
        )
        resp.raise_for_status()
        payload = resp.json()
        return payload.get("redirectURL")

    async def _poll_idmapping_results(
        self,
        client: httpx.AsyncClient,
        job_id: str,
    ) -> list[dict[str, str]]:
        """Poll until an ID-mapping job finishes; return result list."""
        deadline = time.monotonic() + IDMAPPING_POLL_MAX
        status_url = UNIPROT_IDMAPPING_STATUS_URL.format(job_id=job_id)

        while time.monotonic() < deadline:
            await asyncio.sleep(IDMAPPING_POLL_INTERVAL)
            try:
                status_resp = await client.get(status_url, timeout=15)
                status_resp.raise_for_status()
                payload = status_resp.json()

                # UniProt may return {"jobStatus": "RUNNING"} while the job is
                # still in progress, but once the job is ready the same endpoint
                # may already return the results payload directly.
                job_status = payload.get("jobStatus")
                if job_status in {"NEW", "RUNNING"}:
                    continue
                if job_status in {"ERROR", "FAILURE", "NOT_FOUND"}:
                    logger.warning("ID-mapping job %s failed with status %s", job_id, job_status)
                    return []
                if "results" in payload or "failedIds" in payload:
                    return payload.get("results", [])

                results_url = await self._get_idmapping_results_url(client, job_id)
                if not results_url:
                    logger.warning("No redirect URL returned for ID-mapping job %s", job_id)
                    return []

                results_resp = await client.get(results_url, timeout=60)
                results_resp.raise_for_status()
                results_data = results_resp.json()
                return results_data.get("results", [])
            except Exception as exc:
                logger.warning("Error polling ID-mapping job %s: %s", job_id, exc)

        logger.warning("ID-mapping job %s timed out", job_id)
        return []

    async def _map_non_uniprot_ids(
        self,
        client: httpx.AsyncClient,
        refseq_ids: list[str],
        ensembl_ids: list[str],
    ) -> dict[str, str]:
        """Map RefSeq and Ensembl IDs to UniProt accessions.

        Returns
        -------
        dict mapping original non-UniProt ID → UniProt accession.
        """
        mapped: dict[str, str] = {}

        for batch_ids, db_name in [
            (refseq_ids, "RefSeq_Protein"),
            (ensembl_ids, "Ensembl_Protein"),
        ]:
            for i in range(0, len(batch_ids), UNIPROT_BATCH_SIZE):
                batch = batch_ids[i : i + UNIPROT_BATCH_SIZE]
                job_id = await self._submit_idmapping_job(client, batch, from_db=db_name)
                if not job_id:
                    continue
                results = await self._poll_idmapping_results(client, job_id)
                for r in results:
                    orig = r.get("from", "")
                    entry = r.get("to", {})
                    if isinstance(entry, dict):
                        acc = entry.get("primaryAccession", "")
                    else:
                        acc = str(entry)
                    if orig and acc:
                        mapped[orig] = acc
                await asyncio.sleep(self.rate_limit_delay)

        return mapped

    async def fetch_sequences_from_uniprot(
        self,
        protein_ids: list[str],
        id_types: dict[str, str] | None = None,
    ) -> dict[str, str]:
        """Fetch FASTA sequences for a list of protein IDs.

        Ensembl protein IDs (``ENS*P...``) are fetched directly from the
        Ensembl REST sequence endpoint – the UniProt ID-mapping step is
        bypassed entirely for these IDs.  RefSeq IDs are still mapped via
        UniProt ID-mapping before fetching.  Plain UniProt accessions are
        fetched directly from UniProt.

        Parameters
        ----------
        protein_ids:
            List of cleaned protein accessions.
        id_types:
            Optional mapping of accession → id_type ('uniprot', 'refseq',
            'ensembl', 'other').  When provided, non-UniProt IDs are handled
            according to their type (see above).

        Returns
        -------
        dict mapping accession → FASTA text (header + sequence lines).
        """
        if id_types is None:
            id_types = {}

        sequences: dict[str, str] = {}

        async with httpx.AsyncClient(follow_redirects=True) as client:
            # ── 1. Fetch Ensembl sequences DIRECTLY from Ensembl REST API ─────
            # Ensembl protein IDs are resolved here without UniProt mapping.
            ensembl_ids = [p for p in protein_ids if id_types.get(p) == "ensembl"]
            if ensembl_ids:
                logger.info(
                    "Fetching %d Ensembl sequence(s) directly from Ensembl REST API",
                    len(ensembl_ids),
                )
                for idx, eid in enumerate(ensembl_ids):
                    fasta = await self._fetch_ensembl_fasta(client, eid)
                    if fasta:
                        sequences[eid] = fasta
                        logger.debug("Fetched Ensembl sequence: %s", eid)
                    else:
                        logger.warning("No Ensembl sequence found for %s", eid)
                    # Ensembl REST API does not support batching, so we rate-limit
                    # between individual requests (unlike the batched UniProt path).
                    if idx < len(ensembl_ids) - 1:
                        await asyncio.sleep(self.rate_limit_delay)
                logger.info(
                    "Ensembl fetch complete: %d downloaded, %d missing",
                    sum(1 for eid in ensembl_ids if eid in sequences),
                    sum(1 for eid in ensembl_ids if eid not in sequences),
                )

            # ── 2. Map RefSeq IDs to UniProt (Ensembl excluded) ───────────────
            refseq = [p for p in protein_ids if id_types.get(p) == "refseq"]
            id_map: dict[str, str] = {}  # original → uniprot accession
            if refseq:
                id_map = await self._map_non_uniprot_ids(client, refseq, [])

            # Build the set of UniProt accessions to fetch
            fetch_targets: dict[str, str] = {}  # uniprot_acc → original_id
            for pid in protein_ids:
                ptype = id_types.get(pid, "uniprot")
                if ptype == "ensembl":
                    continue  # already fetched above
                if ptype == "refseq":
                    uniprot_acc = id_map.get(pid)
                    if uniprot_acc:
                        fetch_targets[uniprot_acc] = pid
                    else:
                        logger.warning("No UniProt mapping found for %s", pid)
                elif ptype == "other":
                    logger.warning("Skipping unsupported ID type for: %s", pid)
                else:
                    fetch_targets[pid] = pid

            # ── 3. Fetch FASTA in batches from UniProt ────────────────────────
            accessions = list(fetch_targets.keys())
            for i in range(0, len(accessions), UNIPROT_BATCH_SIZE):
                batch = accessions[i : i + UNIPROT_BATCH_SIZE]
                tasks = [self._fetch_one_fasta(client, acc) for acc in batch]
                raw_results = await asyncio.gather(*tasks)
                for acc, fasta_text in zip(batch, raw_results, strict=True):
                    original = fetch_targets[acc]
                    if fasta_text:
                        sequences[original] = fasta_text
                    else:
                        logger.warning("No sequence found for %s", acc)
                if i + UNIPROT_BATCH_SIZE < len(accessions):
                    await asyncio.sleep(self.rate_limit_delay)

        return sequences

    # ------------------------------------------------------------------
    # Main processing entry point
    # ------------------------------------------------------------------

    def process_site_table(
        self,
        input_file: Path,
        output_dir: Path,
    ) -> dict[str, Any]:
        """Process a MaxQuant Phospho (STY)Sites.txt file.

        Parameters
        ----------
        input_file:
            Path to the MaxQuant site table.
        output_dir:
            Directory where output files are written (created if absent).

        Returns
        -------
        dict with keys:
            - ``fasta_path``      : Path to the cleaned FASTA file
            - ``sites_path``      : Path to the cleaned site table
            - ``mapping_path``    : Path to the ID-mapping CSV
            - ``report_path``     : Path to the JSON report
            - ``report``          : :class:`ProcessingReport` instance
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        report = ProcessingReport(input_file=str(input_file))

        logger.info("Reading site table: %s", input_file)
        try:
            df = pd.read_csv(input_file, sep="\t", low_memory=False)
        except Exception as exc:
            msg = f"Failed to read input file: {exc}"
            logger.error(msg)
            report.errors.append(msg)
            return self._write_outputs(output_dir, pd.DataFrame(), {}, {}, report)

        report.total_rows = len(df)
        logger.info("Loaded %d rows", report.total_rows)

        if "Proteins" not in df.columns:
            msg = "'Proteins' column not found in site table"
            logger.error(msg)
            report.errors.append(msg)
            return self._write_outputs(output_dir, df, {}, {}, report)

        # ── 1. Clean and classify all IDs ─────────────────────────────────────
        all_infos: list[ProteinIDInfo] = []
        for cell in df["Proteins"].dropna():
            all_infos.extend(self.clean_protein_ids(str(cell)))

        report.total_protein_entries = len(all_infos)
        report.unique_original_ids = len({i.original_id for i in all_infos})
        report.contaminants_removed = sum(1 for i in all_infos if i.is_contaminant)
        report.reverse_removed = sum(1 for i in all_infos if i.is_reverse)

        # Keep only non-contaminant, non-reverse entries; de-duplicate
        kept = {
            i.cleaned_id: i
            for i in all_infos
            if not i.is_contaminant and not i.is_reverse and i.cleaned_id
        }
        report.unique_cleaned_ids = len(kept)
        report.uniprot_ids = sum(1 for i in kept.values() if i.id_type == "uniprot")
        report.refseq_ids = sum(1 for i in kept.values() if i.id_type == "refseq")
        report.ensembl_ids = sum(1 for i in kept.values() if i.id_type == "ensembl")
        report.other_ids = sum(1 for i in kept.values() if i.id_type == "other")

        logger.info(
            "Unique IDs after cleaning: %d  (UniProt=%d  RefSeq=%d  Ensembl=%d  other=%d)",
            report.unique_cleaned_ids,
            report.uniprot_ids,
            report.refseq_ids,
            report.ensembl_ids,
            report.other_ids,
        )

        # ── 2. Normalise Proteins column in site table ─────────────────────────
        df = df.copy()
        df["Proteins"] = df["Proteins"].apply(
            lambda x: self._clean_proteins_column(str(x)) if pd.notna(x) else ""
        )

        # ── 3. Build ID-mapping table ──────────────────────────────────────────
        # original_id → cleaned_id (all entries, including contaminants)
        mapping: dict[str, str] = {i.original_id: i.cleaned_id for i in all_infos}

        # ── 4. Download FASTA sequences ────────────────────────────────────────
        id_types = {cid: info.id_type for cid, info in kept.items()}
        sequences: dict[str, str] = {}
        try:
            sequences = asyncio.run(
                self.fetch_sequences_from_uniprot(
                    list(kept.keys()),
                    id_types=id_types,
                )
            )
        except Exception as exc:
            msg = f"Error during sequence download: {exc}"
            logger.error(msg)
            report.errors.append(msg)

        report.sequences_downloaded = len(sequences)
        report.sequences_missing = [
            cid for cid in kept if cid not in sequences
        ]

        logger.info(
            "Downloaded %d sequences; %d missing",
            report.sequences_downloaded,
            len(report.sequences_missing),
        )

        return self._write_outputs(output_dir, df, sequences, mapping, report)

    # ------------------------------------------------------------------
    # File writers
    # ------------------------------------------------------------------

    def _write_outputs(
        self,
        output_dir: Path,
        df: pd.DataFrame,
        sequences: dict[str, str],
        mapping: dict[str, str],
        report: ProcessingReport,
    ) -> dict[str, Any]:
        """Write all four output files and return paths."""
        fasta_path = output_dir / "cleaned_proteins.fasta"
        sites_path = output_dir / "cleaned_sites.txt"
        mapping_path = output_dir / "id_mapping.csv"
        report_path = output_dir / "processing_report.json"

        # FASTA
        with fasta_path.open("w", encoding="utf-8") as fh:
            for acc, fasta_text in sequences.items():
                lines = fasta_text.strip().splitlines()
                if not lines:
                    continue
                fh.write(f">{acc}\n")
                for seq_line in lines[1:]:
                    fh.write(seq_line.strip() + "\n")
        logger.info("Wrote FASTA: %s", fasta_path)

        # Cleaned site table
        if not df.empty:
            df.to_csv(sites_path, sep="\t", index=False)
        else:
            sites_path.touch()
        logger.info("Wrote cleaned site table: %s", sites_path)

        # ID mapping CSV
        mapping_rows = [
            {"original_id": orig, "cleaned_id": cleaned}
            for orig, cleaned in mapping.items()
        ]
        mapping_df = pd.DataFrame(mapping_rows, columns=["original_id", "cleaned_id"])
        mapping_df.to_csv(mapping_path, index=False)
        logger.info("Wrote ID mapping: %s", mapping_path)

        # Processing report JSON
        with report_path.open("w", encoding="utf-8") as fh:
            json.dump(report.to_dict(), fh, indent=2)
        logger.info("Wrote processing report: %s", report_path)

        return {
            "fasta_path": fasta_path,
            "sites_path": sites_path,
            "mapping_path": mapping_path,
            "report_path": report_path,
            "report": report,
        }