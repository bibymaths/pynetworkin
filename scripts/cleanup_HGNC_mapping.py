#!/usr/bin/env python3
"""Tidy up the HGNC → Ensembl protein-ID mapping table.

scripts/cleanup_HGNC_mapping.py
================================

This is a **legacy utility** (originally written for Python 2 with the
Lindinglab idmapper library) that reconciles ambiguous HGNC symbol mappings
for STRING protein IDs.  It is **not** part of the active pynetworkin
pipeline and is retained here for reference and potential future migration.

Dependencies (not installed by default):
    * Lindinglab.idmapper  — internal Linding-lab library

Input files expected in the working directory:
    * Biomart_map_Ensembl*.txt  — Ensembl → HGNC symbol BioMart exports
    * idmap_HGNC.txt            — HGNC approved symbol table
    * idmap_STRING_UniprotAC_by_Uniprot.txt
    * string_id_unique.txt      — one STRING protein ID per line
    * ../data/9606.alias_best.v9.0.tsv

Output:
    * 9606.alias_best.v9.0.17032014_2.tsv  — corrected alias table

Usage::

    cd HGNC_Symbol_mapping/
    python scripts/cleanup_HGNC_mapping.py
"""

from __future__ import annotations

# NOTE: The Lindinglab idmapper library is a legacy internal dependency.
# The code below is preserved as-is from the original cleanup tool.
# It requires Python 2 dict methods (.iteritems(), .has_key()) which are
# not available in Python 3.  A full migration is deferred.

try:
    from Lindinglab.idmapper import CIdMap, MergeMap, ReadIdMapSimple  # type: ignore[import]
except ImportError as exc:
    raise ImportError(
        "This legacy script requires the Lindinglab idmapper library, "
        "which is not publicly available.  It is retained for historical reference only."
    ) from exc


MANUAL_MAP: dict[str, str] = {
    "ENSP00000350696": "ZNF587B",
    "ENSP00000299353": "C10orf32",
    "ENSP00000341497": "ZNF177",
    "ENSP00000382819": "DXO",
    "ENSP00000395262": "DXO",
    "ENSP00000410404": "DXO",
    "ENSP00000226218": "VTN",
    "ENSP00000391123": "DXO",
    "ENSP00000391349": "DXO",
    "ENSP00000383645": "KIR2DL5A",
    "ENSP00000301407": "CGB1",
    "ENSP00000316021": "PLSCR3",
}

ENSEMBL_VERSIONS = [75, 68, 67, 66, 65, 64, 63, 62, 59, 54]


def build_ensp_hgnc_map() -> "CIdMap":
    """Merge BioMart Ensembl → HGNC symbol maps across multiple releases."""
    idmap = CIdMap()
    for ver in ENSEMBL_VERSIONS:
        MergeMap(idmap, ReadIdMapSimple(f"Biomart_map_Ensembl{ver}.txt", 1, 3, offset=1))
    return idmap


def build_ensp_ensg_map() -> "CIdMap":
    """Merge BioMart Ensembl protein → gene maps across multiple releases."""
    idmap = CIdMap()
    for ver in ENSEMBL_VERSIONS:
        MergeMap(idmap, ReadIdMapSimple(f"Biomart_map_Ensembl{ver}.txt", 1, 2, offset=1))
    return idmap


def run_cleanup() -> None:
    """Perform the HGNC alias cleanup and write the corrected output file."""
    idmap_ensp_hgnc = build_ensp_hgnc_map()
    idmap_ensp_ensg = build_ensp_ensg_map()

    idmap_hgnc_ensg_hgnc = ReadIdMapSimple("idmap_HGNC.txt", 8, 2, offset=1)
    MergeMap(idmap_hgnc_ensg_hgnc, ReadIdMapSimple("idmap_HGNC.txt", 10, 2, offset=1))

    approved_symbol_names = set(
        ReadIdMapSimple("idmap_HGNC.txt", 2, 1, offset=1).keys()
    )

    idmap_final: CIdMap = CIdMap()
    idmap_final_ambiguous: CIdMap = CIdMap()

    for ensp, hgncs in idmap_ensp_hgnc.iteritems():  # noqa: B007
        if len(hgncs) > 1:
            if idmap_ensp_ensg.has_key(ensp):  # noqa: W601
                hgncs_from_hgnc: set[str] = set()
                ensgs = idmap_ensp_ensg.GetIds(ensp)
                for ensg in ensgs:
                    if idmap_hgnc_ensg_hgnc.has_key(ensg):  # noqa: W601
                        hgncs_from_hgnc.update(idmap_hgnc_ensg_hgnc.GetIds(ensg))

                hgncs_both = hgncs.intersection(hgncs_from_hgnc)
                if len(hgncs_both) == 0:
                    if ensp in MANUAL_MAP:
                        idmap_final[ensp] = {MANUAL_MAP[ensp]}
                    elif len(approved_symbol_names.intersection(hgncs)) == 1:
                        idmap_final[ensp] = approved_symbol_names.intersection(hgncs)
                    else:
                        print(f"No mapping in HGNC: {ensp} {ensgs} {hgncs}")  # noqa: T201
                        idmap_final_ambiguous[ensp] = hgncs
                elif len(hgncs_both) > 1:
                    if ensp in MANUAL_MAP:
                        idmap_final[ensp] = {MANUAL_MAP[ensp]}
                    else:
                        print(f"{ensp} {hgncs_both}")  # noqa: T201
                        idmap_final_ambiguous[ensp] = hgncs_both
                else:
                    idmap_final[ensp] = hgncs_both
            else:
                raise RuntimeError(f"No Ensg mapped for {ensp}")
        else:
            idmap_final[ensp] = hgncs

    out_path = "9606.alias_best.v9.0.17032014_2.tsv"
    with open(out_path, "w") as fo:
        for ensp in idmap_final.keys():
            fo.write(
                f"9606\t{ensp}\t{idmap_final.GetFirstId(ensp)}"
                "\tHGNC symbol (Biomart, HGNC, Uniprot)\n"
            )

        with open("../data/9606.alias_best.v9.0.tsv") as f:
            for line in f:
                org, ensp, name, source = line.split()
                if ensp not in idmap_final:
                    if source in {
                        "BioMart_HUGO",
                        "Ensembl_HGNC_curated_gene",
                        "Ensembl_HGNC",
                    }:
                        fo.write(line)
                    elif ensp in idmap_final_ambiguous:
                        fo.write(
                            f"9606\t{ensp}\t{idmap_final_ambiguous.GetFirstId(ensp)}"
                            "\tHGNC symbol (Biomart, HGNC, Uniprot)\n"
                        )
                    else:
                        fo.write(line)


if __name__ == "__main__":
    run_cleanup()
