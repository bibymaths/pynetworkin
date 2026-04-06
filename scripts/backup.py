#!/usr/bin/env python3

"""
NetworKIN(tm), (C) 2005,2006,2007,2013.
Drs Rune Linding, Lars Juhl Jensen, Heiko Horn & Jinho Kim

Usage: ./networkin.py Organism FastaFile SitesFile

If no sites file is given NetworKIN will predict on all T/S/Y residues
in the given sequences.
"""

import sys
import os
import subprocess
import re
import tempfile
import glob
import platform
import time
import gzip
import random
from optparse import OptionParser
from pathlib import Path

import numpy as np
import pandas as pd

# Custom/Local Imports
from logger import logger
from likelihood import ReadConversionTableBin, ConvertScore2L
from motif_scoring import score_sequences
from inputs.phosphosites import fetch_phosphosite
from inputs.string_network import fetch_string_network
from recovery import recover_false_negatives
from output import write_output

# Weighting parameter, 0=only motif, 1=only STRING
ALPHAS = {"9606": 0.85, "4932": 0.65}
dSpeciesName = {"9606": "human", "4932": "yeast"}
dPenalty = {
    "9606": {"hub penalty": 100, "length penalty": 800},
    "4932": {"hub penalty": 170, "length penalty": 1000}
}

NETWORKIN_SITE_FILE = 1
PROTEOME_DISCOVERER_SITE_FILE = 2
MAX_QUANT_DIRECT_OUTPUT_FILE = 3
LEGACY_SITE_FILE = 4
MS_MCMC_FILE = 5

global options


################################################################################
#                                                                              #
#                             Core Functions                                   #
#                                                                              #
################################################################################

class CSheet(list):
    pass


def myPopen(cmd):
    if platform.system() == 'Windows':
        logger.info('WINDOWS')
        # Wrap the command to run in a Unix-like shell using WSL
        command = f'wsl bash -c "{cmd}"'
    else:
        command = cmd
    try:
        pipe = subprocess.Popen(command, shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = pipe.communicate()
        # Decode the byte strings to utf-8
        stdout = stdout.decode('utf-8')
        stderr = stderr.decode('utf-8')
        if pipe.returncode != 0:
            # If the command failed, print the stderr and raise an exception
            error_message = f"ERROR executing: {repr(command)}\n{stderr}"
            logger.error(error_message)
            raise subprocess.CalledProcessError(pipe.returncode, command, output=error_message)
        else:
            # If the command succeeded, return the stdout
            return stdout
    except subprocess.CalledProcessError as e:
        logger.error(f"Command '{e.cmd}' returned non-zero exit status {e.returncode}")
        return e.output
    except Exception as e:
        error_message = f"ERROR executing: {repr(command)}\n{str(e)}"
        logger.error(error_message)
        return error_message


def readFasta(fastafile):
    id_seq = {}
    aminoacids = re.compile('[^ACDEFGHIKLMNPQRSTVWYXB]')
    data = fastafile.readlines()
    fastafile.close()
    seq = ''
    for line in data:
        if line[0] != ';':
            if line[0] == '>':
                line = line.strip()
                id = line[1:].split('_', 1)[0]
                seq = ''
            else:
                seq += aminoacids.sub('', line)
                if len(seq) > 0:
                    id_seq[id] = seq
        else:
            continue
    return id_seq


def CheckInputType(sitesfile):
    with open(sitesfile, 'r') as f:
        line = f.readline()

    if '\t' in line:
        tokens = line.split('\t')
    elif ',' in line:
        tokens = line.split(',')
    elif ' ' in line:
        tokens = line.split(' ')

    if len(tokens) == 3:
        return NETWORKIN_SITE_FILE
    elif len(tokens) == 2:
        return PROTEOME_DISCOVERER_SITE_FILE
    elif tokens[0] == "Proteins" and tokens[4] == "Leading":
        return MAX_QUANT_DIRECT_OUTPUT_FILE
    elif tokens[1] == "phospho":
        return RUNES_SITE_FILE
    elif sitesfile.startswith('MS'):
        return MS_MCMC_FILE
    else:
        logger.error("Unknown format of site file")
        sys.exit()


def readPhosphoSites(sitesfile):
    id_pos_res = {}
    try:
        with open(sitesfile, 'r') as f:
            data = f.readlines()
            for line in data:
                tokens = line.split('\t')
                id = tokens[0]
                try:
                    pos = int(tokens[1])
                except Exception as e:
                    logger.error(line.strip())
                    raise e
                try:
                    res = tokens[2][:-1]
                except:
                    res = ""

                if id in id_pos_res:
                    id_pos_res[id][pos] = res
                else:
                    id_pos_res[id] = {pos: res}
    except IOError:
        logger.error(f"Could not open site file: {sitesfile}")
        sys.exit()

    return id_pos_res


def readPhosphoSitesProteomeDiscoverer(fastafile, sitesfile):
    # store all fasta sequences
    fasta = open(fastafile).readlines()
    fastadict = {}
    for line in fasta:
        if line[0] == ">":
            ensp = line.split()[0].strip(">")[0:15]
            fastadict[ensp] = ""
        else:
            seq = line.strip()
            fastadict[ensp] = seq

    # store all peptides in dictionary, with parent protein as key
    peptides = open(sitesfile).readlines()
    peptidedict = {}
    for line in peptides[1:]:
        line = line.strip()
        if '\t' in line:
            tokens = line.split('\t')
        elif ',' in line:
            tokens = line.split(',')
        protID = tokens[0][1:-1]
        peptide = tokens[1][1:-1]

        if protID in peptidedict.keys():
            if peptide not in peptidedict[protID]:
                peptidedict[protID][peptide] = ""
        else:
            peptidedict[protID] = {peptide: ""}

    # now map peptide onto full length sequence, then get absolute phosphosite locations
    id_pos_res = {}
    c_sequence_not_found = 0
    for protID in peptidedict.keys():
        for peptide in peptidedict[protID]:
            UPPERpeptide = peptide.upper()
            sequence = fastadict[protID]

            try:
                peptideindex = sequence.index(UPPERpeptide)
            except:
                c_sequence_not_found += 1

            x = 0
            for letter in peptide:
                if letter.islower():
                    if letter in ["s", "t", "y"]:
                        phoslocation = peptideindex + x + 1
                        phosresidue = sequence[phoslocation - 1]

                        if protID in id_pos_res.keys():
                            id_pos_res[protID][phoslocation] = phosresidue
                        else:
                            id_pos_res[protID] = {phoslocation: phosresidue}
                x += 1

    logger.info(f"sequences not found in fasta: {c_sequence_not_found}")
    return id_pos_res


def ReadSheet(fname, offset=0):
    l = CSheet()
    f = open(fname, 'r')

    for i in range(offset):
        f.readline()

    columns = f.readline().strip().split('\t')
    l.columns = columns

    for line in f.readlines():
        instance = {}
        l.append(instance)
        fields = line.strip().split('\t')

        if len(fields) > len(columns):
            logger.error("Error in data file")
            logger.info(str(columns))
            logger.info(str(fields))
            raise

        for i in range(len(columns)):
            try:
                instance[columns[i]] = fields[i]
            except IndexError:
                instance[columns[i]] = ''

    f.close()
    return l


def readPhosphoSitesMaxQuant(fname, only_leading=False):
    id_pos_res = {}
    phosphosites = ReadSheet(fname)

    for site in phosphosites:
        Ids = site["Proteins"].split(';')
        positions = [int(x) for x in site["Positions within proteins"].split(';')]
        aa = site["Amino acid"]
        leading_protein_ids = site["Leading proteins"].split(';')

        for i in range(len(Ids)):
            Id = Ids[i]
            pos = positions[i]

            if only_leading and Id not in leading_protein_ids:
                continue

            if Id in id_pos_res:
                id_pos_res[Id][pos] = aa
            else:
                id_pos_res[Id] = {pos: aa}

    return id_pos_res


def readRunessitesfile(sitesfile):
    id_pos_res = {}
    try:
        with open(sitesfile, 'r') as f:
            data = f.readlines()
            for line in data:
                tokens = line.split(' ')
                id = tokens[3]
                res_pos = tokens[2]
                pos = int(res_pos[2:])
                res = res_pos[0]
                if id in id_pos_res:
                    id_pos_res[id][pos] = res
                else:
                    id_pos_res[id] = {pos: res}
    except IOError:
        logger.error(f"Could not open site file: {sitesfile}")
        sys.exit()

    return id_pos_res


def readMCMCssitesfile(sitesfile):
    id_pos_res = {}
    try:
        with open(sitesfile, 'r') as f:
            data = f.readlines()
            for line in data:
                tokens = line.split(' ')
                id = tokens[0]
                # Note: `res_pos` reference follows original script logic
                if (res_pos != '') or (res_pos != 'M_'):
                    res_pos = tokens[2]
                    pos = int(res_pos[2:])
                    res = res_pos[0]
                    if id in id_pos_res:
                        id_pos_res[id][pos] = res
                    else:
                        id_pos_res[id] = {pos: res}
    except IOError:
        logger.error(f"Could not open site file: {sitesfile}")
        sys.exit()

    return id_pos_res


def readAliasFiles(organism, datadir, map_group_to_domain):
    alias_hash = {}
    desc_hash = {}
    name_hash = {}
    names_list = []

    for tree in map_group_to_domain:
        for pred in map_group_to_domain[tree]:
            for name in map_group_to_domain[tree][pred]:
                if name not in names_list:
                    names_list.append(name)

    try:
        desc_db = myPopen(f'gzip -cd {datadir}/{organism}.text_best.v9.0.tsv.gz')
        desc_db = desc_db.split('\n')
        for line in desc_db:
            (taxID, seqID, desc) = line.split('\t')[:3]
            desc_hash[seqID] = desc
    except:
        logger.info(f"No descriptions available for organism: '{organism}'")

    try:
        name_db = myPopen(f'gzip -cd {datadir}/{organism}.protein.aliases.v12.0.txt.gz')
        name_db = name_db.split('\n')
        name_hash = {}
        for line in name_db:
            (seqID, alias, source) = line.split('\t')
            key = seqID[5:]
            if key in name_hash.keys():
                if alias in name_hash[key]:
                    continue
                else:
                    name_hash[key].append(alias)
            else:
                name_hash[key] = [alias]
    except:
        logger.info(f"No names available for organism: '{organism}'")

    for id in name_hash.keys():
        flag = False
        for name in name_hash[id]:
            if name in names_list:
                alias_hash[id] = name
                flag = True
        if not flag:
            alias_hash[id] = id

    return alias_hash, desc_hash, name_hash


def ReadLines(fname):
    with open(fname) as f:
        lines = f.readlines()
    return lines


def WriteString(fname, s):
    with open(fname, 'w') as f:
        f.write(s)


def mapPeptides2STRING(blastDir, organism, fastafilename, id_pos_res, id_seq, number_of_processes, datadir, fast=False,
                       leave_intermediates=False):
    logger.info("Mapping using blast")
    incoming2string = {}
    string2incoming = {}

    blast_tmpfile = tempfile.NamedTemporaryFile(mode='w', delete=False)
    if id_pos_res == {}:
        for id in id_seq:
            blast_tmpfile.write('>' + id + '\n' + id_seq[id] + '\n')
    else:
        for id in id_pos_res:
            try:
                blast_tmpfile.write('>' + id + '\n' + id_seq[id] + '\n')
            except:
                logger.warning(f"No sequence available for '{id}'")
    blast_tmpfile.flush()

    blastDB = os.path.join(datadir, f"{organism}.protein.sequences.v12.0.fa").replace(' ', '\\ ')

    if not os.path.isfile(blastDB + '.pin'):
        command = f"{blastDir.rsplit('/', 1)[0]}/makeblastdb -in {blastDB} -out {blastDB} -parse_seqids -dbtype prot"
        logger.info(f"Looks like blast database is not initialized, trying to run:\n{command}")
        myPopen(command)

    command = f"{blastDir}blastp -num_threads {number_of_processes} -evalue 1e-10 -db {blastDB} -query {blast_tmpfile.name} -outfmt 6 | sort -k12gr"

    fn_blast_output = f"{fastafilename}.{organism}.blast.out"
    if fast and os.path.isfile(fn_blast_output):
        blast_out = ReadLines(fn_blast_output)
    else:
        blast_out = myPopen(command)
        if leave_intermediates:
            WriteString(fn_blast_output, "".join(blast_out))

    if isinstance(blast_out, str):
        blast_out = blast_out.split('\n')

    for line in blast_out:
        if line == '':
            continue
        tokens = line.split('\t')
        incoming = tokens[0]

        if incoming not in incoming2string:
            string = tokens[1].replace(f"{organism}.", "")
            identity = float(tokens[2])
            evalue = float(tokens[-2])

            if string in string2incoming:
                logger.warning(f"Best hit for {incoming} is not reciprocal")
            if identity < 90:
                logger.warning(f"Best hit for {incoming} has only {identity:.2f} %identity")
            if evalue > 1e-40:
                logger.warning(f"Best hit for {incoming} has high E-value {evalue:.2e}")

            if incoming in incoming2string:
                incoming2string[incoming][string] = True
            else:
                incoming2string[incoming] = {string: True}

            if string in string2incoming:
                string2incoming[string][incoming] = True
            else:
                string2incoming[string] = {incoming: True}

    return incoming2string, string2incoming


def mapRandom(id_seq):
    logger.info("Mapping random")
    incoming2string = {}
    string2incoming = {}

    stringIDs = []
    with open(f'{options.datadir}/{organism}.protein.sequences.fa', 'r') as file:
        data = file.readlines()

    for line in data:
        if line[0] == '>':
            name = line[1:-1]
            stringIDs.append(name)
            string2incoming[name] = {}

    max_idx = len(stringIDs) - 1

    for incoming in id_seq:
        if incoming not in incoming2string:
            rand_int = random.randint(0, max_idx)
            string = stringIDs[rand_int]
            incoming2string[incoming] = {string: True}
            if string in string2incoming:
                string2incoming[string][incoming] = True
            else:
                string2incoming[string] = {incoming: True}

    return incoming2string, string2incoming


def mapOne2one(id_seq):
    logger.info("Mapping one2one")
    incoming2string = {}
    string2incoming = {}
    for incoming in id_seq:
        incoming2string[incoming] = {incoming: True}
        string2incoming[incoming] = {incoming: True}
    return incoming2string, string2incoming


def mapFromFile(filename):
    logger.info("Mapping using external mapping file")
    incoming2string = {}
    string2incoming = {}

    command = f"cat {filename}"
    try:
        mappingFile = myPopen(command)
        mappingFile = mappingFile.split('\n')
    except:
        logger.error(f"Going to sleep, crashed with '{command}'")
        time.sleep(3600)
        return incoming2string, string2incoming

    for line in mappingFile:
        if re.match('^#', line) or line == '':
            continue

        line = line.strip()
        tokens = line.split('\t')

        incoming = tokens[0]
        string = tokens[1]

        if incoming in incoming2string:
            incoming2string[incoming][string] = True
        else:
            incoming2string[incoming] = {string: True}

        if string in string2incoming:
            string2incoming[string][incoming] = True
        else:
            string2incoming[string] = {incoming: True}

    return incoming2string, string2incoming


def loadSTRINGdata(string2incoming, datadir, alias_hash, number_of_processes):
    fn_bestpath = os.path.join(datadir, "string_data", f"{organism}.bestpath_0340_0950.v9.tsv.gz")

    if not os.path.isfile(fn_bestpath):
        logger.warning(f"Best path file does not exist: {fn_bestpath}")

    tree_pred_string_data = {}

    f = gzip.open(fn_bestpath)
    f = f.read().decode('utf-8')
    f = f.split('\n')

    c = 0
    for line in f:
        c += 1
        if line == '':
            continue

        tokens = line.split('\t')
        if len(tokens) == 8:
            (tree, group, name, string1, string2, stringscore, stringscore_indirect, path) = tokens
        elif len(tokens) == 7:
            (tree, group, name, string1, string2, stringscore, stringscore_indirect) = tokens
        elif len(tokens) == 6:
            (tree, group, name, string1, string2, stringscore) = tokens
            stringscore = float(stringscore) / 1000
            string1 = string1.replace(f'{organism}.', '')
            string2 = string2.replace(f'{organism}.', '')
        elif len(tokens) == 3:
            (string1, string2, stringscore) = tokens
            stringscore = float(stringscore) / 1000
            string1 = string1.replace(f'{organism}.', '')
            string2 = string2.replace(f'{organism}.', '')
            name = alias_hash.get(string1, 'notdef')

        if string2 in string2incoming:
            if string2 in tree_pred_string_data:
                tree_pred_string_data[string2][string1] = {"_name": name}
            else:
                tree_pred_string_data[string2] = {string1: {"_name": name}}

            if options.path == "direct":
                tree_pred_string_data[string2][string1]["_score"] = float(stringscore)
            elif options.path == "indirect":
                tree_pred_string_data[string2][string1]["_score"] = float(stringscore_indirect)
            else:
                raise ValueError("Path information should be either direct or indirect.")

            if len(tokens) == 9:
                tree_pred_string_data[string2][string1]["_path"] = path
            elif len(tokens) in [6, 3]:
                tree_pred_string_data[string2][string1]["_path"] = 'notdef'

    logger.info(f"network edges: {c}")
    return tree_pred_string_data


def InsertValueIntoMultiLevelDict(d, keys, value):
    for i in range(len(keys) - 1):
        if not (keys[i] in d.keys()):
            d[keys[i]] = {}
        d = d[keys[i]]

    if not (keys[-1] in d.keys()):
        d[keys[-1]] = []
    d[keys[-1]].append(value)


def ReadGroup2DomainMap(path_group2domain_map):
    map_group2domain = {}
    use_curated = True

    if use_curated:
        with open("data/group_human_protein_name_map_curated.tsv", "r") as f:
            for line in f.readlines():
                tokens = line.split()
                name = tokens[3]
                InsertValueIntoMultiLevelDict(map_group2domain, tokens[:2], name)
    else:
        with open(path_group2domain_map, "r") as f:
            for line in f.readlines():
                tokens = line.split()
                name = tokens[2]
                InsertValueIntoMultiLevelDict(map_group2domain, tokens[:2], name)

    return map_group2domain


def SetValueIntoMultiLevelDict(d, keys, value):
    for i in range(len(keys) - 1):
        if keys[i] not in d:
            d[keys[i]] = {}
        d = d[keys[i]]

    if (keys[-1] in d) and type(d[keys[-1]]) != type(value):
        logger.warning("Caution: multi-dict already has value and try to assign a value of different type")

    if (keys[-1] in d):
        if d[keys[-1]] != value:
            logger.info(f"This operation replaces a value ({' '.join(map(str, keys))})")

    d[keys[-1]] = value


def filter_and_rank_predictions(predictions, min_networkin=2.0, min_motif=0.05, top_k=5):
    df = pd.DataFrame(predictions)
    if df.empty:
        return predictions

    df = df[
        (df["NetworKIN score"] > min_networkin) &
        (df["Motif probability"] > min_motif)
        ].copy()

    df = df.sort_values(
        ["Name", "Position", "NetworKIN score"],
        ascending=[True, True, False],
    )

    df = df.groupby(["Name", "Position"], as_index=False).head(top_k)
    return df.to_dict(orient="records")


def printResult(id_pos_tree_pred, tree_pred_string_data, incoming2string, string_alias, name_hash, string_desc,
                organism, mode, dir_likelihood_conversion_tbl, map_group2domain, res_dir, fn_fasta):
    ALPHA = ALPHAS[organism]
    species = dSpeciesName[organism]

    csv_filename = os.path.join(res_dir, f"{Path(fn_fasta).name}.result.tsv")
    predictions = []

    dLRConvTbl = {}
    for fname in glob.glob(os.path.join(dir_likelihood_conversion_tbl, "conversion_tbl_*_smooth*")):
        match = re.findall(r"conversion_tbl_([a-z]+)_smooth_([a-z]+)_([A-Z0-9]+)_([a-zA-Z0-9_/-]+)",
                           os.path.basename(os.path.splitext(fname)[0]))
        if not match:
            continue

        score_src, species_of_conversion_table, tree, player_name = match[0]

        if score_src == "string":
            score_src = "string"
        else:
            score_src = "motif"

        if species_of_conversion_table != species:
            continue

        conversion_tbl = ReadConversionTableBin(fname)
        SetValueIntoMultiLevelDict(dLRConvTbl, [species_of_conversion_table, tree, player_name, score_src],
                                   conversion_tbl)

        if options.verbose:
            logger.info(f"Conversion table {species_of_conversion_table} {tree} {player_name} {score_src}")

    c_np = 0
    c_map = 0
    c_notmap = 0
    c_string = 0
    c_stringscore = 0
    c_if = 0
    c_else = 0
    c_res = 0
    pred_not_mapped = []
    name_not_mapped = []

    for id in id_pos_tree_pred:
        c_np += 1
        if id in incoming2string:
            c_map += 1
            for pos in id_pos_tree_pred[id]:
                for tree in id_pos_tree_pred[id][pos]:
                    score_results = {}
                    for pred in id_pos_tree_pred[id][pos][tree]:
                        for string1 in incoming2string[id]:
                            bestName1 = string_alias.get(string1, 'notdef')
                            desc1 = string_desc.get(string1, 'notdef')

                            if string1 in tree_pred_string_data:
                                c_string += 1
                                (res, peptide, motifScore) = id_pos_tree_pred[id][pos][tree][pred]
                                for string2 in tree_pred_string_data[string1]:
                                    bestName2 = string_alias.get(string2, 'notdef')
                                    desc2 = string_desc.get(string2, 'notdef')

                                    stringScore = tree_pred_string_data[string1][string2]["_score"]
                                    c_stringscore += 1
                                    path = 'notdef'
                                    name = tree_pred_string_data[string1][string2]["_name"]

                                    try:
                                        map_names = map_group2domain[tree][pred]
                                    except:
                                        pred_not_mapped.append(pred)

                                    if (tree not in map_group2domain.keys()) or (
                                            pred not in map_group2domain[tree].keys()) or (
                                            name not in map_group2domain[tree][pred]):
                                        if (tree in map_group2domain) and (pred in map_group2domain[tree]) and (
                                                name not in map_group2domain[tree][pred]):
                                            name_not_mapped.append(name)

                                        if options.string_for_uncovered:
                                            if species == "human":
                                                if tree in ["1433", "BRCT", "WW", "PTB", "WD40", "FHA"]:
                                                    conversion_tbl_string = dLRConvTbl[species]["SH2"]["general"][
                                                        "string"]
                                                else:
                                                    conversion_tbl_string = dLRConvTbl[species][tree]["general"][
                                                        "string"]
                                            elif species == "yeast":
                                                conversion_tbl_string = dLRConvTbl[species][tree]["general"]["string"]
                                            else:
                                                raise ValueError("This species is not supported")

                                            likelihood_motif = 1
                                            likelihood_string = ConvertScore2L(stringScore, conversion_tbl_string)
                                            unified_likelihood = likelihood_motif * likelihood_string
                                            networkinScore = unified_likelihood

                                            row = None
                                            if networkinScore >= 0.02:
                                                c_if += 1
                                                row = {
                                                    "Name": id,
                                                    "Position": res + str(pos),
                                                    "Tree": tree,
                                                    "Motif Group": pred,
                                                    "Kinase/Phosphatase/Phospho-binding domain": name,
                                                    "NetworKIN score": float(format(networkinScore, ".4f")),
                                                    "Motif probability": float(format(motifScore, ".4f")),
                                                    "STRING score": float(format(stringScore, ".4f")),
                                                    "Target STRING ID": string1,
                                                    "Kinase STRING ID": string2,
                                                    "Target Name": bestName1,
                                                    "Kinase Name": bestName2,
                                                    "Target description": desc1,
                                                    "Kinase description": desc2,
                                                    "Peptide sequence window": peptide,
                                                    "Intermediate nodes": path,
                                                    "recovered": False,
                                                    "recovery_method": "",
                                                }
                                        else:
                                            continue
                                    else:
                                        c_else += 1
                                        if species == "human":
                                            if tree in ["1433", "BRCT", "WW", "PTB", "WD40", "FHA"]:
                                                conversion_tbl_motif = dLRConvTbl[species]["SH2"]["general"]["motif"]
                                                conversion_tbl_string = dLRConvTbl[species]["SH2"]["general"]["string"]
                                            else:
                                                if name in dLRConvTbl[species][tree]:
                                                    conversion_tbl_motif = dLRConvTbl[species][tree][name]["motif"]
                                                    conversion_tbl_string = dLRConvTbl[species][tree][name]["string"]
                                                else:
                                                    conversion_tbl_motif = dLRConvTbl[species][tree]["general"]["motif"]
                                                    conversion_tbl_string = dLRConvTbl[species][tree]["general"][
                                                        "string"]
                                        elif species == "yeast":
                                            if name in dLRConvTbl[species][tree]:
                                                conversion_tbl_motif = dLRConvTbl[species][tree][name]["motif"]
                                                conversion_tbl_string = dLRConvTbl[species][tree][name]["string"]
                                            else:
                                                conversion_tbl_motif = dLRConvTbl[species][tree]["general"]["motif"]
                                                conversion_tbl_string = dLRConvTbl[species][tree]["general"]["string"]
                                        else:
                                            raise ValueError("This species is not supported")

                                        likelihood_motif = ConvertScore2L(motifScore, conversion_tbl_motif)
                                        likelihood_string = ConvertScore2L(stringScore, conversion_tbl_string)
                                        unified_likelihood = likelihood_motif * likelihood_string
                                        networkinScore = unified_likelihood

                                        c_res += 1
                                        row = {
                                            "Name": id,
                                            "Position": res + str(pos),
                                            "Tree": tree,
                                            "Motif Group": pred,
                                            "Kinase/Phosphatase/Phospho-binding domain": name,
                                            "NetworKIN score": float(format(networkinScore, ".4f")),
                                            "Motif probability": float(format(motifScore, ".4f")),
                                            "STRING score": float(format(stringScore, ".4f")),
                                            "Target STRING ID": string1,
                                            "Kinase STRING ID": string2,
                                            "Target Name": bestName1,
                                            "Kinase Name": bestName2,
                                            "Target description": desc1,
                                            "Kinase description": desc2,
                                            "Peptide sequence window": peptide,
                                            "Intermediate nodes": path,
                                            "recovered": False,
                                            "recovery_method": "",
                                        }

                                    if row is not None:
                                        predictions.append(row)

                                    if networkinScore not in score_results:
                                        score_results[networkinScore] = []
                                    score_results[networkinScore].append(row)
        else:
            c_notmap += 1

    logger.warning(f"names not mapped: {len(np.unique(name_not_mapped))}")

    # --- False-negative recovery ---
    motif_score_dict = {}
    for prot_id, pos_data in id_pos_tree_pred.items():
        if prot_id not in incoming2string:
            continue
        for string1 in incoming2string[prot_id]:
            if string1 not in tree_pred_string_data:
                continue
            for pos, tree_data in pos_data.items():
                for tree, pred_data in tree_data.items():
                    for pred_name, (res, peptide, ms) in pred_data.items():
                        for string2 in tree_pred_string_data[string1]:
                            key = (string2, string1)
                            motif_score_dict[key] = max(ms, motif_score_dict.get(key, -1.0))

    all_string_pairs = []
    protein_ids = set()
    for sub, kins in tree_pred_string_data.items():
        protein_ids.add(sub)
        for kin in kins:
            protein_ids.add(kin)
            all_string_pairs.append((kin, sub))

    protein_list = sorted(protein_ids)
    node_index = {pid: i for i, pid in enumerate(protein_list)}
    n = len(protein_list)

    dist_matrix = np.full((n, n), np.inf, dtype=np.float32)
    np.fill_diagonal(dist_matrix, 0.0)
    for sub, kins in tree_pred_string_data.items():
        si = node_index[sub]
        for kin, data in kins.items():
            ki = node_index[kin]
            score = data.get("_score", 0.0)
            if score > 0:
                d = float(1.0 / score - 1.0)
                if d < dist_matrix[ki, si]:
                    dist_matrix[ki, si] = d

    recovered_preds = recover_false_negatives(
        candidates=all_string_pairs,
        dist_matrix=dist_matrix,
        node_index=node_index,
        motif_scores=motif_score_dict,
    )

    for r in recovered_preds:
        sub_id = r["substrate_uniprot"]
        kin_id = r["kinase_id"]
        sub_name = string_alias.get(sub_id, sub_id)
        kin_name = string_alias.get(kin_id, kin_id)
        sub_desc = string_desc.get(sub_id, 'notdef')
        kin_desc = string_desc.get(kin_id, 'notdef')

        predictions.append({
            "Name": sub_name,
            "Position": "",
            "Tree": "",
            "Motif Group": "",
            "Kinase/Phosphatase/Phospho-binding domain": kin_name,
            "NetworKIN score": float(format(r["networkin_score"], ".4f")),
            "Motif probability": float(format(r["motif_score"], ".4f")),
            "STRING score": float(format(r["context_score"], ".4f")),
            "Target STRING ID": sub_id,
            "Kinase STRING ID": kin_id,
            "Target Name": sub_name,
            "Kinase Name": kin_name,
            "Target description": sub_desc,
            "Kinase description": kin_desc,
            "Peptide sequence window": "",
            "Intermediate nodes": "notdef",
            "recovered": r["recovered"],
            "recovery_method": r["recovery_method"],
        })

    predictions = filter_and_rank_predictions(predictions)
    write_output(predictions, csv_filename)

    return c_np, c_map, c_notmap, c_string, c_stringscore, c_if, c_else, c_res


def Main():
    res_dir = 'results'
    if not os.path.exists(res_dir):
        os.makedirs(res_dir)

    tmpdir = 'tmp'
    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)

    logger.info("Reading fasta input file")
    id_seq = readFasta(fastafile)
    if options.verbose:
        logger.info(f"{len(id_seq.keys())} sequences loaded")

    if sitesfile:
        logger.info("Reading phosphosite file")
        input_type = CheckInputType(sitesfile)
        if input_type == NETWORKIN_SITE_FILE:
            id_pos_res = readPhosphoSites(sitesfile)
        elif input_type == PROTEOME_DISCOVERER_SITE_FILE:
            id_pos_res = readPhosphoSitesProteomeDiscoverer(fn_fasta, sitesfile)
        elif input_type == MAX_QUANT_DIRECT_OUTPUT_FILE:
            id_pos_res = readPhosphoSitesMaxQuant(sitesfile)
        elif input_type == LEGACY_SITE_FILE:
            id_pos_res = readRunessitesfile(sitesfile)
        elif input_type == MS_MCMC_FILE:
            id_pos_res = readMCMCssitesfile(sitesfile)
    else:
        id_pos_res = {}

    logger.info("Fetching live phosphorylation site reference data")
    phosphosite_df = fetch_phosphosite(refresh=options.refresh)  # noqa: F841

    if organism == "9606":
        path_group2domain_map = os.path.join(options.datadir, "group_human_protein_name_map.tsv")
    elif organism == "4932":
        path_group2domain_map = os.path.join(options.datadir, "group_yeast_KIN.tsv")

    logger.info("Loading aliases and descriptions")
    map_group2domain = ReadGroup2DomainMap(path_group2domain_map)
    (string_alias, string_desc, name_hash) = readAliasFiles(args[0], options.datadir, map_group2domain)

    incoming2string, string2incoming = mapPeptides2STRING(
        blastDir, organism, fastafile.name, id_pos_res, id_seq, options.threads, options.datadir
    )

    logger.info("Fetching live STRING network data")
    all_uniprot_ids = list(id_seq.keys())
    string_df = fetch_string_network(proteins=all_uniprot_ids, refresh=options.refresh)  # noqa: F841

    logger.info("Loading STRING network")
    tree_pred_string_data = loadSTRINGdata(string2incoming, options.datadir, string_alias, options.threads)

    logger.info("Running motif scorer")
    id_pos_tree_pred = score_sequences(id_seq, id_pos_res)

    logger.info("Writing results")
    logger.info(
        "#Name\tPosition\tTree\tMotif Group\tKinase/Phosphatase/Phospho-binding domain\tNetworKIN score\tMotif probability\tSTRING score\tTarget STRING ID\tKinase/Phosphatase/Phospho-binding domain STRING ID\tTarget description\tKinase/Phosphatase/Phospho-binding domain description\tTarget Name\tKinase/Phosphatase/Phospho-binding domain Name\tPeptide sequence window\tIntermediate nodes")

    if options.path == "direct":
        dir_likelihood_conversion_tbl = os.path.join(options.datadir, "likelihood_conversion_table_direct")
    elif options.path == "indirect":
        dir_likelihood_conversion_tbl = os.path.join(options.datadir, "likelihood_conversion_table_indirect")
    else:
        raise ValueError("Path information should be either direct or indirect.")

    c_np, c_map, c_notmap, c_string, c_stringscore, c_if, c_else, c_res = printResult(
        id_pos_tree_pred, tree_pred_string_data, incoming2string, string_alias, name_hash, string_desc,
        args[0], options.mode, dir_likelihood_conversion_tbl, map_group2domain, res_dir, fn_fasta
    )

    summary_str = (
        f"c_np = {c_np}\n"
        f"c_map = {c_map}\n"
        f"c_notmap = {c_notmap}\n"
        f"c_string = {c_string}\n"
        f"c_stringscore = {c_stringscore}\n"
        f"c_if = {c_if}\n"
        f"c_else = {c_else}\n"
        f"c_res = {c_res}"
    )
    logger.info(summary_str)

    return


if __name__ == '__main__':
    try:
        blastDir = os.environ['BLAST_PATH']
    except KeyError:
        blastDir = ""

    usage = "usage: %prog [options] organism FASTA-file [sites-file]"
    parser = OptionParser(usage=usage, version="%prog 3.0")
    parser.add_option("-b", "--blast", dest="blast", default=blastDir,
                      help="set the directory for the BLAST binaries (formatdb and blastall), overwrites the 'BLAST_PATH' environmental variable. [ENV: %default]")
    parser.add_option("-m", "--mode", dest="mode", default=False,
                      help="if set to 'network', gives only one best scoring result for each site. In case of multiple candidate kinases with the same core, the selection hapens randomly. [default: %default]")
    parser.add_option("-p", "--path", dest="path", default="direct",
                      help="NetworKIN uses both direct and indirect paths. Otherwise, it uses only indirect paths. [default: %default]")
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true",
                      help="print out everything [default: %default]")
    parser.add_option("-f", "--fast", dest="fast", default=False, action="store_true",
                      help="Speed up by using the intermediate files of previous run [default: %default]")
    parser.add_option("-l", "--leave", dest="leave", default=False, action="store_true",
                      help="leave intermediate files [default: %default]")
    parser.add_option("-u", "--uncovered", dest="string_for_uncovered", default=False, action="store_true",
                      help="Use STRING likelihood for uncovered Kinases [default: %default]")
    parser.add_option("-t", "--threads", dest="threads", default=1, type="int",
                      help="number of available threads/CPUs. Also leads to less memory usage, as result files are read sequentially [default: %default]")
    parser.add_option("--nt", dest="active_threads", default=2, type="int",
                      help="number of active threads at a time")
    parser.add_option("-c", "--compress", dest="compress", default=True,
                      help="compress temporary result files, saves discspace [default: %default]")
    parser.add_option("-d", "--data", dest="datadir",
                      default=os.path.join(os.path.split(os.path.realpath(sys.argv[0]))[0], 'data'),
                      help="location for the additional files like the pre-computed STRING network, STRING sequence database etc. [default: %default]")
    parser.add_option("--refresh", dest="refresh", default=False, action="store_true",
                      help="force re-download of live data (PhosphoSitePlus, STRING) even if a valid 7-day cache exists [default: %default]")

    (options, args) = parser.parse_args()

    if options.active_threads < 2:
        parser.error("Number of active thread (--nt) is less than 2")

    try:
        organism = args[0]
    except IndexError:
        parser.error("Organism not defined!")

    try:
        fn_fasta = args[1]
        fn_blast_output = f"{fn_fasta}.{organism}.blast.out"
        fastafile = open(fn_fasta, 'r')
    except IndexError:
        logger.error(str(args))
        parser.error("FASTA-file not defined!")
    except FileNotFoundError:
        parser.error(f"Cannot open FASTA file: {args[1]}")

    try:
        sitesfile = args[2]
    except IndexError:
        sitesfile = False

    if options.blast:
        blastDir = options.blast

    if options.verbose:
        logger.info(f"\nPredicting using parameters as follows:\nOrganism:\t{organism}\nFastaFile:\t{fn_fasta}")
        logger.info(f"Threads:\t{options.threads}\nCompress:\t{options.compress}")
        if options.string_for_uncovered:
            logger.info("Use STRING likelihood when kinases are not covered by the motif atlas.")
        if sitesfile:
            logger.info(f"Sitesfile:\t{sitesfile}")
        else:
            logger.info("No sites-file given, predicting on all S,T,Y residues.")
        if options.mode:
            logger.info(f"Mode:\t\t{options.mode}")
        logger.info(f"Blast dir: {blastDir}")

    Main()
