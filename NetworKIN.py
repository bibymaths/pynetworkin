#!/usr/bin/env python3

"""
NetworKIN(tm), (C) 2005,2006,2007,2013.
Drs Rune Linding, Lars Juhl Jensen, Heiko Horn & Jinho Kim

Usage: ./networkin.py Organism FastaFile SitesFile

If no sites file is given NetworKIN will predict on all T/S/Y residues
in the given sequences.
"""

### Coding strategy & TODO
# FIX: filtering when no sites given
# 2) add code for prediction of downstream recognition site (14-3-3 first)
# 3) add code for predicting downstream binding module (14-3-3 first)
# 4) finnish up output format and CLI options
# 5) add code for doing 'pathway' analysis
# 6) predict phenotypes
###

# Changelog
# 21.10.06: tested blast -a 8 option, no differences on proteome
# 21.10.06: Filtercode broken
#			-fixed
# 28.01.07: Testing autophosphorylation (self == 1 in update script)
# 30.07.07: Working on v1.5 milestone
# 07.05.08: Initiated 2.5 w. scaling factor
# 11.07.13: New scoring scheme (Bayesian)
# 03.02.14: Several minor bug fixes. blast sorting bug fix (k12nr -> k12gr)
#




import sys, os, subprocess, re, tempfile, random, operator
import threading
from optparse import OptionParser
import multiprocessing
import glob
import csv
from string import *
from likelihood import ReadConversionTableBin
from likelihood import ConvertScore2L
from motif_scoring import score_sequences
from inputs.phosphosites import fetch_phosphosite
from inputs.string_network import fetch_string_network
from recovery import recover_false_negatives
from output import write_output
import platform
from itertools import chain
import numpy as np

# debugging
import time

# Weighting parameter, 0=only motif, 1=only STRING
# estimated while benchmarking an then hardcoded here
# feel free to play with it, but it is your own responsibility
ALPHAS = {"9606": 0.85, "4932": 0.65}
dSpeciesName = {"9606": "human", "4932": "yeast"}
dPenalty = {"9606": {"hub penalty": 100, "length penalty": 800}, "4932": {"hub penalty": 170, "length penalty": 1000}}

NETWORKIN_SITE_FILE = 1
PROTEOME_DISCOVERER_SITE_FILE = 2
MAX_QUANT_DIRECT_OUTPUT_FILE = 3
RUNES_SITE_FILE = 4
MS_MCMC_FILE = 5

global options


# Temporary directory to store the files
# tempfile.tempdir= '/tmp'

# Location of NetworKIN input files (string network, alias file etc.)
# datadir = sys.argv[0].rsplit("/", 1)[0]+'/data'
# global DATADIR
# DATADIR = ""

# Number of threads used for motif scoring and BLAST
# setting it to a high number can also help in case NetworKIN uses too much memory
# -> the more files, the less memory usage
# global NUMBER_OF_PROCESSES
# NUMBER_OF_PROCESSES = "1";

# Should the analysis be limited to specific trees?
# just comment it out if you want to predict on all
# limitTrees = ['SH2']
# limitTrees = ['KIN', "SH2", "PTP"]
# limitTrees = ['KIN']

################################################################################
#                                                                              #
#                             the code starts here                             #
#                                                                              #
################################################################################

class CSheet(list):
    pass


'''
# Run system binary
def myPopen(cmd):
	try:
		pipe = subprocess.Popen(cmd, shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdout = pipe.stdout.readlines()
	except:
		sys.stderr.write('ERROR executing: '+`cmd`+'\n')
		sys.exit()
	else:
		return stdout 
'''


### myPopen:
def myPopen(cmd):
    if platform.system() == 'Windows':
        sys.stderr.write('WINDOWS\n')
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
            error_message = 'ERROR executing: ' + repr(command) + '\n' + stderr
            sys.stderr.write(error_message + '\n')
            raise subprocess.CalledProcessError(pipe.returncode, command, output=error_message)
        else:
            # If the command succeeded, return the stdout
            return stdout
    except subprocess.CalledProcessError as e:
        # Handle the subprocess.CalledProcessError exception
        print("Error: Command '{}' returned non-zero exit status {}".format(e.cmd, e.returncode))
        return e.output  # Return the error message
    except Exception as e:
        # Handle other exceptions
        error_message = 'ERROR executing: ' + repr(command) + '\n' + str(e) + '\n'
        sys.stderr.write(error_message)
        return error_message


# Read sequences from fasta file
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
                # id = line[1:].split(' ', 1)[0]
                id = line[1:].split('_', 1)[0]
                #				id = line[1:-1].split(' ', 1)[0]
                seq = ''
            else:
                seq += aminoacids.sub('', line)
                if len(seq) > 0:
                    id_seq[id] = seq
        else:
            continue
    return id_seq


def CheckInputType(sitesfile):
    f = open(sitesfile, 'r')
    line = f.readline()
    f.close()
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
        sys.stderr("Unknown format of site file")
        sys.exit()


# Read phosphorylation sites from tsv file
# id -> position -> residue
def readPhosphoSites(sitesfile):
    id_pos_res = {}
    f = open(sitesfile, 'r')
    if f:
        data = f.readlines()
        f.close()
        for line in data:
            tokens = line.split('\t')
            id = tokens[0]
            try:
                pos = int(tokens[1])
            except:
                sys.stderr.write(line)
                raise
            try:
                res = tokens[2][:-1]
            except:
                res = ""
            if id in id_pos_res:
                id_pos_res[id][pos] = res
            else:
                id_pos_res[id] = {pos: res}
    else:
        sys.stderr.write("Could not open site file: %s" % sitesfile)
        sys.exit()

    return id_pos_res


def readPhosphoSitesProteomeDiscoverer(fastafile, sitesfile):
    # This function takes a tab-separated list of ProteinID \t Peptide, in ProteomeDiscoverer format (i.e. phosphosites listed as
    # lowercase), maps the peptide to the full-length protein sequence and conducts the NetworKIN search.
    #

    # store all fasta sequences
    fasta = open(fastafile).readlines()
    fastadict = {}
    for line in fasta:
        if line[0] == ">":  # use this as each new entry in FASTA file starts with '>'
            ensp = line.split()[0].strip(">")[0:15]  # strip the '>' symbol from the line before copying value into dictionary
            # print ensp
            fastadict[ensp] = ""
        # retrieve the AA sequence belonging to individual ESPN identifiers and store as value for ESPN key
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
            if peptide in peptidedict[protID]:
                pass
            else:
                peptidedict[protID][peptide] = ""
        else:
            peptidedict[protID] = {peptide: ""}

    # now map peptide onto full length sequence, then get absolute phosphosite locations
    id_pos_res = {}
    c_sequence_not_found=0
    for protID in peptidedict.keys():
        for peptide in peptidedict[protID]:

            # first make peptide upper case to match to FASTA sequence
            UPPERpeptide = peptide.upper()
            # print UPPERpeptide

            # get parent protein sequence
            sequence = fastadict[protID]
            try:
                peptideindex = sequence.index(UPPERpeptide)
            except:
                c_sequence_not_found+=1
            # print peptideindex

            x = 0
            for letter in peptide:
                if letter.islower():
                    if letter == "s" or letter == "t" or letter == "y":

                        phoslocation = peptideindex + x + 1
                        phosresidue = sequence[phoslocation - 1]

                        if protID in id_pos_res.keys():
                            id_pos_res[protID][phoslocation] = phosresidue
                        else:
                            id_pos_res[protID] = {phoslocation: phosresidue}

                    else:
                        # non phosho modification
                        pass
                else:
                    pass

                x += 1
    sys.stdout.write('sequences not found in fasta: '+str(c_sequence_not_found))
    return id_pos_res


# return a list of dictionaries
# Usage l[index][column name]
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
            sys.stderr.write("Error in data file")
            sys.stderr.write(columns)
            sys.stderr.write(fields)
            raise

        for i in range(len(columns)):
            try:
                instance[columns[i]] = fields[i]
            except IndexError:
                instance[columns[i]] = ''

    # sys.stderr.write("No. of entries in file input: %d" % len(l))
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

            if only_leading and not Id in leading_protein_ids:
                continue

            if Id in id_pos_res:
                id_pos_res[Id][pos] = aa
            else:
                id_pos_res[Id] = {pos: aa}

    return id_pos_res

def readRunessitesfile(sitesfile):
    id_pos_res = {}
    f = open(sitesfile, 'r')
    if f:
        data = f.readlines()
        f.close()
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
    else:
        sys.stderr.write("Could not open site file: %s" % sitesfile)
        sys.exit()

    return id_pos_res

def readMCMCssitesfile(sitesfile):
    id_pos_res = {}
    f = open(sitesfile, 'r')
    if f:
        data = f.readlines()
        f.close()
        for line in data:
            tokens = line.split(' ')
            id = tokens[0]
            if (res_pos!= '') or (res_pos!= 'M_'):
                res_pos = tokens[2]
                pos = int(res_pos[2:])
                res = res_pos[0]
                if id in id_pos_res:
                    id_pos_res[id][pos] = res
                else:
                    id_pos_res[id] = {pos: res}
    else:
        sys.stderr.write("Could not open site file: %s" % sitesfile)
        sys.exit()

    return id_pos_res


'''
#Alias hashes
def readAliasFiles(organism, datadir):
	alias_hash = {}
	desc_hash = {}

	# Read alias db
	try:
		alias_db = myPopen('gzip -cd %s/%s.alias_best.tsv.gz'%(datadir, organism))
		for line in alias_db:
			(taxID, seqID, alias) = line.strip().split('\t')[:3]
			alias_hash[seqID] = alias
	except:
		sys.stderr.write("No aliases available for organism: '%s'\n"%organism)

	# Read desc db
	try:
		desc_db = myPopen('gzip -cd %s/%s.text_best.tsv.gz'%(datadir, organism))
		for line in desc_db:
			(taxID, seqID, desc) = line.split('\t')[:3]
			desc_hash[seqID] = desc
	except:
		sys.stderr.write("No descriptions available for organism: '%s'\n"%organism)

	return alias_hash, desc_hash
'''


###  my Alias hashes
def readAliasFiles(organism, datadir,map_group_to_domain):
    alias_hash = {}
    desc_hash = {}
    name_hash = {}
    names_list=[]
    for tree in map_group_to_domain:
        for pred in map_group_to_domain[tree]:
            #names_list.append(map_group_to_domain[tree][pred])
            for name in map_group_to_domain[tree][pred]:
                if name not in names_list:
                    names_list.append(name)
    # Read desc db
    try:
        desc_db = myPopen('gzip -cd %s/%s.text_best.v9.0.tsv.gz' % (datadir, organism))
        # desc_db = myPopen('gzip -cd %s/%s.protein.info.v12.0.txt.gz'%(datadir, organism))
        desc_db = desc_db.split('\n')
        for line in desc_db:
            # print(line)
            (taxID, seqID, desc) = line.split('\t')[:3]
            desc_hash[seqID] = desc
    except:
        sys.stderr.write("No descriptions available for organism: '%s'\n" % organism)

    try:
        # desc_db = myPopen('gzip -cd %s/%s.text_best.v9.0.tsv.gz'%(datadir, organism))
        name_db = myPopen('gzip -cd %s/%s.protein.aliases.v12.0.txt.gz' % (datadir, organism))
        name_db = name_db.split('\n')
        name_hash={}
        for line in name_db:
            # print(line)
            (seqID, alias, source) = line.split('\t')
            key=seqID[5:]
            if key in name_hash.keys():
                if alias in name_hash[key]:
                    continue
                else:
                    name_hash[key].append(alias)
            else:
                name_hash[key]=[alias]
    except:
        sys.stderr.write("No names available for organism: '%s'\n" % organism)
    for id in name_hash.keys():
        flag=False
        for name in name_hash[id]:
            if name in names_list:
                alias_hash[id]=name
                flag=True
        if not flag:
            alias_hash[id]=id

    return alias_hash, desc_hash, name_hash


def ReadLines(fname):
    f = open(fname)
    lines = f.readlines()
    f.close()
    return lines


def WriteString(fname, s):
    f = open(fname, 'w')
    f.write(s)
    f.close()



def mapPeptides2STRING(blastDir, organism, fastafilename, id_pos_res, id_seq, number_of_processes, datadir, fast=False,
                       leave_intermediates=False):
    sys.stderr.write("Mapping using blast\n")
    incoming2string = {}
    string2incoming = {}

    # Speedup, only blast sequences with site specified
    blast_tmpfile = tempfile.NamedTemporaryFile(mode='w',delete=False)
    if id_pos_res == {}:
        for id in id_seq:
            blast_tmpfile.write('>' + id + '\n' + id_seq[id] + '\n')
    else:
        for id in id_pos_res:
            try:
                blast_tmpfile.write('>' + id + '\n' + id_seq[id] + '\n')
            except:
                sys.stderr.write("No sequence available for '%s'\n" % id)
    blast_tmpfile.flush()

    blastDB = os.path.join(datadir, "%s.protein.sequences.v12.0.fa" % (organism)).replace(' ', '\\ ')

    # Check if blast database is actually initialized, if not: do it
    if not os.path.isfile(blastDB + '.pin'):
        command = f"{blastDir.rsplit('/', 1)[0]}/makeblastdb -in {blastDB} -out {blastDB} -parse_seqids -dbtype prot"
        sys.stderr.write("Looks like blast database is not initialized, trying to run:\n%s\n" % command)
        myPopen(command)

    #command = "%s -F F -a %s -p blastp -e 1e-10 -m 8 -d %s -i %s | sort -k12gr" % (blastDir, number_of_processes, blastDB, blast_tmpfile.name)
    command = "%sblastp -num_threads %s -evalue 1e-10 -db %s -query %s -outfmt 6 | sort -k12gr" % (blastDir, number_of_processes, blastDB, blast_tmpfile.name)

    # to save time - jhkim
    if fast and os.path.isfile(fn_blast_output):
        blast_out = ReadLines(fn_blast_output)
    else:
        blast_out = myPopen(command)
        if leave_intermediates:
            WriteString(fn_blast_output, "".join(blast_out))

    blast_out=blast_out.split('\n')
    for line in blast_out:
        #sys.stderr.write(line+'\n')
        if line == '':
            continue
        tokens = line.split('\t')
        incoming = tokens[0]

        if incoming not in incoming2string:
            # get rid of organism prefix
            string = tokens[1].replace("%s." % organism, "")
            identity = float(tokens[2])
            evalue = float(tokens[-2])

            if string in string2incoming:
                sys.stderr.write('Best hit for ' + incoming + ' is not reciprocal\n')
            if identity < 90:
                #sys.stderr.write('Best hit for ' + incoming + ' has only ' + format.fix(identity, 2) + ' %identity\n')
                sys.stderr.write('Best hit for ' + incoming + ' has only {:.2f} %identity\n'.format(identity))
            if evalue > 1e-40:
                #sys.stderr.write('Best hit for ' + incoming + ' has high E-value ' + format.sci(evalue, 2) + ' \n')
                sys.stderr.write('Best hit for ' + incoming + ' has high E-value {:.2e} \n'.format(evalue))
            if incoming in incoming2string:
                incoming2string[incoming][string] = True
            else:
                incoming2string[incoming] = {string: True}
            if string in string2incoming:
                string2incoming[string][incoming] = True
            else:
                string2incoming[string] = {incoming: True}
        else:
            pass
    return incoming2string, string2incoming

# Random mapping of incoming identifiers to string
def mapRandom(id_seq):
    sys.stderr.write("Mapping random\n")
    import random

    incoming2string = {}
    string2incoming = {}

    stringIDs = []
    file = open('%s/%s.protein.sequences.fa' % (datadir, organism), 'r')
    data = file.readlines()
    file.close()

    for line in data:
        if line[0] == '>':
            name = line[1:-1]
            stringIDs.append(name)
            string2incoming[name] = {}
    max = len(stringIDs) - 1

    for incoming in id_seq:
        if incoming not in incoming2string:
            int = random.randint(0, max)
            string = stringIDs[int]
            incoming2string[incoming] = {string: True}
            if string in string2incoming:
                string2incoming[string][incoming] = True
            else:
                string2incoming[string] = {incoming: True}
    return incoming2string, string2incoming


# In case we run NetworKIN on the same sequence set we use in STRING, we can skip the mapping by blast
def mapOne2one(id_seq):
    sys.stderr.write("Mapping one2one\n")
    incoming2string = {}
    string2incoming = {}
    for incoming in id_seq:
        incoming2string[incoming] = {incoming: True}
        string2incoming[incoming] = {incoming: True}
    return incoming2string, string2incoming


# In case we have a better mapping then we can expect from blasting, we can use an external file to do so
# file format:
# incoming ID -> STRING ID
def mapFromFile(filename):
    sys.stderr.write("Mapping using external mapping file\n")
    incoming2string = {}
    string2incoming = {}

    command = "cat %s" % (filename)
    try:
        mappingFile = myPopen(command)
    except:
        sys.stderr.write("Going to sleep, crashed with '%s'\n" % command)
        time.sleep(3600)

    for line in mappingFile:
        if (re.match('^#', line)):
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


# Load the precalculated STRING network file
def loadSTRINGdata(string2incoming, datadir,alias_hash, number_of_processes):
    # command = 'gzip -cd %s/%s.string_data.tsv.gz'%(datadir, organism)

    # fn_bestpath = "%s/%s.string_000_%04d_%04d.tsv.gz" % (os.path.join(datadir, "string_data"), organism, dPenalty[organism]["hub penalty"], dPenalty[organism]["length penalty"])
    fn_bestpath = "%s/%s.bestpath_0340_0950.v9.tsv.gz" % (os.path.join(datadir, "string_data"), organism)
    #fn_bestpath = "%s/%s.protein.links.v12.0.txt.gz" % (os.path.join(datadir, "string_data"), organism)

    if not os.path.isfile(fn_bestpath):
        sys.stderr.write("Best path file does not exist: %s" % fn_bestpath)

    '''	    
	command = "gzip -cd %s" % fn_bestpath

	try:
		data = myPopen(command)
	except:
		sys.stderr.write("Error loading STRING data using '%s', sleeping fo 1h.\n"%command)
		time.sleep(3600)
	'''

    tree_pred_string_data = {}

    # for memory efficiency
    import gzip
    f = gzip.open(fn_bestpath)
    f = f.read().decode('utf-8')
    #print(f)
    f = f.split('\n')
    #print(f)
    c=0
    for line in f:
        #print(line)
        c+=1
        if line == '':
            continue
        # line = line.strip()
        #tokens = line.split(' ')
        tokens = line.split('\t')
        #print(len(tokens))
        if len(tokens) == 8:
            (tree, group, name, string1, string2, stringscore, stringscore_indirect, path) = tokens
        if len(tokens) == 7:
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
            name = alias_hash[string1]
        '''
		if len(tokens) == 8:
			(tree, group, name, string1, string2, stringscore, stringscore_indirect, path) = tokens
		elif len(tokens) == 7:
			(tree, group, name, string1, string2, stringscore, stringscore_indirect) = tokens
			path = ""	# path to itself,  we will miss the path information
		elif len(tokens) == 6:
			(name, string1, string2, stringscore, stringscore_indirect, path) = tokens
		elif len(tokens) == 5:
			(name, string1, string2, stringscore, stringscore_indirect) = tokens
			path = ""	# path to itself,  we will miss the path information
		'''

        if string2 in string2incoming:
            if string2 in tree_pred_string_data:
                tree_pred_string_data[string2][string1] = {"_name": name}
            else:
                tree_pred_string_data[string2] = {string1: {"_name": name}}
            if options.path == "direct":
                tree_pred_string_data[string2][string1]["_score"] = float(stringscore)
            elif options.path == "indirect":
                tree_pred_string_data[string2][string1]["_score"] = float(stringscore_indirect)  # Use indirect path
            else:
                raise "Path information should be either direct or indirect."
            if len(tokens) == 9:
                tree_pred_string_data[string2][string1]["_path"] = path
            elif (len(tokens) == 6) or (len(tokens) == 3):
                tree_pred_string_data[string2][string1]["_path"] = 'notdef'
                '''
        elif string1 in string2incoming:
            if string1 in tree_pred_string_data:
                tree_pred_string_data[string1][string2] = {"_name": name}
            else:
                tree_pred_string_data[string1] = {string2: {"_name": name}}

            if options.path == "direct":
                tree_pred_string_data[string1][string2]["_score"] = float(stringscore)
            elif options.path == "indirect":
                tree_pred_string_data[string1][string2]["_score"] = float(stringscore_indirect)  # Use indirect path
            else:
                raise "Path information should be either direct or indirect."
            if len(tokens) == 9:
                tree_pred_string_data[string1][string2]["_path"] = path
            elif (len(tokens) == 6) or (len(tokens) == 3):
                tree_pred_string_data[string1][string2]["_path"] = 'notdef'
                '''
        else:
            pass

    # f.close()
    sys.stderr.write("network edges: %s \n" % c)
    return tree_pred_string_data


def InsertValueIntoMultiLevelDict(d, keys, value):
    for i in range(len(keys) - 1):
        # if not d.has_key(keys[i]):
        if not (keys[i] in d.keys()):
            d[keys[i]] = {}
        d = d[keys[i]]

    if not (keys[-1] in d.keys()):
        d[keys[-1]] = []
    d[keys[-1]].append(value)


def ReadGroup2DomainMap(path_group2domain_map):
    map_group2domain = {}  # KIN   group   name
    use_hannos=True
    if use_hannos:
        f = open("data/hanno_group_human_protein_name_map.tsv", "r")
        for line in f.readlines():
            tokens = line.split()
            name = tokens[3]
            #name = tokens[2]
            InsertValueIntoMultiLevelDict(map_group2domain, tokens[:2], name)
        f.close()
    else:
        f = open(path_group2domain_map, "r")
        for line in f.readlines():
            tokens = line.split()
            name = tokens[2]
            InsertValueIntoMultiLevelDict(map_group2domain, tokens[:2], name)
        f.close()
    return map_group2domain


def SetValueIntoMultiLevelDict(d, keys, value):
    for i in range(len(keys) - 1):
        # if not d.has_key(keys[i]):
        if keys[i] not in d:
            d[keys[i]] = {}
        d = d[keys[i]]
    if (keys[-1] in d) and type(d[keys[-1]]) != type(value):
        sys.stderr.write("Causion: multi-dict already has value and try to assign a value of different type")
        pass
    if (keys[-1] in d):
        if d[keys[-1]] != value:
            # sys.stderr.write("This operation replaces a value (%s)" % " ".join(map(lambda(x):str(x), keys)))
            sys.stderr.write("This operation replaces a value ({})".format(' '.join(map(str, keys))))
    d[keys[-1]] = value



def printResult(id_pos_tree_pred, tree_pred_string_data, incoming2string, string_alias,name_hash, string_desc, organism, mode,
                dir_likelihood_conversion_tbl, map_group2domain,res_dir, fn_fasta):
    ALPHA = ALPHAS[organism]
    species = dSpeciesName[organism]

    fas = open(fn_fasta, 'r')
    csv_filename = os.path.join(res_dir, f'{fas.name}.result.tsv')
    predictions = []

    dLRConvTbl = {}
    for fname in glob.glob(os.path.join(dir_likelihood_conversion_tbl, "conversion_tbl_*_smooth*")):
        score_src, species_of_conversion_table, tree, player_name = \
        re.findall("conversion_tbl_([a-z]+)_smooth_([a-z]+)_([A-Z0-9]+)_([a-zA-Z0-9_/-]+)",
                   os.path.basename(os.path.splitext(fname)[0]))[0]
        # Normalise the legacy filename tag to a uniform in-memory key.
        # Conversion table files are named either "…_string_…" (STRING co-expression
        # likelihood) or "…_<motif-scorer>_…" (motif probability likelihood).
        # Map all non-"string" tags to "motif" so downstream lookups are scorer-agnostic.
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
            sys.stderr.write("Conversion table %s %s %s %s\n" % (
            species_of_conversion_table, tree, player_name, score_src))

    c_np = 0
    c_map = 0
    c_notmap = 0
    c_string = 0
    c_stringscore=0
    c_if=0
    c_else=0
    c_res = 0
    pred_not_mapped=[]
    name_not_mapped=[]
    # For each scored protein ID
    for id in id_pos_tree_pred:
        c_np += 1
        # We have a mapping to STRING
        if id in incoming2string:
            c_map+=1
            # For each predicted position
            for pos in id_pos_tree_pred[id]:
                # For each of the trees (KIN, SH@ etc.)
                for tree in id_pos_tree_pred[id][pos]:
                    score_results = {}
                    # For each single classifier KIN eg.
                    for pred in id_pos_tree_pred[id][pos][tree]:
                        #print(pred)
                        # For each mapped sequence
                        for string1 in incoming2string[id]:
                            if string1 in string_alias:
                                bestName1 = string_alias[string1]
                            else:
                                bestName1 = 'notdef'
                            if string1 in string_desc:
                                desc1 = string_desc[string1]
                            else:
                                desc1 = 'notdef'
                            # if string in string network
                            if string1 in tree_pred_string_data:
                                c_string+=1
                                (res, peptide, motifScore) = id_pos_tree_pred[id][pos][tree][pred]
                                for string2 in tree_pred_string_data[string1]:

                                    if string2 in string_alias:
                                        bestName2 = string_alias[string2]
                                    else:
                                        bestName2 = 'notdef'
                                    if string2 in string_desc:
                                        desc2 = string_desc[string2]
                                    else:
                                        desc2 = 'notdef'
                                    stringScore = tree_pred_string_data[string1][string2]["_score"]
                                    c_stringscore+=1
                                    #path = tree_pred_string_data[string1][string2]["_path"]
                                    path='notdef'
                                    name = tree_pred_string_data[string1][string2]["_name"]  # string2 = kinase
                                    '''
                                    names=name_hash[string2]
                                    #names=[name]
                                    '''
                                    try:
                                        map_names=map_group2domain[tree][pred]
                                    except:
                                        pred_not_mapped.append(pred)



                                    '''
                                    flag=False
                                    for mn in map_names:
                                        pattern = rf"(?i)\b{re.escape(mn)}\w*\b"
                                        l = [n for n in names if re.search(pattern,n)]
                                    for n in names:
                                        pattern = rf"(?i)\b{re.escape(n)}\w*\b"
                                        l1 = [mn for mn in map_names if re.search(pattern,mn)]
                                    if len(l)>0 or len(l1)>0 or any(n in names for n in map_names):
                                        flag=True
                                    '''
                                    #flag=any(n in names for n in map_names)
                                    #if not flag:
                                    #    name_not_mapped.append(name)
                                    # sys.stderr.write("%s %s %s\n" % (tree, pred, name))

                                    if (tree not in map_group2domain.keys()) or (pred not in map_group2domain[tree].keys()) or (name not in map_group2domain[tree][pred]):
                                        #if (tree not in map_group2domain.keys()) or (pred not in map_group2domain[tree].keys()) or not flag:
                                        if name not in map_group2domain[tree][pred]:
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
                                                raise "This species is not supported"

                                            likelihood_motif = 1
                                            likelihood_string = ConvertScore2L(stringScore, conversion_tbl_string)
                                            unified_likelihood = likelihood_motif * likelihood_string
                                            networkinScore = unified_likelihood

                                            # NetworKIN result
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
                                        c_else+=1
                                        # sys.stderr.write("%s %s\n" % (string1, string2))

                                        if species == "human":
                                            if tree in ["1433", "BRCT", "WW", "PTB", "WD40", "FHA"]:
                                                conversion_tbl_motif = dLRConvTbl[species]["SH2"]["general"][
                                                    "motif"]
                                                conversion_tbl_string = dLRConvTbl[species]["SH2"]["general"]["string"]
                                            else:
                                                if name in dLRConvTbl[species][tree]:
                                                    conversion_tbl_motif = dLRConvTbl[species][tree][name][
                                                        "motif"]
                                                    conversion_tbl_string = dLRConvTbl[species][tree][name]["string"]
                                                else:
                                                    conversion_tbl_motif = dLRConvTbl[species][tree]["general"][
                                                        "motif"]
                                                    conversion_tbl_string = dLRConvTbl[species][tree]["general"][
                                                        "string"]
                                        elif species == "yeast":
                                            if name in dLRConvTbl[species][tree]:
                                                conversion_tbl_motif = dLRConvTbl[species][tree][name][
                                                    "motif"]
                                                conversion_tbl_string = dLRConvTbl[species][tree][name]["string"]
                                            else:
                                                conversion_tbl_motif = dLRConvTbl[species][tree]["general"][
                                                    "motif"]
                                                conversion_tbl_string = dLRConvTbl[species][tree]["general"]["string"]
                                        else:
                                            raise "This species is not supported"

                                        likelihood_motif = ConvertScore2L(motifScore,
                                                                               conversion_tbl_motif)
                                        likelihood_string = ConvertScore2L(stringScore, conversion_tbl_string)
                                        unified_likelihood = likelihood_motif * likelihood_string
                                        networkinScore = unified_likelihood

                                        # NetworKIN result
                                        c_res+=1
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
                    '''               
                    if mode == 'network':
                        highestScore = sorted(score_results.keys(), reverse=True)[0]

                        if len(score_results[highestScore]) > 1:
                            index = random.randint(0, len(score_results[highestScore]) - 1)
                            sys.stdout.write((score_results[highestScore][index]))
                        else:
                            sys.stdout.write((score_results[highestScore][0]))
                        pass
                    else:
                        for score in sorted(score_results.keys(), reverse=True):
                            sys.stdout.write("".join(score_results[score]))
                    '''
        else:
            c_notmap+=1

        #print('preds not mapped:')
        #print(pred_not_mapped)

        print('names not mapped:')
        print(len(np.unique(name_not_mapped)))

    # --- False-negative recovery ---
    # Build motif_score_dict: {(kinase_string_id, substrate_string_id): motif_score}
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

    # Build node_index and distance matrix from STRING network data.
    # tree_pred_string_data already holds best-path STRING scores, so we use
    # d = 1/score - 1 to convert each score to a distance without re-running
    # Floyd–Warshall on the raw edge graph.
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

    # ORIGINAL: csv.writer rows written inline
    # NEW: write all predictions at once via output module
    write_output(predictions, csv_filename)

    return c_np,c_map,c_notmap,c_string,c_stringscore,c_if,c_else,c_res


# MAIN
def Main():
    res_dir = 'results'
    if not os.path.exists(res_dir):
        os.makedirs(res_dir)

    tmpdir = 'tmp'
    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)

    sys.stderr.write("Reading fasta input file\n")
    id_seq = readFasta(fastafile)
    if options.verbose:
        sys.stderr.write("%s sequences loaded\n" % len(id_seq.keys()))

    if sitesfile:
        sys.stderr.write("Reading phosphosite file\n")
        input_type = CheckInputType(sitesfile)
        if input_type == NETWORKIN_SITE_FILE:
            id_pos_res = readPhosphoSites(sitesfile)
        elif input_type == PROTEOME_DISCOVERER_SITE_FILE:
            id_pos_res = readPhosphoSitesProteomeDiscoverer(fn_fasta, sitesfile)
        elif input_type == MAX_QUANT_DIRECT_OUTPUT_FILE:
            id_pos_res = readPhosphoSitesMaxQuant(sitesfile)
        elif input_type == RUNES_SITE_FILE:
            id_pos_res = readRunessitesfile(sitesfile)
        elif input_type == MS_MCMC_FILE:
            id_pos_res = readMCMCssitesfile(sitesfile)
    else:
        id_pos_res = {}

    # ORIGINAL: sites loaded exclusively from user-provided flat file (phospho.tsv / *.dat dumps).
    # NEW: fetch live reference phosphorylation-site data from PhosphoSitePlus (7-day cache).
    # phosphosite_df is available here for downstream enrichment or validation against id_pos_res.
    sys.stderr.write("Fetching live phosphorylation site reference data\n")
    phosphosite_df = fetch_phosphosite(refresh=options.refresh)  # noqa: F841


    if organism == "9606":
        path_group2domain_map = os.path.join(options.datadir, "group_human_protein_name_map.tsv")
    elif organism == "4932":
        path_group2domain_map = os.path.join(options.datadir, "group_yeast_KIN.tsv")

    sys.stderr.write("Loading aliases and descriptions\n")
    map_group2domain = ReadGroup2DomainMap(path_group2domain_map)
    (string_alias, string_desc, name_hash) = readAliasFiles(args[0], options.datadir, map_group2domain);

    # Default way of mapping using BLAST
    # incoming2string, string2incoming = mapPeptides2STRING(blastDir, organism, fastafile.name, id_pos_res, id_seq, options.threads, options.datadir, options.fast, options.leave)
    incoming2string, string2incoming = mapPeptides2STRING(blastDir, organism, fastafile.name, id_pos_res, id_seq,
                                                          options.threads, options.datadir)


    # Hack for random mapping to proteins
    # incoming2string, string2incoming = mapRandom(id_seq)

    # Use if a mapping file for the input can be provided
    # incoming2string, string2incoming = mapFromFile("/home/red1/hhorn/projects2/2012_03_22_Jesper/ensmusp_ensp.tsv")

    # Used if only ensembl of the same version used
    # incoming2string, string2incoming = mapOne2one(id_seq)

    # ORIGINAL: tree_pred_string_data loaded exclusively from pre-processed local gzip file:
    #   fn_bestpath = "%s/%s.bestpath_0340_0950.v9.tsv.gz" % (os.path.join(datadir, "string_data"), organism)
    # NEW: fetch live STRING functional association network via REST API (7-day cache).
    # string_df (protein_a, protein_b, combined_score) is available here for downstream use.
    # The pre-processed local file is still loaded below for the pipeline's internal scoring format.
    sys.stderr.write("Fetching live STRING network data\n")
    all_uniprot_ids = list(id_seq.keys())
    string_df = fetch_string_network(proteins=all_uniprot_ids, refresh=options.refresh)  # noqa: F841

    # Load the STRING network data
    sys.stderr.write("Loading STRING network\n")
    tree_pred_string_data = loadSTRINGdata(string2incoming, options.datadir,string_alias, options.threads)
    # Run motif scoring
    sys.stderr.write("Running motif scorer\n")
    id_pos_tree_pred = score_sequences(id_seq, id_pos_res)

    # Writing result to STDOUT
    sys.stderr.write("Writing results\n")
    sys.stdout.write(
        "#Name\tPosition\tTree\tMotif Group\tKinase/Phosphatase/Phospho-binding domain\tNetworKIN score\tMotif probability\tSTRING score\tTarget STRING ID\tKinase/Phosphatase/Phospho-binding domain STRING ID\tTarget description\tKinase/Phosphatase/Phospho-binding domain description\tTarget Name\tKinase/Phosphatase/Phospho-binding domain Name\tPeptide sequence window\tIntermediate nodes\n")
    if options.path == "direct":
        dir_likelihood_conversion_tbl = os.path.join(options.datadir, "likelihood_conversion_table_direct")
    elif options.path == "indirect":
        dir_likelihood_conversion_tbl = os.path.join(options.datadir, "likelihood_conversion_table_indirect")
    else:
        raise "Path information should be either direct or indirect."
    #print(id_pos_tree_pred)
    #c_np,c_map,c_notmap,c_string,c_res = printResult(id_pos_tree_pred, tree_pred_string_data, incoming2string, string_alias, string_desc, string_names,
    #            args[0], options.mode, dir_likelihood_conversion_tbl, map_group2domain, res_dir, fn_fasta)
    #print(id_pos_tree_pred['SRC'][128]['KIN']['Abl_group']
    c_np,c_map,c_notmap,c_string,c_stringscore,c_if,c_else,c_res =printResult(id_pos_tree_pred, tree_pred_string_data, incoming2string, string_alias,name_hash, string_desc, args[0],
                options.mode, dir_likelihood_conversion_tbl, map_group2domain, res_dir, fn_fasta)

    sys.stdout.write('c_np = '+str(c_np)+'\n'+'c_map = '+str(c_map)+'\n'+'c_notmap = '+str(c_notmap)+'\n'+'c_string = '+str(c_string)+'\n'+'c_stringscore = '+str(c_stringscore)+'\n'+'c_if = '+str(c_if)+'\n'+'c_else = '+str(c_else)+'\n'+'c_res = '+str(c_res)+'\n')

    return


if __name__ == '__main__':
    # BLAST
    try:
        blastDir = os.environ['BLAST_PATH']
    except:
        blastDir = ""
    # Binary option removed; motif scoring now uses the Python atlas scorer

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
    # if os.environ.has_key("TMPDIR"):
    #	parser.add_option("--tmp", dest="tmpdir", default=os.environ["TMPDIR"],
    #									help="location for the temporary files [default: %default]")
    # else:
    #	print >> sys.stderr, "TMPDIR environmental variable is not defined. Please define this variable or specify tmp directory by using --tmp command line option"
    #	sys.exit()

    global options
    (options, args) = parser.parse_args()

    if options.active_threads < 2:
        parser.error("Number of active thread (--nt) is less than 2")
    # tempfile.tempdir= options.tmpdir

    # ORGANISM
    try:
        organism = args[0]
    except:
        parser.error("Organism not defined!")

    # SEQUENCE FILE
    try:
        fn_fasta = args[1]
        fn_blast_output = "%s.%s.blast.out" % (fn_fasta, organism)
        fastafile = open(fn_fasta, 'r')
    except:
        sys.stderr.write("%s" % args)
        parser.error("FASTA-file not defined!")

    # SITES FILE
    try:
        sitesfile = args[2]
    except:
        sitesfile = False

    # BLAST
    if options.blast:
        blastDir = options.blast

    # Show runtime parameters
    if (options.verbose):
        sys.stderr.write(
            '\nPredicting using parameters as follows:\nOrganism:\t%s\nFastaFile:\t%s\n' % (organism, fn_fasta))
        sys.stderr.write('Threads:\t%s\nCompress:\t%s\n' % (options.threads, options.compress))
        if options.string_for_uncovered:
            sys.stderr.write("Use STRING likelihood when kinases are not covered by the motif atlas.\n")
        if sitesfile:
            sys.stderr.write('Sitesfile:\t%s\n' % sitesfile)
        else:
            sys.stderr.write('No sites-file given, predicting on all S,T,Y residues.\n')
        if options.mode:
            sys.stderr.write('Mode:\t\t%s\n' % options.mode)
        sys.stderr.write("Blast dir: %s" % blastDir)
        sys.stderr.write('\n')
    Main()
