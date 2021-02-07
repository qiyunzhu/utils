#!/usr/bin/env python3
"""Extract hierarchical mappings and descriptions from the MetaCyc database.

Usage:
    python me.py metacyc_data_dir <protein.list>

Notes:
    A local copy of the MetaCyc database is needed for running this script.
    - metacyc_data_dir: The "data" directory in a MetaCyc release.
    - protein.list (optional): A list of protein IDs. Only those proteins and
      their higher-level relationships will be included.
    Tested and working with MetaCyc release 23.0.

References:
    The latest MetaCyc paper (Caspi et al., 2020)
    https://academic.oup.com/nar/article/48/D1/D445/5581728
"""

import sys
import re


__author__ = 'Qiyun Zhu'
__license__ = 'BSD-3-Clause'
__version__ = '0.0.1-dev'
__email__ = 'qiyunzhu@gmail.com'


# metacyc directory
mcdir = None

# encoding of MetaCyc database
mckw = {'encoding': 'Windows-1252'}


def parse_dat(name, ids=None, keys=[]):
    """Parse MetaCyc .dat file.

    Parameters
    ----------
    name : str
        Stem filename of target .dat file.
    ids : set of str, optional
        Entries to include in extraction.
    keys : list of str, optional
        Metadata keys to extract from each entry.

    Returns
    -------
    dict of dict
        Extracted entry to metadata dictionary.
    """
    data = {}
    with open(f'{mcdir}/{name}.dat', 'r', **mckw) as f:
        cid = None
        for line in f:
            line = line.rstrip()
            if line.startswith('UNIQUE-ID - '):
                uid = line[12:]
                if not ids or uid in ids:
                    cid = uid
                    data[cid] = {}
                else:
                    cid = None
                    continue
            if line.startswith('COMMON-NAME - '):
                if 'name' in data[cid]:
                    raise ValueError(f'Duplicate name for {cid}.')
                data[cid]['name'] = line[14:]
            for key in keys:
                if line.startswith(key.upper() + ' - '):
                    data[cid].setdefault(key, []).append(line[len(key) + 3:])
    return data


def write_names(data, fname):
    """Write entry-to-name mapping.

    Parameters
    ----------
    data : dict
        Main data structure.
    fname : str
        Output file name.
    """
    with open(fname, 'w') as f:
        for uid, datum in sorted(data.items()):
            if 'name' in datum:
                print(uid, datum['name'], sep='\t', file=f)


def write_map(data, key, fname):
    """Write entry-to-parent(s) mapping.

    Parameters
    ----------
    data : dict
        Main data structure.
    key : str
        Metadata key.
    fname : str
        Output file name.

    Returns
    -------
    set of str
        Parent entries.
    """
    pids = set()
    with open(fname, 'w') as f:
        for uid, datum in sorted(data.items()):
            if key in datum:
                print(uid, '\t'.join(datum[key]), sep='\t', file=f)
                pids.update(datum[key])
    return pids


def parse_pwy_genes(pwys=None, genes=None):
    """Parse MetaCyc pwy-genes.dat file.

    Parameters
    ----------
    pwys : set of str, optional
        Pathways to include.
    genes : set of str, optional
        Genes to include.

    Returns
    -------
    dict of list
        Pathway-to-genes dictionary.
    dict of list
        Gene-to-pathway(s) dictionary.

    Notes
    -----
    Pathway-to-genes dictionary may contain genes that are not in the provided
    gene list (because the dictionary reflects compositions of whole pathways).
    """
    pwy2genes = {}
    with open(f'{mcdir}/pwy-genes.dat', 'r', **mckw) as f:
        # merge lines and trim beginning and end
        dat = f.read().replace('\n', '')[20:-4]
    for rec in dat.split(') ('):
        # occassionall there are '|' pairs surrounding certain entries, and
        # with these entries one needs to replace spaces with dashes
        for x in re.findall(r'\|[^\|]+\|', rec):
            rec = rec.replace(x, x.strip('|').replace(' ', '-'))
        x = rec.split()
        if pwys and x[0] not in pwys:
            continue
        pwy2genes[x[0]] = x[1:]
    gene2pwy = {}
    for pwy, gs in pwy2genes.items():
        for g in gs:
            if genes and g not in genes:
                continue
            gene2pwy.setdefault(g, []).append(pwy)
    return pwy2genes, gene2pwy


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)

    # metacyc directory
    global mcdir
    mcdir = sys.argv[1]

    # proteins to include (optional)
    pts = None
    if len(sys.argv) > 2:
        with open(sys.argv[2], 'r') as f:
            pts = set(f.read().splitlines())

    # parse protein data
    data = parse_dat('proteins', pts, ['catalyzes', 'gene', 'go-terms'])
    print(f'Proteins: {len(data)}.')
    write_names(data, 'protein_name.txt')
    genes = write_map(data, 'gene', 'protein-to-gene.txt')
    enzrxns = write_map(data, 'catalyzes', 'protein-to-enzrxn.txt')
    goes = write_map(data, 'go-terms', 'protein-to-go.txt')
    print(f'GO terms: {len(goes)}.')

    # parse gene data
    data = parse_dat('genes', pts and genes, [])
    print(f'Genes: {len(data)}.')
    write_names(data, 'gene_name.txt')

    # parse enzymatic reaction data
    data = parse_dat('enzrxns', pts and enzrxns, ['reaction', 'regulated-by'])
    print(f'Enzymatic reactions: {len(data)}.')
    write_names(data, 'enzrxn_name.txt')
    reactions = write_map(data, 'reaction', 'enzrxn-to-reaction.txt')
    regulations = write_map(data, 'regulated-by', 'enzrxn-to-regulation.txt')

    # parse regulation data
    data = parse_dat('regulation', pts and regulations, ['regulator'])
    print(f'Regulations: {len(data)}.')
    regulators = write_map(data, 'regulator', 'regulation-to-regulator.txt')
    print(f'Regulators: {len(regulators)}.')

    # parse reaction data
    data = parse_dat('reactions', pts and reactions, [
        'ec-number', 'in-pathway', 'left', 'right'])
    print(f'Reactions: {len(data)}.')
    write_names(data, 'reaction_name.txt')
    ecs = write_map(data, 'ec-number', 'reaction-to-ec.txt')
    print(f'EC numbers: {len(ecs)}.')
    pathways = write_map(data, 'in-pathway', 'reaction-to-pathway.txt')
    lefts = write_map(data, 'left', 'reaction-to-left_compound.txt')
    print(f'Compounds consumed: {len(lefts)}.')
    rights = write_map(data, 'right', 'reaction-to-right_compound.txt')
    print(f'Compounds produced: {len(rights)}.')

    # parse compound data
    compounds = regulators.union(lefts).union(rights)
    data = parse_dat('compounds', pts and compounds, [
        'types', 'inchi', 'smiles'])
    print(f'Compounds: {len(data)}.')
    write_names(data, 'compound_name.txt')
    ctypes = write_map(data, 'types', 'compound_type.txt')
    print(f'Compound types: {len(ctypes)}.')
    write_map(data, 'inchi', 'compound_inchi.txt')
    write_map(data, 'smiles', 'compound_smiles.txt')

    # parse pathway data
    data = parse_dat('pathways', pts and pathways, [
        'types', 'reaction-list', 'reaction-layout', 'super-pathways',
        'taxonomic-range'])
    print(f'Pathways: {len(data)}.')
    write_names(data, 'pathway_name.txt')
    ptypes = write_map(data, 'types', 'pathway_type.txt')
    print(f'Pathway types: {len(ptypes)}.')
    write_map(data, 'super-pathways', 'pathway-to-super_pathway.txt')
    write_map(data, 'taxonomic-range', 'pathway-to-taxonomic_range.txt')
    write_map(data, 'reaction-list', 'pathway-to-reaction_list.txt')
    write_map(data, 'reaction-layout', 'pathway-to-reaction_layout.txt')

    # parse class data
    classes = set().union(*[
        pts, genes, enzrxns, goes, reactions, regulations, ecs, pathways,
        compounds, ctypes, ptypes]) if pts else None
    data = parse_dat('classes', classes, [])
    print(f'All classes: {len(data)}.')
    write_names(data, 'all_class_name.txt')

    # parse pathway-to-genes data
    pwy2genes, gene2pwy = parse_pwy_genes(pts and pathways, pts and genes)
    print(f'Pathway-to-gene_list mappings: {len(pwy2genes)}.')
    with open(f'pathway-to-gene_list.txt', 'w') as f:
        for pwy, genes in sorted(pwy2genes.items()):
            print(pwy, '\t'.join(genes), sep='\t', file=f)
    with open(f'gene-to-pathway.txt', 'w') as f:
        for gene, pwy in sorted(gene2pwy.items()):
            print(gene, '\t'.join(pwy), sep='\t', file=f)


if __name__ == '__main__':
    main()
