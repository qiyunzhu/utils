#!/usr/bin/env python3
"""Sample, download and organize a reference genome database from NCBI RefSeq.

Usage:
    python me.py [options] -o output_dir

Notes:
    This script is ported from the following code, also written by me, with
    modifications.

    - https://github.com/qiyunlab/HGTector/blob/master/hgtector/database.py

    License: BSD 3-Clause

    - https://github.com/qiyunlab/HGTector/blob/master/LICENSE
"""

import sys
import re
import argparse
import subprocess as sp

from os import remove, makedirs
from os.path import join, isfile, isdir, getsize, basename
from shutil import rmtree
from time import sleep
import gzip
import tarfile
from tempfile import mkdtemp
import datetime

import pandas as pd


__author__ = 'Qiyun Zhu'
__license__ = 'BSD-3-Clause'
__version__ = '0.0.1-dev'
__email__ = 'qiyunzhu@gmail.com'

usage = """
  %(prog)s [options] -o OUTPUT_DIR
"""

description = """Examples:
  %(prog)s --default -o .
  %(prog)s --output <output_dir> \\
    --cats microbe --sample 1 --rank species_latin --above \\
    --reference --represent --typemater
"""

epilog = """
Sample, download and organize a reference genome database from NCBI RefSeq.
 """

server = 'rsync://ftp.ncbi.nlm.nih.gov'

arguments = [
    'Protocols',
    ['-p|--protocol', 'Apply a pre-defined protocol.', {'choices': [
        'basic', 'standard', 'extended', 'extreme']}],

    'Basic',
    ['-o|--output',  'Output directory.', {'required': True}],
    ['-c|--cats',    'Enter one or more of the following RefSeq genome '
                     'categories, separated by comma: archaea, bacteria, '
                     'fungi, invertebrate, plant, protozoa, vertebrate_'
                     'mammalian, vertebrate_other, viral. Enter "all" for '
                     'all categories, or "microbe" (default) for archaea, '
                     'bacteria, fungi, protozoa and viral.',
                     {'default': 'microbe'}],

    'Custom list',
    ['-t|--taxids',  'Provide a custom list of TaxIDs. Only genomes under '
                     'these taxonomic groups will be included.'],
    ['-g|--genoids', 'Provide a custom list of genome IDs, which are NCBI '
                     'assembly accessions (e.g., "GCF_000123456.1"), with or '
                     'without version number, or simply "G000123456".'],
    ['-e|--exclude', 'Exclude instead of include TaxIDs and genomes defined '
                     'by custom lists.', {'action': 'store_true'}],

    'Taxon sampling',
    ['-s|--sample',  'Sample up to this number of genomes per taxonomic '
                     'group at given rank (0 for none, omit for all).',
                     {'type': int}],
    ['-r|--rank',    'Taxonomic rank at which sampling will be performed '
                     '(default: species).', {'default': 'species'}],
    ['--above',      'Sampling will also be performed on ranks from the '
                     'designated one to phylum.',
                     {'action': 'store_true'}],

    'Genome sampling',
    ['--genbank',    'Also search GenBank if a genome is not found in RefSeq. '
                     'Otherwise only include RefSeq genomes',
                     {'action': 'store_true'}],
    ['--complete',   'Only include complete genomes or chromosomes.',
                     {'action': 'store_true'}],
    ['--reference',  'Include NCBI-defined reference genomes.',
                     {'action': 'store_true'}],
    ['--represent',  'Include NCBI-defined representative genomes.',
                     {'action': 'store_true'}],
    ['--typemater',  'Include NCBI-defined type material genomes.',
                     {'action': 'store_true'}],

    'Name filter',
    ['--capital',    'Organism name must be capitalized.',
                     {'action': 'store_true'}],
    ['--block',      'Ignore organism names containing any of these words. '
                     'Default: "unknown,uncultured,unidentified,unclassified,'
                     'unresolved,environmental,synthetic".',
                     {'default': 'unknown,uncultured,unidentified,'
                      'unclassified,unresolved,environmental,synthetic'}],
    ['--latin',      'Genomes must have Latinate species names, e.g. "Vibrio '
                     'cholerae", but not "Vibrio sp. 123".',
                     {'action': 'store_true'}],

    'Download',
    ['--filext',     'Genomic data file types to download (comma-separated). '
                     'Default: "genomic.fna,genomic.gff".',
                     {'default': 'genomic.fna,genomic.gff'}],
    ['--manual',     'Export URLs of sampled genomes without downloading them'
                     ', so the user may download them manually, then resume '
                     'the database pipeline.', {'action': 'store_true'}],
    ['--overwrite',  'Overwrite existing files with newly downloaded ones. '
                     'Othewise, use existing files whenever available, i.e., '
                     'resume an interrupted run.', {'action': 'store_true'}],
    ['--retries',    'Number of trials for downloading each data file '
                     '(default: 3).', {'type': int, 'default': 3}],
    ['--delay',      'Seconds between two download trials (default: 5).',
                     {'type': int, 'default': 5}],
    ['--timeout',    'Seconds before program gives up waiting (default: 60).',
                     {'type': int, 'default': 60}],
]


def main():
    """Main workflow.
    """
    Database()(parse_args(arguments))


def parse_args(args):
    """Parse command-line arguments.

    Parameters
    ----------
    args : list
        Command-line arguments.

    Returns
    -------
    dict
        arguments
    """
    cli = argparse.ArgumentParser(
        usage=usage, description=description, epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    for arg in args:
        if isinstance(arg, str):
            par = cli.add_argument_group(arg)
        else:
            flags = arg[0].split('|')
            kwargs = arg[2] if len(arg) > 2 else {}
            par.add_argument(*flags, **kwargs, help=arg[1])
    for arg in cli._actions:
        arg.metavar = '\b' if arg.required else ''
    if len(sys.argv) == 1:
        cli.print_help(sys.stderr)
        sys.exit(1)
    return cli.parse_args()


class Database(object):

    def __init__(self):
        self.capital = arguments
        self.description = description

    def __call__(self, args):
        print('Database building started at {}.'.format(
            datetime.datetime.now()))

        # create temporary directory
        self.tmpdir = mkdtemp()

        # read and validate arguments
        self.set_parameters(args)

        # retrieve taxonomy database
        self.retrieve_taxdump()

        # retrieve genome assembly summary
        self.retrieve_summary()

        # retrieve genome categories
        self.retrieve_categories()

        # filter genomes
        self.filter_genomes()

        # sort genomes by quality
        self.sort_genomes()

        # identify taxonomy of genomes
        self.filter_by_taxonomy()

        # sample genomes by taxonomy
        self.sample_by_taxonomy()

        # sample genomes by quality
        self.sample_by_quality()

        # filter genomes to sampled one
        self.filter_to_sampled()

        # download genomes
        self.download_genomes()

        # extract genomes
        self.extract_genomes()

        # identify genome lineages
        self.genome_lineages()

        # build taxonomy database
        self.build_taxdump()

        # extract genome annotations
        self.extract_genes()

        # write genome metadata
        self.genome_metadata()

        # clean up temporary directory
        rmtree(self.tmpdir)

        print('Database building finished at {}.'.format(
            datetime.datetime.now()))

    def set_parameters(self, args):
        """Workflow for validating and setting arguments.

        Parameters
        ----------
        args : dict
            command-line arguments
        """
        # load arguments
        for key, val in vars(args).items():
            setattr(self, key, val)

        # default protocols
        protocol = self.protocol
        if protocol == 'basic':
            print('Basic protocol: Sample one complete genome per genus.')
            self.sample = 1
            self.rank = 'genus'
            self.complete = True
        elif protocol == 'standard':
            print('Standard protocol: Collect all NCBI-defined reference and '
                  'representative genomes.')
            self.sample = 0
            self.reference = True
            self.represent = True
        elif protocol == 'extended':
            print('Extended protocol: Sample one genome per species that has '
                  'a Latinate name, as well as higher ranks, also include all'
                  ' reference, representative, and type material genomes.')
            self.sample = 1
            self.rank = 'species_latin'
            self.above = True
            self.reference = True
            self.represent = True
            self.typemater = True
        elif protocol == 'extreme':
            print('Extreme protocol: Sample one genome per species and '
                  'higher ranks from both RefSeq and GenBank, also include '
                  'all reference, representative, and type material genomes.')
            self.sample = 1
            self.rank = 'species'
            self.above = True
            self.genbank = True
            self.reference = True
            self.represent = True
            self.typemater = True

        # create target directories
        makedirs(self.output, exist_ok=True)
        self.down = join(self.output, 'download')
        makedirs(self.down, exist_ok=True)

    def retrieve_taxdump(self):
        """Retrieve NCBI taxdump.
        """
        fname = 'taxdump.tar.gz'
        rfile = f'pub/taxonomy/{fname}'
        lfile = join(self.down, fname)

        # download taxdump
        if not check_local_file(lfile, self.overwrite):
            print('Downloading NCBI taxonomy database...', end='', flush=True)
            run_command(f'rsync -Ltq {server}/{rfile} {self.down}')
            print(' done.')

        # read taxdump
        print('Reading NCBI taxonomy database...', end='', flush=True)
        with tarfile.open(lfile, 'r:gz') as f:
            f.extract('names.dmp', self.tmpdir)
            f.extract('nodes.dmp', self.tmpdir)
        self.taxdump = read_taxdump(self.tmpdir)
        print(' done.')
        print(f'  Total number of TaxIDs: {len(self.taxdump)}.')

    def retrieve_summary(self, genbank=False):
        """Retrieve genome assembly summary.
        """

        def get_summary(target):
            key = target.lower()
            fname = f'assembly_summary_{key}.txt'
            rfile = f'genomes/{key}/{fname}'
            lfile = join(self.down, fname)

            # download summary
            if not check_local_file(lfile, self.overwrite):
                print(f'Downloading {target} assembly summary...', end='',
                      flush=True)
                run_command(f'rsync -Ltq {server}/{rfile} {self.down}')
                print(' done.')

            # read summary
            print(f'Reading {target} assembly summary...', end='', flush=True)
            df = pd.read_table(lfile, skiprows=1, low_memory=False)
            print(' done.')
            return df

        self.df = get_summary('RefSeq')
        if self.genbank:
            self.df = pd.concat([self.df, get_summary('GenBank')])
        print(f'  Total number of genomes: {self.df.shape[0]}.')

    def retrieve_categories(self):
        """Retrieve genome categories.
        """
        # parse categories
        cats = self.cats
        if cats == 'all':
            return
        elif cats == 'microbe':
            self.cats = ['archaea', 'bacteria', 'fungi', 'protozoa', 'viral']
        else:
            self.cats = cats.split(',')
        print(f'Genome categories: {", ".join(self.cats)}')

        def get_categories(target):
            key = target.lower()

            # validate categories
            ec, out = run_command(
                f'rsync --list-only --no-motd {server}/genomes/{key}/')
            cats = [line.rsplit(None, 1)[-1] for line in out if
                    line.startswith('d') and not line.endswith('.')]
            for cat in self.cats:
                if cat not in cats:
                    raise ValueError(
                        f'"{cat}" is not a valid {target} genome category.')

            # get genome list per category
            print(f'Downloading genome list per {target} category...')
            makedirs(join(self.down, 'cats'), exist_ok=True)
            ldir = join(self.down, 'cats')
            fname = 'assembly_summary.txt'
            file_ = join(self.tmpdir, fname)

            asms = []
            for cat in self.cats:
                lfile = join(ldir, f'{key}_{cat}.txt')

                # use local file
                if check_local_file(lfile):
                    with open(lfile, 'r') as f:
                        asms_ = f.read().splitlines()

                # download remote file
                else:
                    rfile = f'genomes/{key}/{cat}/{fname}'
                    run_command(f'rsync -Ltq {server}/{rfile} {self.tmpdir}')
                    with open(file_, 'r') as f:
                        asms_ = [x.split('\t', 1)[0] for x in f.read().
                                 splitlines() if not x.startswith('#')]
                    with open(lfile, 'w') as f:
                        f.write(''.join([x + '\n' for x in asms_]))

                print(f'  {cat}: {len(asms_)}')
                asms += asms_

            if isfile(file_):
                remove(file_)
            print('Done.')
            return asms

        # filter genomes by category
        asmset = set(get_categories('RefSeq'))
        if self.genbank:
            asmset.update(get_categories('GenBank'))
        self.df = self.df[self.df['# assembly_accession'].isin(asmset)]
        print(f'  Total number of genomes in categories: {self.df.shape[0]}.')

    def filter_genomes(self):
        """Filter genomes based on genome metadata.
        """
        print('Filtering genomes...')
        n = self.df.shape[0]

        def report_diff(msg):
            nonlocal n
            n_ = self.df.shape[0]
            if n_ < n:
                print('  ' + msg.format(n - n_))
            n = n_

        # complete genomes only
        if self.complete:
            self.df = self.df[self.df['assembly_level'].isin(
                {'Complete Genome', 'Chromosome'})]
            report_diff('Dropped {} non-complete genomes.')

        # non-redundant genome IDs
        # typically not necessary, just in case
        self.df.rename(columns={'# assembly_accession': 'accession'},
                       inplace=True)
        self.df['accnov'] = self.df['accession'].str.split('.', 1).str[0]
        self.df['genome'] = 'G' + self.df['accnov'].str.split('_', 1).str[-1]
        self.df.drop_duplicates(subset=['genome'], inplace=True)

        # include/exclude genome Ids
        if self.genoids:
            self.genoids = set(list_from_param(self.genoids))
            print(f'{"Ex" if self.exclude else "In"}cluding '
                  f'{len(self.genoids)} custom genome IDs...')
            self.df = self.df[(self.df['accession'].isin(self.genoids) |
                               self.df['accnov'].isin(self.genoids) |
                               self.df['genome'].isin(self.genoids)) !=
                              self.exclude]
            report_diff('Dropped {} genomes.')

        # genomes without download link
        self.df.query('ftp_path != "na"', inplace=True)
        report_diff('Dropped {} genomes without download link.')
        print('Done.')

    def sort_genomes(self):
        """Sort genomes by quality as informed by metadata.
        """
        # sort by reference > representative > type material > other
        self.df['rc_seq'] = self.df.apply(
            lambda x: 0 if x['refseq_category'] == 'reference genome'
            else (1 if x['refseq_category'] == 'representative genome'
                  else (2 if pd.notnull(x['relation_to_type_material'])
                  else 3)), axis=1)

        # sort by complete > scaffold > contig
        self.df['al_seq'] = self.df['assembly_level'].map(
            {'Chromosome': 0, 'Complete Genome': 0, 'Scaffold': 1,
             'Contig': 2})

        # sort genomes by three criteria
        self.df.sort_values(by=['rc_seq', 'al_seq', 'genome'], inplace=True)

        # clean up temporary columns
        self.df.drop(columns=['al_seq', 'rc_seq'], inplace=True)

    def filter_by_taxonomy(self):
        """Filter genomes by taxonomy.
        """
        print('Filtering genomes by taxonomy...')
        n = self.df.shape[0]

        def report_diff(msg):
            nonlocal n
            n_ = self.df.shape[0]
            if n_ < n:
                print('  ' + msg.format(n - n_))
            n = n_

        # remove non-capitalized organism names
        if self.capital:
            self.df = self.df[self.df['organism_name'].apply(is_capital)]
            report_diff('Dropped {} genomes without captalized organism name.')

        # block certain words in organism names
        if self.block:
            self.block = list_from_param(self.block)
            self.df = self.df[~self.df['organism_name'].apply(
                contain_words, args=(self.block,))]
            report_diff('Dropped {} genomes with one or more blocked words in '
                        'organism name.')

        # remove original species information
        self.df.drop(columns=['species_taxid'], inplace=True)

        # drop genomes whose taxIds are not in taxdump
        self.df.dropna(subset=['taxid'], inplace=True)
        self.df['taxid'] = self.df['taxid'].astype(str)
        self.df = self.df[self.df['taxid'].isin(self.taxdump)]
        report_diff('Dropped {} genomes without valid taxId.')

        # assign genomes to species (represented by taxId not name)
        self.df['species'] = self.df['taxid'].apply(
            taxid_at_rank, rank='species', taxdump=self.taxdump)

        # drop genomes without species taxId
        self.df.dropna(subset=['species'], inplace=True)
        report_diff('Dropped {} genomes without valid species taxId.')

        # drop genomes without Latinate species name
        self.df['latin'] = self.df['species'].apply(
            lambda x: is_latin(self.taxdump[x]['name']))
        if self.latin:
            self.df.query('latin', inplace=True)
            report_diff('Dropped {} genomes without Latinate species name.')
        print('Done.')

        # include/exclude taxIds
        if self.taxids:
            self.taxids = set(list_from_param(self.taxids))
            print(f'{"Ex" if self.exclude else "In"}cluding '
                  f'{len(self.taxids)} custom TaxIDs...')

            self.df = self.df[self.df['taxid'].apply(
                lambda x: is_ancestral(x, self.taxids, self.taxdump)
                != self.exclude)]
            report_diff('Dropped {} genomes.')

    def sample_by_taxonomy(self):
        """Sample genomes at designated taxonomic rank.
        """
        self.selected, n = set(), 0
        rank, sample, latin = self.rank, self.sample, False
        if sample is None:
            self.selected = set(self.df['genome'])
            return
        print('Sampling genomes based on taxonomy...')

        # Latinate species names
        if rank == 'species_latin':
            rank, latin = 'species', True
        print(f'Up to {sample} genome(s) will be sampled per {rank}' +
              (' (Latinate names only) ' if latin else '') + '.')

        # assign genomes to given rank
        if rank not in self.df.columns:
            self.df[rank] = self.df['taxid'].apply(
                taxid_at_rank, rank=rank, taxdump=self.taxdump)
        df_ = self.df.dropna(subset=[rank])

        # keep Latinate species names only
        if latin:
            df_ = df_.query('latin')

        # select up to given number of genomes of each taxonomic group
        self.selected = set(df_.groupby(rank).head(sample)['genome'])
        n = df_[rank].unique().shape[0]
        print(f'  Sampled {len(self.selected)} genomes from {n} '
              f'{rank_plural(rank)}.')

        # sample at ranks above the current one
        if not self.above:
            return

        # determine ranks to sample at
        ranks = ['phylum', 'class', 'order', 'family', 'genus', 'species']
        if rank not in ranks:
            raise ValueError(
                f'Cannot determine taxonomic ranks above "{rank}", because it '
                'is not among the six standard ranks: {", ".join(ranks)}.')
        ranks = ranks[:ranks.index(rank)][::-1]
        print(f'Sampling will also be performed on: {", ".join(ranks)}.')

        for r in ranks:
            self.df['tmp_rank'] = self.df['taxid'].apply(
                taxid_at_rank, rank=r, taxdump=self.taxdump)
            df_ = self.df.dropna(subset=['tmp_rank'])

            # calculate number of genomes yet to be sampled per group
            goal = pd.Series(sample, df_['tmp_rank'].unique())
            done = df_.query('genome in @self.selected')[
                'tmp_rank'].value_counts()
            left = goal.subtract(done, fill_value=0)
            left = left[left > 0].astype(int)

            # sample specific number of genomes per group
            s = []
            for tid, num in left.items():
                s.extend(df_.query(f'tmp_rank == "{tid}"').head(num)['genome'])
            print(f'  Sampled {len(s)} more genomes from {left.shape[0]} '
                  f'{rank_plural(r)}.')
            self.selected.update(s)

        print(f'  Sampled a total of {len(self.selected)} genomes at '
              f'{rank} and above.')
        self.df.drop(columns=['tmp_rank'], inplace=True)

    def sample_by_quality(self):
        """Sample genomes according to quality categories.
        """
        # add reference genomes
        if self.reference:
            df_ = self.df.query('refseq_category == "reference genome"')
            print(f'Include {df_.shape[0]} reference genomes.')
            self.selected.update(df_['genome'])

        # add representative genomes
        if self.represent:
            df_ = self.df.query('refseq_category == "representative genome"')
            print(f'Include {df_.shape[0]} representative genomes.')
            self.selected.update(df_['genome'])

        # add type material genomes
        if self.typemater:
            df_ = self.df[self.df['relation_to_type_material'].notna()]
            print(f'Include {df_.shape[0]} type material genomes.')
            self.selected.update(df_['genome'])

    def filter_to_sampled(self):
        """Filter genomes to sampled ones.
        """
        n = len(self.selected)
        if n == 0:
            raise ValueError('No genome is retained after sampling.')
        self.df.query('genome in @self.selected', inplace=True)
        self.df.sort_values('genome', inplace=True)
        print(f'Total number of sampled genomes: {n}.')

    def download_genomes(self):
        """Download genomes from NCBI.
        """
        if self.manual:
            fname = 'urls.txt'
            self.df['ftp_path'].to_csv(
                join(self.output, 'urls.txt'), header=False, index=False)
            print(f'URLs of sampled genomes written to {fname}. You may '
                  'manually download them later.')
            sys.exit(0)

        print('Downloading non-redundant genomic data from NCBI...',
              flush=True)

        failed = []
        filext = self.filext.split(',')
        ldirs = {}
        for ext in filext:
            ldir = join(self.down, ext.replace('.', '_'))
            makedirs(ldir, exist_ok=True)
            ldirs[ext] = ldir

        for row in self.df.itertuples():
            g = row.genome
            rdir = row.ftp_path.split('/', 3)[-1]
            stem = rdir.rsplit('/', 1)[-1]
            for ext in filext:
                ldir = ldirs[ext]
                fname = f'{stem}_{ext}.gz'
                lfile = join(ldir, fname)
                if check_local_file(lfile):
                    continue
                for i in range(self.retries):
                    ec, _ = run_command(
                        f'rsync -Ltq {server}/{rdir}/{fname} {ldir}')
                    if ec == 0:
                        print('  ' + g, flush=True)
                        break
                    else:
                        sleep(self.delay)
                if ec > 0:
                    print(f'WARNING: Cannot retrieve {fname}.')
                    failed.append(g)
                    break
        print('Done.')

        # drop genomes that cannot be retrieved
        if len(failed):
            print('Failed to retrieve the following genomes:')
            print('  ' + ', '.join(failed))
            failed = set(failed)
            self.df.query('genome in @failures', inplace=True)

    def extract_genomes(self):
        """Extract genome sequences.

        Notes
        -----
        Write all nucleotide sequences into db.fna.
        Write genome to length map to length.map.
        """
        ldir = join(self.down, 'genomic_fna')
        if not isdir(ldir):
            return
        print('Extracting downloaded genomic data...', end='', flush=True)
        fname = 'db.fna'
        fout = open(join(self.output, fname), 'w')
        fout_write = fout.write

        gap = 'N' * 20  # fill contig gaps with 20 "N"s
        gap_join = gap.join

        # extract individual genome sequences
        nucl2l, g2nucl, g2len = {}, {}, {}
        for row in self.df.itertuples():
            g = row.genome
            stem = row.ftp_path.rsplit('/', 1)[-1]
            lfile = join(ldir, f'{stem}_genomic.fna.gz')

            with gzip.open(lfile, 'rb') as f:
                try:
                    content = f.read().decode().splitlines()
                except (OSError, EOFError, TypeError):
                    print(f' skipping corrupted file {stem}.', end='',
                          flush=True)
                    continue

            nucls, seqs = [], []
            nucls_append, seqs_append = nucls.append, seqs.append
            nucl, seq = None, ''  # current nucleotide and sequence
            for line in content:
                if line.startswith('>'):
                    if nucl and seq:
                        nucls_append(nucl)
                        seqs_append(seq)
                        nucl2l[nucl] = len(seq)
                    nucl = line[1:].split(None, 1)[0]
                    seq = ''
                else:
                    seq += line
            if nucl and seq:
                nucls_append(nucl)
                seqs_append(seq)
                nucl2l[nucl] = len(seq)

            g2nucl[g], g2len[g] = nucls, sum(map(len, seqs))
            fout_write(f'>{g}\n{gap_join(seqs)}\n')

        fout.close()
        print(' done.')
        print(f'Combined genome sequences written to {fname}.')

        self.nucl2l = nucl2l
        self.g2nucl = g2nucl
        self.df['contigs'] = self.df['genome'].map({
            k: len(v) for k, v in g2nucl.items()})
        self.df['length'] = self.df['genome'].map(g2len)
        self.df.dropna(subset=['length'], inplace=True)

        # write genome-to-length map
        fname = 'length.map'
        self.df[['genome', 'length']].to_csv(
            join(self.output, fname), sep='\t', header=False, index=False)
        print(f'Genome-to-length map written to {fname}.')

    def genome_lineages(self):
        """Generate lineage information for genomes.

        Notes
        -----
        Write genome lineages to lineages.txt.
        """
        # write genome-to-taxId map
        fname = 'taxid.map'
        self.df[['genome', 'taxid']].to_csv(
            join(self.output, fname), sep='\t', header=False, index=False)
        print(f'Genome-to-taxonomy ID map written to {fname}.')

        # identify taxa at standard ranks
        ranks = ['superkingdom', 'kingdom', 'phylum', 'class', 'order',
                 'family', 'genus', 'species']
        self.df[ranks[:-1]] = self.df['taxid'].apply(lambda x: pd.Series(
            taxids_at_ranks(x, ranks[:-1], self.taxdump)))

        # report number of taxa represented at each rank
        print('Number of taxonomic groups represented:')
        for rank in ranks:
            print(f'  {rank}: {self.df[rank].nunique()}.')

        # merge superkingdom and kingdom
        self.df['kingdom'] = self.df[['superkingdom', 'kingdom']].apply(
            lambda x: x[1] if x[1] else x[0], axis=1)

        # generate lineage string
        self.df['lineage'] = self.df[ranks[1:]].fillna('').apply(
            lambda col: col.apply(lambda val, name: '{}__{}'.format(
                name[0], self.taxdump[val]['name'] if val else ''),
                args=(col.name,))).apply('; '.join, axis=1)

        # write table
        fname = 'lineages.txt'
        self.df[['genome', 'lineage']].to_csv(
            join(self.output, fname), sep='\t', header=False, index=False)
        print(f'Genome lineages written to {fname}.')

    def build_taxdump(self):
        """Build taxonomy database.

        Notes
        -----
        Write refined taxonomy database to taxdump/nodes.dmp and names.dmp.
        """
        taxdump = self.taxdump
        tids = self.df['taxid'].unique().tolist()

        # identify all ancestral taxIds
        ancs = set(tids)
        ancs_add = ancs.add
        for tid in tids:
            cid = tid
            while True:
                pid = taxdump[cid]['parent']
                if cid == pid or pid in ancs:
                    break
                ancs_add(pid)
                cid = pid

        # shrink taxdump files
        fnames = ('nodes.dmp', 'names.dmp')
        for fname in fnames:
            fo = open(join(self.output, fname), 'w')
            fi = open(join(self.tmpdir, fname), 'r')
            for line in fi:
                row = line.rstrip('\r\n').replace('\t|', '').split('\t')
                if row[0] in ancs:
                    if fname == 'nodes.dmp' or 'scientific name' in row:
                        fo.write(line)
            fi.close()
            fo.close()
        print('Taxonomy database written to {}.'.format(' and '.join(fnames)))

    def extract_genes(self):
        """Build taxonomy database.

        Notes
        -----
        Write gene coordinates to coords.txt.
        Write mappings to gene.map, protein.map, product.map.
        """
        ldir = join(self.down, 'genomic_gff')
        if not isdir(ldir):
            return
        print('Extracting genome annotations...', end='', flush=True)

        nucl2l = self.nucl2l
        g2nucl = self.g2nucl

        # read genome annotations
        # write gene coordinates
        fname = 'coords.txt'
        fo = open(join(self.output, fname), 'w')
        fo_write = fo.write
        data = []
        g2n, g2aa = {}, {}

        for row in self.df.itertuples():
            g = row.genome
            fo_write(f'>{g}\n')
            stem = row.ftp_path.rsplit('/', 1)[-1]
            lfile = join(ldir, f'{stem}_genomic.gff.gz')

            with gzip.open(lfile, 'rb') as f:
                try:
                    content = f.read().decode().splitlines()
                except (OSError, EOFError, TypeError):
                    print(f' skipping corrupted file {stem}.', end='',
                          flush=True)
                    continue

            # calculate nucleotide positions in concatenated genome
            nucl2pos, pos = {}, 0
            for nucl in g2nucl[g]:
                nucl2pos[nucl] = pos
                pos += nucl2l[nucl] + 20

            # parse GFF file
            genes = parse_gff(content)

            # calculate gene coordinates in concatenated genome
            nucl, pos, aa = None, 0, 0
            g2n[g] = len(genes)
            for gene in genes:
                if gene[0] != nucl:
                    nucl = gene[0]
                    pos = nucl2pos[nucl]
                beg, end = int(gene[1]), int(gene[2])
                fo_write(f'{gene[4]}\t{beg + pos}\t{end + pos}\n')
                aa += end - beg + 1

            g2aa[g] = aa
            data.extend([[g] + x for x in genes])
        fo.close()
        self.df['proteins'] = self.df['genome'].map(g2n)
        self.df['protein_length'] = self.df['genome'].map(g2aa)
        print(' done.')
        print(f'Gene coordinates written to {fname}.')

        # convert to data frame
        cols = ['genome', 'chromosome', 'start', 'end', 'strand', 'locus_tag',
                'protein_id', 'gene', 'product']
        df_ = pd.DataFrame(data, columns=cols)

        # write locus tag-to-gene symbol map
        fname = 'gene.map'
        df_[['locus_tag', 'gene']].dropna().to_csv(
            join(self.output, fname), sep='\t', header=False, index=False)
        print(f'Locus-to-gene map written to {fname}.')

        # write locus tag-to-protein Id map
        fname = 'protein.map'
        df_[['locus_tag', 'protein_id']].dropna().to_csv(
            join(self.output, fname), sep='\t', header=False, index=False)
        print(f'Locus-to-protein map written to {fname}.')

        # write locus tag-to-protein Id map
        fname = 'product.map'
        df_[['locus_tag', 'product']].dropna().to_csv(
            join(self.output, fname), sep='\t', header=False, index=False)
        print(f'Locus-to-product map written to {fname}.')

    def genome_metadata(self):
        """Write and report genome metadata.

        Notes
        -----
        Write genome metadata to metadata.tsv.
        """
        cols = ['asm_name', 'accession', 'bioproject', 'biosample',
                'assembly_level', 'taxid', 'organism_name',
                'infraspecific_name', 'isolate', 'ftp_path']
        if 'proteins' in self.df.columns:
            cols = ['proteins', 'protein_length'] + cols
        if 'contigs' in self.df.columns:
            cols = ['contigs', 'length'] + cols

        fname = 'metadata.tsv'
        self.df.set_index('genome', inplace=True)
        self.df[cols].to_csv(join(self.output, fname), sep='\t')
        print(f'Genome metadata written to {fname}.')


def check_local_file(file, overwrite=False):
    """Check existing local file.

    Parameters
    ----------
    file : str
        expected local file path and name
    overwrite : bool, optional
        whether overwrite existing file (default: False)

    Returns
    -------
    bool
        whether file exists
    """
    if isfile(file):
        if getsize(file) == 0:
            remove(file)
        elif overwrite:
            print(f'  WARNING: existing local file {basename(file)} will '
                  'be overwritten.')
            remove(file)
        else:
            print(f'  Using local file {basename(file)}.')
            return True
    return False


def run_command(cmd, capture=True, merge=True):
    """Run an external command and retrieve screen output.

    Parameters
    ----------
    cmd : str
        command to execute
    capture : bool, optional
        whether capture screen output (stdout) (default: True)
    merge : bool, optional
        whether merge stderr into stdout (default: True)

    Returns
    -------
    int
        exit code
    list or None
        screen output split by line, or None if not capture
    """
    res = sp.run(cmd, shell=True,
                 stdout=(sp.PIPE if capture else None),
                 stderr=(sp.STDOUT if merge else sp.DEVNULL))
    return (res.returncode,
            res.stdout.decode('utf-8').splitlines() if capture else None)


def list_from_param(param):
    """Read list of entries from file or string.

    Parameters
    ----------
    param : str
        parameter which may be a file path or a comma-delimited string

    Returns
    -------
    list
        list of entries
    """
    if not param:
        return []
    elif isinstance(param, list):
        return param
    elif isinstance(param, str):
        if isfile(param):
            with open(param, 'r') as f:
                return f.read().splitlines()
        else:
            return param.split(',')


def read_taxdump(dir_):
    """Read NCBI taxdump or compatible taxonomy systems.

    Parameters
    ----------
    dir_ : str
        directory containing nodes.dmp and names.dmp

    Returns
    -------
    dict
        taxonomy database
    """
    taxdump = {}
    with open(join(dir_, 'nodes.dmp'), 'r') as f:
        for line in f:
            x = line.rstrip('\r\n').replace('\t|', '').split('\t')
            taxdump[x[0]] = {'parent': x[1], 'rank': x[2]}
    with open(join(dir_, 'names.dmp'), 'r') as f:
        for line in f:
            x = line.rstrip('\r\n').replace('\t|', '').split('\t')
            if len(x) < 4 or x[3] == 'scientific name':
                try:
                    taxdump[x[0]]['name'] = x[1]
                except KeyError:
                    pass
    return taxdump


def parse_gff(fh):
    """Extract information from a GFF file.

    Parameters
    ----------
    fh : file handle
        GFF file to be parsed

    Returns
    -------
    list of list
        chromosome, start, end, strand, locus tag, protein id, gene, product

    Notes
    -----
    Only "CDS" features will be extracted.
    """
    res = []
    res_append = res.append
    for line in fh:
        if line[0] == '#':
            continue
        row = line.rstrip('\r\n').split('\t')
        if row[2] != 'CDS':
            continue
        cds = [row[0], row[3], row[4], row[6]]
        attr = dict(x.split('=') for x in row[8].split(';'))
        cds.extend([attr[x] if x in attr else None for x in (
            'locus_tag', 'protein_id', 'gene', 'product')])
        res_append(cds)
    return res


def is_latin(name):
    """Check if a species name is Latinate.

    Parameters
    ----------
    name : str
        species name to check

    Returns
    -------
    bool
        whether species name is Latinate
    """
    if name == '':
        return False
    elif name.count(' ') != 1:
        return False
    if name[0] == '[':
        i = name.find(']')
        if i == -1:
            return False
        name = name[1:i] + name[i + 1:]
    name = name.replace(' ', '')
    if not name.istitle():
        return False
    elif not name.isalpha():
        return False
    return True


def is_capital(name):
    """Check if a taxon name is capitalized.

    Parameters
    ----------
    name : str
        taxon name to check

    Returns
    -------
    bool
        whether taxon name is capitalized

    Notes
    -----
    A valid taxon name may start with "[".
    """
    try:
        return name.lstrip('[')[0].isupper()
    except IndexError:
        return False


def contain_words(text, words):
    """Check if a string contains any of given words

    Parameters
    ----------
    text : str
        query string
    words : list of str
        target words

    Returns
    -------
    bool
        whether string contains at least one of given words
    """
    return re.search(r'\b{}\b'.format('|'.join(words)), text,
                     re.IGNORECASE) is not None


def _get_taxon(tid, taxdump):
    """Get information of a given taxId from taxonomy database.

    Parameters
    ----------
    tid : str
        taxId to query
    taxdump : dict
        taxonomy database

    Returns
    -------
    dict
        information of taxon

    Raises
    ------
    ValueError
        If taxId is not found in taxonomy database.
    """
    try:
        return taxdump[tid]
    except KeyError:
        raise ValueError(f'TaxID {tid} is not found in taxonomy database.')


def get_lineage(qid, taxdump):
    """Get taxIds of self and all ancestral hierarchies in order.

    Parameters
    ----------
    qid : str
        query taxId
    taxdump : dict
        taxonomy database

    Returns
    -------
    list of str
        taxIds in order (low-to-hight)
    """
    cid = qid
    pid = ''
    lineage = [qid]
    while True:
        taxon = _get_taxon(cid, taxdump)
        pid = taxon['parent']
        if pid == cid or pid == '0':
            break
        lineage.append(pid)
        cid = pid
    return lineage


def is_ancestral(qid, tids, taxdump):
    """Loop up the taxdump hierarchies for a particular taxId.

    Parameters
    ----------
    qid : str
        query taxId
    tids : set of str
        target taxIds
    taxdump : dict
        taxonomy database

    Returns
    -------
    bool
        whether any target taxId is reached
    """
    cid = qid
    pid = ''
    while True:
        if cid in tids:
            return True
        pid = _get_taxon(cid, taxdump)['parent']
        if pid == cid or pid == '0':
            break
        cid = pid
    return False


def taxid_at_rank(qid, rank, taxdump):
    """Find taxId at certain rank for a query taxId.

    Parameters
    ----------
    qid : str
        query taxId
    rank : str
        target rank
    taxdump : dict
        taxonomy database

    Returns
    -------
    str or None
        taxId at target rank, or None if not found
    """
    cid = qid
    pid = ''
    while True:
        taxon = _get_taxon(cid, taxdump)
        if taxon['rank'] == rank:
            return cid
        pid = taxon['parent']
        if pid == cid or pid == '0':
            break
        cid = pid
    return None


def taxids_at_ranks(qid, ranks, taxdump):
    """Find taxId at certain rank for a query taxId.

    Parameters
    ----------
    qid : str
        query taxId
    rank : list of str
        target ranks
    taxdump : dict
        taxonomy database

    Returns
    -------
    dict of str or None
        taxIds at target ranks, or None if not found
    """
    cid = qid
    pid = ''
    res = {x: None for x in ranks}
    rankset = set(ranks)
    while True:
        taxon = _get_taxon(cid, taxdump)
        rank = taxon['rank']
        if rank in rankset:
            res[rank] = cid
        pid = taxon['parent']
        if pid == cid or pid == '0':
            break
        cid = pid
    return res


def find_lca(tids, taxdump):
    """Find the lowest common ancestor (LCA) of given taxIds.

    Parameters
    ----------
    tids : iterable of str
        query taxIds
    taxdump : dict
        taxonomy database

    Returns
    -------
    str
        taxId of LCA
    """
    cal = None  # common ancestral lineage (CAL)
    for tid in tids:
        lineage = get_lineage(tid, taxdump)[::-1]
        n = len(lineage)

        # let current lineage be CAL
        if cal is None:
            cal = lineage
            continue

        # attempt to find LCA between current CAL and current lineage
        idx = None
        for i, tid_ in enumerate(cal):
            if i >= n or tid_ != lineage[i]:
                break
            idx = i + 1

        # LCA not found (i.e., even roots are different)
        if idx is None:
            raise ValueError('Cannot find LCA of TaxIDs in database.')

        # reduce CAL
        if idx < len(cal):
            cal = cal[:idx]

    # LCA is the lowest taxId of CAL
    return cal[-1]


def rank_plural(rank):
    """Convert taxonomic rank to plural form.

    Parameters
    ----------
    rank : str
        taxonomic rank

    Returns
    -------
    str
        taxonomic rank in plural form
    """
    if rank.endswith('phylum'):
        return rank[:-2] + 'a'
    elif rank.endswith('class'):
        return rank + 'es'
    elif rank.endswith('family'):
        return rank[:-1] + 'ies'
    elif rank.endswith('genus'):
        return rank[:-2] + 'era'
    elif rank.endswith('species'):
        return rank
    else:
        return rank + 's'


if __name__ == "__main__":
    main()
