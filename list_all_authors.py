#!/usr/bin/env python3

"""
PubMed author list generator

Description:
    Generates a list of authors with the PMID and date of the last co-authored
    paper of each author as well as the number of papers co-authored from an
    Medline-formatted PubMed paper list.
    Useful for batch-populating the NSF collaborators and other affiliations
    (COA) table.

Usage:
    Navigate to: NCBI - My Bibliography - select all - Manage citations - Export
    file (MEDLINE).
    Then execute: python me.py medline.txt > table.tsv

Notes:
    The PMID, date and affiliation correspond to the last co-authored paper as
    determined by the data when it was recorded in PubMed.
    When multiple affiliations are present for one author in one paper, only the
    first affiliation is retained.
    Human review is recommended after the author list is generated.
"""

import fileinput

__author__ = 'Qiyun Zhu'
__license__ = 'BSD-3-Clause'
__version__ = '0.0.1-dev'
__email__ = 'qiyunzhu@gmail.com'


def main():
    """Main workflow.
    """
    # author data
    data = {}

    # current record
    pmid = None
    date = None
    authors = []  # author - department pairs

    # parse individual Medline
    def parse_medline(key, value):
        nonlocal pmid
        nonlocal date
        if key == 'PMID':
            pmid = value
        elif key == 'FAU':  # full author name
            authors.append([value, None])
        elif key == 'AD':  # author department
            if authors[-1][1] is None:
                authors[-1][1] = value
        elif key == 'PHST' and value.endswith('[pubmed]'):
            date = value.split()[0]

    # add one publication record to author data
    def add_record():
        nonlocal pmid
        nonlocal date
        for author, dept in authors:

            # create new record
            if author not in data:
                data[author] = {
                    'dept': dept,
                    'count': 1,
                    'date': date,
                    'pmid': pmid}

            # update existing record
            else:
                rec = data[author]
                rec['count'] += 1
                if date > rec['date']:
                    rec['dept'] = dept
                    rec['date'] = date
                    rec['pmid'] = pmid

        # reset
        pmid = None
        date = None
        authors.clear()

    # cache
    key, value = None, None

    # parse Medline file
    for line in fileinput.input():
        line = line.rstrip()

        # skip empty line (only occurs between records)
        if not line:
            parse_medline(key, value)
            key, value = None, None
            add_record()

        # append continuous line
        elif line.startswith(' '):
            value += ' ' + line.lstrip()

        # extract key / value from new line
        else:
            parse_medline(key, value)
            key, _, value = line.partition('-')
            if not value:
                raise ValueError(f'Invalid Medline: {line}.')
            key = key.rstrip()
            value = value.lstrip()

    # write author list
    print('Author', 'Department', 'Count', 'Date', 'PMID', sep='\t')
    for author, rec in sorted(data.items()):
        print(author, rec['dept'], str(rec['count']), rec['date'], rec['pmid'],
              sep='\t')


if __name__ == "__main__":
    main()
