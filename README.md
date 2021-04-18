# My Utilities

This repository hosts a collection of single-file scripts I wrote that may be useful to people. They are typically written in Python 3, with a simple command-line interface including help information.


## [list_all_authors.py](list_all_authors.py)

Generates a list of authors with their information on the last paper of each author as well as the number of papers co-authored from a Medline-formatted PubMed paper list. Useful for populating the collaborators and other affiliations (COA) table for grant agencies.


## [kegg_query.py](kegg_query.py)

Fetch hierarchical mappings and descriptions of a given list of KEGG
Orthology (KO) entries from the [KEGG server](https://www.kegg.jp/kegg/) using the official [KEGG API](https://www.kegg.jp/kegg/rest/).


## [metacyc_build.py](metacyc_build.py)

Extract hierarchical mappings and descriptions from a local copy of the [MetaCyc](https://metacyc.org/) database.


## [go_build.py](go_build.py)

Extract hierarchical mappings and descriptions from the [Gene Ontology](http://geneontology.org/) (GO) database.
