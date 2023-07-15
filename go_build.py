#!/usr/bin/env python3
"""Extract hierarchical mappings and descriptions from the Gene Ontology (GO)
database.

Usage:
    python me.py go.obo

Notes:
    The data file "go.obo" can be downloaded from the GO website.
    Tested and working with GO releases 2019-06-09 and 2021-02-01.

Reference:
    The Gene Ontology Consortium, 2021
    https://academic.oup.com/nar/article/49/D1/D325/6027811
"""

import sys


__author__ = 'Qiyun Zhu'
__license__ = 'BSD-3-Clause'
__version__ = '0.0.1-dev'
__email__ = 'qiyunzhu@gmail.com'


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)

    # mappings to generate
    unimaps = {x: {} for x in ('name', 'domain', 'altid', 'replaced')}
    polymaps = {x: {} for x in ('slim', 'parent', 'consider', 'ec', 'kegg',
                                'metacyc', 'reactome', 'rhea')}

    # obsolete terms
    obsoletes = set()

    # helper to add a unique mapping
    def add_unival(key, go, val):
        if go is None:
            return
        if go in unimaps[key]:
            raise ValueError(f'{key} - {go}: {unimaps[key][go]} vs {val}')
        unimaps[key][go] = val

    # helper to add a target to a multiple mapping
    def add_polyval(key, go, val):
        if go is None:
            return
        polymaps[key].setdefault(go, set()).add(val)

    # parsing
    # the file "go.obo" is one of the ontology data release formats
    # the latest release can be downloaded from:
    #   http://purl.obolibrary.org/obo/go.obo
    # the data release is described at:
    #   http://purl.obolibrary.org/obo/go.obo
    # the obo file format is defined at:
    #   http://owlcollab.github.io/oboformat/doc/obo-syntax.html
    f = open(sys.argv[1], 'r')
    go = None
    for line in f:
        line = line.rstrip()

        # start of a GO term
        if not line or line == '[Term]':
            go = None

        # GO term ("GO:#######")
        elif line.startswith('id: GO:'):
            go = line[4:]

        # descriptive name
        elif line.startswith('name: '):
            add_unival('name', go, line[6:])

        # alternative id
        elif line.startswith('alt_id: '):
            add_unival('altid', line[8:], go)

        # domain
        # three domains (namespaces) serve as top-level categories:
        #   molecular function, cellular component, biological process
        # see: http://geneontology.org/docs/ontology-documentation/
        elif line.startswith('namespace: '):
            add_unival('domain', go, line[11:])

        # GO subset (slim)
        elif line.startswith('subset: goslim_'):
            add_polyval('slim', go, line[15:])

        # relation
        # the most basic child-to-parent relation is "is_a"
        # the follow relations trail "relationship" or "intersection_of":
        #   ends_during, happens_during, has_part, negatively_regulates,
        #   occurs_in, part_of, positively_regulates, regulates
        # the official website explains four relations:
        #   is_a, part_of, has_part, regulates
        # see: http://geneontology.org/docs/ontology-relations/
        # the QuickGO website considers the following relations as child-to-
        # parent:
        #   is_a, part_of, regulates, occurs_in, positively_regulates,
        #   negatively_regulates, capable_of, capable_of_part_of
        # example: https://www.ebi.ac.uk/QuickGO/term/GO:0001591
        elif line.startswith('is_a: '):
            add_polyval('parent', go, line[6:].split(' ! ')[0])
        elif line.startswith('relationship: part_of '):
            add_polyval('parent', go, line[22:].split(' ! ')[0])
        elif line.startswith('intersection_of: part_of '):
            add_polyval('parent', go, line[25:].split(' ! ')[0])

        # external reference
        elif line.startswith('xref: '):
            xref = line[6:]
            if xref.startswith('EC:'):
                add_polyval('ec', go, xref[3:].split()[0])
            elif xref.startswith('KEGG_REACTION:'):
                add_polyval('kegg', go, xref[14:].split()[0])
            elif xref.startswith('MetaCyc:'):
                add_polyval('metacyc', go, xref[8:].split()[0])
            elif xref.startswith('Reactome:'):
                add_polyval('reactome', go, xref[9:].split()[0])
            elif xref.startswith('RHEA:'):
                add_polyval('rhea', go, xref[5:].split()[0])

        # obsolete
        # redirected to "replaced_by" and "consider"
        elif line == 'is_obsolete: true':
            obsoletes.add(go)
        elif line.startswith('replaced_by: '):
            if go not in obsoletes:
                raise ValueError(f'Term to replace: {go} is not obsolete.')
            add_unival('replaced', go, line[13:])
        elif line.startswith('consider: '):
            if go not in obsoletes:
                raise ValueError(f'Term to consider: {go} is not obsolete.')
            add_polyval('consider', go, line[10:])

    f.close()

    # write unique mappings
    for key in 'name', 'domain', 'altid', 'replaced':
        with open(f'{key}.map', 'w') as f:
            for go, val in sorted(unimaps[key].items()):
                if go not in obsoletes or key in ('altid', 'replaced'):
                    print(go, val, sep='\t', file=f)

    # write multiple mappings:
    for key in ('slim', 'parent', 'ec', 'kegg', 'metacyc', 'reactome', 'rhea',
                'consider'):
        with open(f'{key}.map', 'w') as f:
            for go, vals in sorted(polymaps[key].items()):
                if go not in obsoletes or key == 'consider':
                    print(go, '\t'.join(sorted(vals)), sep='\t', file=f)

    # obsolete
    with open('obsolete.txt', 'w') as f:
        for go in sorted(obsoletes):
            print(go, file=f)

    # graph
    # the go graph is a directed acyclic graph (DAG), with one child being
    # mapped to one or multiple parents, while all terms eventually grouped
    # into the three domains as roots
    # see: http://geneontology.org/docs/ontology-documentation/
    # this script outputs an DAG in a simple format:
    #   child <tab> parent1 <tab> parent2 <tab> parent3 ...
    graph = {}
    for go, parents in polymaps['parent'].items():
        if go not in obsoletes:
            graph[go] = set(x for x in parents if x not in obsoletes)

    # root
    # perform post-order traversal of the DAG until reaching root(s)
    roots = {}

    def find_root(node):
        if node in roots:
            return roots[node]
        if node in graph:
            root = set().union(*map(find_root, graph[node]))
        else:
            root = {node}
        roots[node] = root
        return root

    for node in graph:
        find_root(node)

    # the root has to be one of the three domains:
    #   GO:0003674 molecular_function
    #   GO:0005575 cellular_component
    #   GO:0008150 biological_process
    domains = ('GO:0003674', 'GO:0005575', 'GO:0008150')

    with open('root.map', 'w') as f:
        for node, root in sorted(roots.items()):
            if root.difference(domains):
                raise ValueError(f'Invalid root for {node}: {root}')
            print(node, '\t'.join(sorted(root)), sep='\t', file=f)

    # slim (subset = generic)
    generic = set(k for k, v in polymaps['slim'].items()
                  if k not in obsoletes and 'generic' in v)
    slims = {}

    def find_slim(node):
        if node in slims:
            return slims[node]
        slim = set()
        if node in generic:
            slim.add(node)
        if node in graph:
            slim.update(set().union(*map(find_slim, graph[node])))
        slims[node] = slim
        return slim

    for node in graph:
        find_slim(node)

    with open('generic.map', 'w') as f:
        for node, slim in sorted(slims.items()):
            if slim:
                print(node, '\t'.join(sorted(slim)), sep='\t', file=f)

    # graph per domain
    # in each graph, all nodes and their parents are under the current domain,
    # all nodes can be traced back to the current domain, if there is a break
    # in the path, discard the node
    for domain in domains:
        name = unimaps['name'][domain]
        graph_ = {}
        for go, parents in graph.items():
            if unimaps['domain'][go] == name:
                graph_[go] = set(x for x in parents if unimaps[
                    'domain'][x] == name)

        # check root
        roots_ = {}

        def find_root_(node):
            if node in roots_:
                return roots_[node]
            if node in graph_:
                root = set().union(*map(find_root_, graph_[node]))
            else:
                root = {node}
            roots_[node] = root
            return root

        for node in graph_:
            find_root_(node)

        # keep continuous paths only
        valid = set()
        for node, root in roots_.items():
            if root == {domain}:
                valid.add(node)

        with open('{}.parent.map'.format(name.split('_')[-1]), 'w') as f:
            for go, parents in sorted(graph_.items()):
                if go not in valid:
                    continue
                parents = parents.intersection(valid)
                if not parents:
                    continue
                print(go, '\t'.join(sorted(parents)), sep='\t', file=f)

        # slim (generic)
        valid = set(k for k, v in unimaps['domain'].items()
                    if k not in obsoletes and v == name)
        with open('{}.generic.map'.format(name.split('_')[-1]), 'w') as f:
            for node, slim in sorted(slims.items()):
                if node not in valid:
                    continue
                slim = slim.intersection(valid) - {domain}
                if slim:
                    print(node, '\t'.join(sorted(slim)), sep='\t', file=f)


if __name__ == '__main__':
    main()
