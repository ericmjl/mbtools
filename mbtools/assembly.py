import networkx as nx
from .primers import gibson_assembly_primers, gibson_sequencing_primers
from .pcr import compute_amplicon
from Bio.SeqRecord import SeqRecord
from collections import defaultdict


class GibsonAssembler(nx.DiGraph):
    """
    Inherits nx.DiGraph.

    Graph definition:
    - nodes: BioPython Seq objects
        - id: seqrecord ID
        - obj: SeqRecord object
    - edges: the order in which parts are to be assembled.
    """
    def __init__(self, sequences):
        super(GibsonAssembler, self).__init__()

        self.add_sequences(sequences)
        self.add_nodes()
        self.add_edges()
        self.compute_assembly_primers()
        self.compute_sequencing_primers()

    def add_sequences(self, sequences):
        self.sequences = []
        for s in sequences:
            assert isinstance(s, SeqRecord),\
                "sequences must be an iterable of BioPython SeqRecord objects."
            self.sequences.append(s)

    def add_nodes(self):
        for i, s in enumerate(self.sequences):
            self.add_node(s.seq, id=s.id, obj=s)

    def add_edges(self):
        for i, s in enumerate(self.sequences):
            if i < len(self.sequences) - 1:
                self.add_edge(s.seq,
                              self.sequences[i+1].seq)
            else:
                self.add_edge(s.seq,
                              self.sequences[0].seq)

    def compute_assembly_primers(self):
        for up, dn in self.edges():
            dn_fw, up_re = gibson_assembly_primers(up, dn)

            self.node[dn]['fw_gibson'] = dn_fw
            self.node[up]['re_gibson'] = up_re

    def compute_sequencing_primers(self):
        for up, dn in self.edges():
            up_3p, dn_5p = gibson_sequencing_primers(up, dn)
            self.node[up]['3p_sequencing'] = up_3p
            self.node[dn]['5p_sequencing'] = dn_5p

    def primers(self):
        """
        Returns a dictionary of part ID and primers.
        """
        primers = defaultdict(dict)
        for n, d in self.nodes(data=True):
            for p in ['fw_gibson', 're_gibson',
                      '3p_sequencing', '5p_sequencing']:
                primers[d['id']][p] = d[p]
        return primers

    def pcr_products(self):
        """
        Returns a dictionary of part ID and PCR product sequences.
        """
        pcr_products = dict()
        for n, d in self.nodes(data=True):
            pcr_products[d['id']] = compute_amplicon(d['fw_gibson'],
                                                     d['re_gibson'],
                                                     n)

        return pcr_products
