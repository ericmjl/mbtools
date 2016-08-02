import networkx as nx
from .primers import gibson_assembly_primers, gibson_sequencing_primers


class GibsonAssembler(nx.DiGraph):
    def __init__(self, sequences):
        super(GibsonAssembler, self).__init__()

        self.sequences = sequences
        self.add_nodes()
        self.add_edges()
        self.compute_assembly_primers()
        self.compute_sequencing_primers()

    def add_nodes(self):
        for s in self.sequences:
            self.add_node(s)

    def add_edges(self):
        for i in range(len(self.sequences)):
            if i < len(self.sequences) - 1:
                self.add_edge(self.sequences[i], self.sequences[i+1])
            else:
                self.add_edge(self.sequences[i], self.sequences[0])

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
