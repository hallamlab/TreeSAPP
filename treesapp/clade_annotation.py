class CladeAnnotation:
    def __init__(self, name: str, key: str):
        """A class to store clade-level annotations for a reference package"""
        self.name = name  # The sub-group feature annotation
        self.feature = key  # The group annotation name e.g. Function, Activity, Substrate
        self.members = set()  # List of reference package leaf nodes
        self.taxa = set()  # List of representative taxa
        self.colour = ""

    def __str__(self) -> str:
        return self.summarise()

    def summarise(self) -> str:
        summary_str = "Annotation '{}' of feature '{}' status:\n".format(self.name, self.feature)
        if len(self.members) == 0:
            summary_str += "\tUninitiated.\n"
        else:
            summary_str += "\t{} leaf nodes.\n".format(len(self.members))
            summary_str += "\t{} different taxa covered.\n".format(len(set(self.taxa)))
        return summary_str

    def get_internal_nodes(self, internal_node_leaf_map: dict) -> list:
        """Returns a list of all leaf nodes that are represented by only the leaves in self.members."""
        internal_nodes = []
        for i_node, leaves in internal_node_leaf_map.items():
            if self.members.issuperset(set(leaves)):
                internal_nodes.append(i_node)

        return internal_nodes

