class CladeAnnotation:
    def __init__(self, name: str, key: str):
        """
        A class to store clade-level annotations for a reference package.
        self.members attribute is a dictionary of reference package leaf nodes as keys,
        and the hierarchy depth of taxa they were annotated with as values.
        self.taxa is a list of representative taxa for this CladeAnnotation.

        :param name: The sub-group feature annotation
        :param key: The group annotation name e.g. Function, Activity, Substrate
        """
        self.name = name
        self.feature = key
        self.members = dict()
        self.taxa = set()
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
            if set(self.members).issuperset(set(leaves)):
                internal_nodes.append(i_node)

        return internal_nodes


def generate_leaf_node_memberships(clade_annotations: list) -> dict:
    leaf_membership = {}
    for clade_anno in clade_annotations:  # type: CladeAnnotation
        for leaf_node in clade_anno.members:
            if leaf_node not in leaf_membership:
                leaf_membership[leaf_node] = []
            leaf_membership[leaf_node].append(clade_anno.name)
    return leaf_membership
