class CladeAnnotation:
    def __init__(self, name: str, key: str):
        """A class to store clade-level annotations for a reference package"""
        self.name = name  # The sub-group feature annotation
        self.feature = key  # The group annotation name e.g. Function, Activity, Substrate
        self.taxa_leaf_map = {}
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

