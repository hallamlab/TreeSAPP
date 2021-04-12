import unittest
import pytest


class CreateRefPkgTester(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp.entrez_utils import EntrezRecord
        er_one = EntrezRecord("seq_1", "seq_1.1")
        er_two = EntrezRecord("seq_2", "seq_2.1")
        er_three = EntrezRecord("seq_3", "seq_3.1")
        er_four = EntrezRecord("seq_4", "seq_4.1")
        self.entrez_record_dict = {'1': er_one, '2': er_two, '3': er_three, '4': er_four}
        return

    def test_present_cluster_rep_options(self):
        from treesapp.create_refpkg import present_cluster_rep_options
        from treesapp.seq_clustering import Cluster
        from treesapp.fasta import Header
        # Make the mock test data
        clust_one = Cluster("seq_3")
        clust_one.members = ["seq_3", "seq_2"]
        clust_two = Cluster("seq_1")
        clust_two.members = ["seq_1", "seq_4"]
        cluster_dict = {'0': clust_one, '1': clust_two}

        head_one = Header("seq_1")
        head_two = Header("seq_2")
        head_three = Header("seq_3")
        head_four = Header("seq_4")
        header_dict = {'1': head_one, '2': head_two, '3': head_three, '4': head_four}
        # Test failure
        present_cluster_rep_options(cluster_dict, self.entrez_record_dict, header_dict,
                                    default=1, each_lineage=True)

        # Test success, but no
        self.assertEqual("seq_3", clust_one.representative)
        self.assertEqual("seq_1", clust_two.representative)
        return

    def test_lineages_to_dict(self):
        from treesapp.create_refpkg import lineages_to_dict
        from treesapp.entrez_utils import EntrezRecord
        lineages_to_dict(fasta_entrez_records=self.entrez_record_dict, taxa_lca=True)

        tmp_entrez_records = self.entrez_record_dict
        tmp_entrez_records['5'] = EntrezRecord("seq_5", '')
        with pytest.raises(AssertionError):
            lineages_to_dict(tmp_entrez_records)
        return


if __name__ == '__main__':
    unittest.main()
