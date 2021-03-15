import unittest
import pytest

from .testing_utils import get_test_data


@pytest.fixture(scope="class")
def jp_dat(request):
    from treesapp.jplace_utils import jplace_parser, demultiplex_pqueries
    request.cls.db = jplace_parser(get_test_data("epa_result.jplace"))
    pqueries = demultiplex_pqueries(request.cls.db)
    request.cls.db.pqueries = pqueries
    return


@pytest.mark.usefixtures("jp_dat")
class JPlaceTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.test_jplace = get_test_data("epa_result.jplace")
        self.placement_dict = {'p': [[0, -56570.7911467148, 1.0, 0.859477616, 1.2276349525]],
                               'n': ['scaffold_52752_c1_1 partial=11;start_type=Edge;rbs_motif=None|McrA|1_204']}
        return

    def test_jplace_parser(self):
        from treesapp.jplace_utils import jplace_parser
        jplace_dat = jplace_parser(self.test_jplace)
        self.assertEqual(12, len(jplace_dat.pqueries))
        self.assertIsInstance(jplace_dat.pqueries, list)
        return

    def test_demultiplex_pqueries(self):
        from treesapp.jplace_utils import demultiplex_pqueries, jplace_parser
        from treesapp.phylo_seq import PQuery
        jplace_dat = jplace_parser(self.test_jplace)
        pqueries = demultiplex_pqueries(jplace_dat)
        for pquery in pqueries:
            self.assertIsInstance(pquery.placements, list)
        self.assertEqual(12, len(pqueries))
        self.assertIsInstance(pqueries[0], PQuery)
        return

    def test_write_jplace(self):
        from os import path, remove
        output_jplace = "./tmp.jplace"
        self.db.write_jplace(output_jplace)
        self.assertTrue(path.isfile(output_jplace))
        remove(output_jplace)
        return

    def test_summarize(self):
        from treesapp.jplace_utils import jplace_parser
        summary = jplace_parser(self.test_jplace).summarize()
        self.assertTrue(summary)
        return

    def test_add_bipartitions(self):
        from treesapp import jplace_utils
        from treesapp import refpkg
        import re
        jplace_dat = jplace_utils.jplace_parser(self.test_jplace)
        mcra_refpkg = refpkg.ReferencePackage()
        mcra_refpkg.f__pkl = get_test_data("refpkgs/McrA_build.pkl")
        mcra_refpkg.slurp()
        test_ete_tree = mcra_refpkg.get_ete_tree()

        # Test with wrong format of bipartitions
        with pytest.raises(TypeError):
            jplace_utils.add_bipartitions(jplace_dat, bipartitions={})

        # Test with a properly formatted list
        jplace_utils.add_bipartitions(jplace_dat, mcra_refpkg.boot_tree)
        bootstrap_vals = []
        end_pos = 0
        while jplace_dat.tree:
            jplace_dat.tree = jplace_dat.tree[end_pos:]
            match = re.search(r"\[([0-9.]+)]", jplace_dat.tree)
            if not match:
                break
            bootstrap_vals.append(match.group())
            end_pos = match.end()
        self.assertEqual(len(test_ete_tree.get_descendants()), len(bootstrap_vals))

        return


if __name__ == '__main__':
    unittest.main()
