import unittest
import os
import shutil

from .testing_utils import get_test_data, get_treesapp_root


class CladeExclusionTester(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp.refpkg import ReferencePackage
        from treesapp.utilities import fetch_executable_path

        # Create and slurp the reference package that will be used for clade exclusions
        self.ref_pkg = ReferencePackage()
        self.puha_pkl = get_test_data(os.path.join("refpkgs", "PuhA_build.pkl"))
        self.ref_pkg.f__json = get_test_data(self.puha_pkl)
        self.ref_pkg.slurp()

        # Create the directory for temporary files
        self.tmp_dir = "./ce_intermediates"
        if os.path.isdir(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
        os.mkdir(self.tmp_dir)

        # Find the executables
        self.exe_map = {}
        self.treesapp_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.sep
        for dep in ["hmmbuild", "hmmalign", "raxml-ng", "mafft", "BMGE.jar"]:
            self.exe_map[dep] = fetch_executable_path(dep, get_treesapp_root())

        self.num_procs = 2
        return

    def tearDown(self) -> None:
        if os.path.isdir(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
        return

    def test_get_testable_lineages_for_rank(self):
        from treesapp import clade_exclusion_evaluator
        from treesapp.refpkg import ReferencePackage
        # Generate the inputs
        mcra_rp = ReferencePackage()
        mcra_rp.f__json = get_test_data(get_test_data(os.path.join("refpkgs", "McrA_build.pkl")))
        mcra_rp.slurp()
        rlm = {leaf.number: leaf.lineage for leaf in mcra_rp.generate_tree_leaf_references_from_refpkg()}
        qlm = {'q1': 'r__Root; d__Bacteria',
               'q2': 'r__Root; d__Archaea; p__Euryarchaeota; c__Methanobacteria; o__Methanobacteriales; f__Methanobacteriaceae',
               'q3': 'r__Root; d__Archaea; p__Euryarchaeota; c__Methanobacteria; o__Methanobacteriales; f__Methanobacteriaceae; g__Methanobacterium; s__Methanobacterium alcaliphilum',
               'q4': 'r__Root; d__Archaea; p__Euryarchaeota',
               'q5': 'r__Root; d__Archaea; p__Euryarchaeota; c__Methanococci; o__Methanococcales',
               'q6': 'r__Root; d__Bacteria; p__Proteobacteria; c__Betaproteobacteria'}

        # Test when a lineage isn't in the taxonomic hierarchy
        lineages = clade_exclusion_evaluator.get_testable_lineages_for_rank(ref_lineage_map=rlm, query_lineage_map=qlm,
                                                                            rank="phylum",
                                                                            taxonomy=mcra_rp.taxa_trie)
        self.assertEqual(1, len(lineages))

        # Test when the ref_lineage_map lineages are not rooted
        rlm = {'1': "d__Archaea; p__Euryarchaeota; c__Methanobacteria; o__Methanobacteriales; f__Methanobacteriaceae; g__Methanobacterium; s__Methanobacterium sp. Maddingley MBC34",
               '2': "d__Archaea; p__Euryarchaeota; c__Halobacteria; o__Halobacteriales; f__Haloarculaceae; g__Halapricum; s__Halapricum salinum",
               '3': "d__Archaea; p__Euryarchaeota; c__Methanococci; o__Methanococcales; f__Methanococcaceae; g__Methanococcus; s__Methanococcus aeolicus",
               '4': "d__Archaea; p__Euryarchaeota; c__Methanomicrobia; o__Methanosarcinales; f__Methanosarcinaceae; g__Methanosarcina"}
        clade_exclusion_evaluator.check_lineage_compatibility(rlm, mcra_rp.taxa_trie)
        lineages = clade_exclusion_evaluator.get_testable_lineages_for_rank(ref_lineage_map=rlm, query_lineage_map=qlm,
                                                                            rank="class",
                                                                            taxonomy=mcra_rp.taxa_trie)
        self.assertEqual(2, len(lineages))

        # Test normal conditions
        lineages = clade_exclusion_evaluator.get_testable_lineages_for_rank(ref_lineage_map=rlm, query_lineage_map=qlm,
                                                                            rank="class",
                                                                            taxonomy=mcra_rp.taxa_trie)
        self.assertEqual(["r__Root; d__Archaea; p__Euryarchaeota; c__Methanobacteria",
                          "r__Root; d__Archaea; p__Euryarchaeota; c__Methanococci"], lineages)
        return

    def test_select_rep_seqs(self):
        from treesapp.clade_exclusion_evaluator import select_rep_seqs
        from treesapp.entrez_utils import EntrezRecord
        # Make some fake Entrez records
        er_one = EntrezRecord('CMV57640', 'CMV57640.1')
        er_one.sequence = "AAAA"
        er_two = EntrezRecord('WP_000245505', 'WP_000245505.1')
        er_two.sequence = "CCCC"
        er_three = EntrezRecord("er_three", "er_three.1")
        er_three.sequence = "GGGG"

        mocked_dedup_assignments = {'r__Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales;'
                                    ' f__Streptococcaceae; g__Streptococcus; s__Streptococcus pneumoniae':
                                        ['CMV57640', 'WP_000245505']}
        mocked_entrez_records = {'1': er_one, '2': er_two, '3': er_three}
        selected_fa_dict = select_rep_seqs(mocked_dedup_assignments, mocked_entrez_records, self.ref_pkg.taxa_trie)
        self.assertEqual(2, len(selected_fa_dict))
        self.assertEqual("AAAA", selected_fa_dict['CMV57640'])

        # Ensure the taxonomic filtering works
        mocked_dedup_assignments.update({"r__Root; d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria;"
                                         " o__Caulobacterales; f__Caulobacteraceae; g__Brevundimonas":
                                             ["er_three"]})
        selected_fa_dict = select_rep_seqs(deduplicated_assignments=mocked_dedup_assignments,
                                           test_sequences=mocked_entrez_records,
                                           taxon_hierarchy=self.ref_pkg.taxa_trie,
                                           target_lineage="r__Root; d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria")
        self.assertEqual(1, len(selected_fa_dict))
        self.assertEqual("er_three", list(selected_fa_dict.keys())[0])
        return

    def test_prep_graftm_ref_files(self):
        from treesapp.clade_exclusion_evaluator import prep_graftm_ref_files
        taxon_str = 'd__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales'
        ce_dir = os.path.join(self.tmp_dir, taxon_str.split("; ")[-1])
        prep_graftm_ref_files(tmp_dir=ce_dir, target_clade=taxon_str,
                              ref_pkg=self.ref_pkg, executables=self.exe_map)
        self.assertTrue(os.path.exists(os.path.join(ce_dir, "PuhA_lineage_ids.txt")))
        self.assertTrue(os.path.exists(os.path.join(ce_dir, "PuhA.fa")))
        self.assertTrue(os.path.exists(os.path.join(ce_dir, "PuhA.mfa")))
        return

    # def test_build_graftm_package(self):
    #     from treesapp.clade_exclusion_evaluator import build_graftm_package
    #     pkg_path = os.path.join(self.tmp_dir, "PuhA.gpkg")
    #     build_graftm_package(gpkg_path=pkg_path, threads=self.num_procs,
    #                          tax_file=get_test_data("PuhA_lineage_ids.txt"),
    #                          fa_file=get_test_data("PuhA.fa"),
    #                          mfa_file=get_test_data("PuhA.mfa"))
    #     self.assertTrue(os.path.isdir(pkg_path))
    #     return


if __name__ == '__main__':
    unittest.main()
