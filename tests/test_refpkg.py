import unittest
import pytest
import os
from shutil import rmtree

from . import testing_utils as utils


@pytest.fixture(scope="class")
def refpkg_class(request):
    from treesapp import refpkg
    request.cls.db = refpkg.ReferencePackage("McrA")
    request.cls.db.f__json = utils.get_test_data(os.path.join("refpkgs", "McrA_build.pkl"))
    request.cls.db.slurp()
    return


@pytest.mark.usefixtures("refpkg_class")
class RefPkgTester(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp.utilities import fetch_executable_path
        self.pkl_path = utils.get_test_data(os.path.join("refpkgs", "McrA_build.pkl"))
        self.new_pkl_path = "./test_write_json" + self.db.refpkg_suffix
        self.disband_path = os.path.join("tests", "_".join([self.db.prefix, self.db.refpkg_code, self.db.date]))
        if os.path.isdir(self.disband_path):
            rmtree(self.disband_path)
        self.intermediates_dir = "refpkg_test_dir"
        if os.path.isdir(self.intermediates_dir):
            rmtree(self.intermediates_dir)
        os.mkdir(self.intermediates_dir)

        # Find the executables
        self.exe_map = {}
        self.treesapp_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.sep
        for dep in ["hmmbuild", "hmmalign", "raxml-ng", "mafft", "BMGE.jar", "FastTree"]:
            self.exe_map[dep] = fetch_executable_path(dep, utils.get_treesapp_root())
        return

    def tearDown(self) -> None:
        if os.path.isdir(self.disband_path):
            rmtree(self.disband_path)
        if os.path.isfile(self.new_pkl_path):
            os.remove(self.new_pkl_path)
        if os.path.isdir(self.intermediates_dir):
            rmtree(self.intermediates_dir)
        return

    def test_band(self):
        self.db.change_file_paths(new_dir=os.path.dirname(self.db.f__json))
        self.db.band()
        return

    def test_disband(self):
        # Basic disband
        self.db.disband(output_dir="./tests/")
        self.assertTrue(os.path.isfile(os.path.join(self.disband_path, "McrA.fa")))

        # Ensure the optional files are written properly too
        self.db.boot_tree = "test_string"
        self.db.disband(output_dir="./tests/")
        self.assertTrue(os.path.isfile(os.path.join(self.disband_path, "McrA_bipart.nwk")))
        return

    def test_remove_taxon_from_lineage_ids(self):
        # Ensure the initial state is as expected
        self.assertEqual(249, self.db.num_seqs)
        self.assertEqual(249, len(self.db.lineage_ids))

        # Test
        self.db.remove_taxon_from_lineage_ids("r__Root; d__Archaea; p__Euryarchaeota; c__Methanobacteria;"
                                              " o__Methanobacteriales; f__Methanobacteriaceae; g__Methanosphaera")
        self.assertEqual(191, self.db.num_seqs)
        self.db.remove_taxon_from_lineage_ids("r__Root; d__Archaea; p__Euryarchaeota; c__Methanobacteria")
        self.assertEqual(144, len(self.db.lineage_ids))

        self.db.remove_taxon_from_lineage_ids("r__Root; d__Archaea")
        self.assertEqual(0, len(self.db.lineage_ids))
        self.assertEqual(0, self.db.num_seqs)

    def test_get_fasta(self):
        refpkg_fa = self.db.get_fasta()
        self.assertEqual(249, len(refpkg_fa.fasta_dict))
        return

    def test_get_ete_tree(self):
        from treesapp.refpkg import ReferencePackage
        rt = self.db.get_ete_tree()
        self.assertEqual(249, len(rt))
        blank = ReferencePackage()
        with pytest.raises(SystemExit):
            blank.get_ete_tree()
        return

    def test_hmm_length(self):
        from treesapp.refpkg import ReferencePackage
        blank = ReferencePackage()
        self.assertEqual(blank.profile, [])
        with pytest.raises(ValueError):
            blank.hmm_length()
        return

    def test_bail(self):
        self.db.bail()

    def test_get_internal_node_leaf_map(self):
        self.db.f__json = self.pkl_path
        self.db.slurp()
        node_map = self.db.get_internal_node_leaf_map()
        self.assertEqual(2*self.db.num_seqs, len(node_map))
        self.assertEqual(self.db.num_seqs, len(node_map[max(node_map.keys())]))
        return

    def test_enumerate_taxonomic_lineages(self):
        from treesapp.refpkg import ReferencePackage
        mock_rp = ReferencePackage()
        self.assertEqual(0, len(mock_rp.enumerate_taxonomic_lineages()))
        taxa_counts = self.db.enumerate_taxonomic_lineages()
        self.assertIsInstance(taxa_counts, dict)
        self.assertTrue("r__Root" in taxa_counts)
        return

    def test_exclude_clade_from_ref_files(self):
        ce_taxon = "c__Methanosarcinales"
        ce_refpkg = self.db.clone(clone_path=os.path.join(self.intermediates_dir, ce_taxon,
                                                          self.db.prefix + self.db.refpkg_suffix))

        # Fail due to invalid target taxon
        with pytest.raises(SystemExit):
            ce_refpkg.exclude_clade_from_ref_files(tmp_dir=os.path.join(self.intermediates_dir, ce_taxon),
                                                   executables=self.exe_map, target_clade=ce_taxon)

        # Work as normal
        ce_taxon = "r__Root; d__Archaea; p__Euryarchaeota; c__Methanomicrobia; o__Methanosarcinales"
        ce_refpkg.exclude_clade_from_ref_files(tmp_dir=os.path.join(self.intermediates_dir, ce_taxon.split("; ")[-1]),
                                               executables=self.exe_map,
                                               target_clade=ce_taxon, fresh=True)
        self.assertEqual(144, ce_refpkg.num_seqs)
        return

    def test_pickle_package(self):
        self.db.f__json = self.new_pkl_path
        self.db.pickle_package()
        self.db.slurp()
        self.assertTrue("McrA" == self.db.prefix)
        return

    def test_slurp(self):
        self.db.f__json = self.pkl_path
        self.db.slurp()
        self.assertEqual("McrA", self.db.prefix)
        self.assertEqual(249, self.db.num_seqs)
        return

    def test_taxonomically_label_tree(self):
        from treesapp.taxonomic_hierarchy import Taxon
        labelled_rt = self.db.taxonomically_label_tree()
        self.assertEqual(249, len(labelled_rt))
        self.assertTrue('taxon' in labelled_rt.features)
        self.assertEqual('Root', labelled_rt.taxon.name)
        self.assertEqual('root', labelled_rt.taxon.rank)
        self.assertTrue(isinstance(labelled_rt.taxon, Taxon))
        return

    def test_rename(self):
        from treesapp.refpkg import rename, ReferencePackage
        self.assertEqual("McrA", self.db.prefix)
        with pytest.raises(SystemExit):
            rename(refpkg=self.db, attributes=['McrAA'], output_dir=self.disband_path, overwrite=False)
        with pytest.raises(SystemExit):
            rename(refpkg=self.db, attributes=['prefix', 'McrAA', 'no'], output_dir=self.disband_path, overwrite=False)

        rename(refpkg=self.db, attributes=['prefix', 'McrA2'], output_dir=self.disband_path, overwrite=False)
        # Load the edited reference package
        renamed_refpkg = ReferencePackage()
        renamed_refpkg.f__json = os.path.join(self.disband_path, "McrA2_build.pkl")
        renamed_refpkg.slurp()
        self.assertEqual("McrA2", renamed_refpkg.prefix)

        # Test remaking in same directory
        wd = os.getcwd()
        os.chdir(self.disband_path)
        rename(renamed_refpkg, attributes=["prefix", "TMP"], output_dir="./", overwrite=True)
        # Move back
        os.chdir(wd)
        return

    def test_view(self):
        from treesapp.refpkg import view
        with pytest.raises(SystemExit):
            view(self.db, ["boot_tree", "f__boot_tree", "annotations"])
        view(self.db, attributes=["num_seqs"])
        return


if __name__ == '__main__':
    unittest.main()
