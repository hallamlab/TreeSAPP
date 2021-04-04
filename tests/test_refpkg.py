import unittest
import pytest
import os
from shutil import rmtree

from . import testing_utils as utils


@pytest.fixture(scope="class")
def refpkg_class(request):
    from treesapp import refpkg
    request.cls.db = refpkg.ReferencePackage("McrA")
    request.cls.db.f__pkl = utils.get_test_data(os.path.join("refpkgs", "McrA_build.pkl"))
    request.cls.db.slurp()
    return


@pytest.mark.usefixtures("refpkg_class")
class RefPkgTester(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp.refpkg import ReferencePackage
        from treesapp.utilities import fetch_executable_path
        from shutil import copyfile
        self.pkl_path = utils.get_test_data(os.path.join("refpkgs", "McrA_build.pkl"))
        self.new_pkl_path = "./test_write_json" + self.db.refpkg_suffix
        self.disband_path = os.path.join("tests", self.db.refpkg_code)
        if os.path.isdir(self.disband_path):
            rmtree(self.disband_path)
        self.intermediates_dir = os.path.join("tests", "refpkg_test_dir") + os.sep
        if os.path.isdir(self.intermediates_dir):
            rmtree(self.intermediates_dir)
        os.mkdir(self.intermediates_dir)

        # Copy the original McrA test reference package to a new directory
        copyfile(self.pkl_path, os.path.join(self.intermediates_dir, os.path.basename(self.db.f__pkl)))
        self.mutable_ref_pkg = ReferencePackage()
        self.mutable_ref_pkg.f__pkl = os.path.join(self.intermediates_dir, os.path.basename(self.db.f__pkl))
        self.mutable_ref_pkg.slurp()

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
        self.db.change_file_paths(new_dir=os.path.dirname(self.db.f__pkl))
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

    def test_map_taxa_to_leaf_nodes(self):
        leaf_map = self.db.map_taxa_to_leaf_nodes(leaf_names=["g__Methanothrix", "g__Methanosarcina"])
        self.assertEqual(4, len(leaf_map["g__Methanothrix"]))
        self.assertEqual(13, len(leaf_map["g__Methanosarcina"]))
        return

    def test_get_internal_node_leaf_map(self):
        self.db.f__pkl = self.pkl_path
        self.db.slurp()
        node_map = self.db.get_internal_node_leaf_map()
        self.assertEqual(2*self.db.num_seqs, len(node_map))
        self.assertEqual(self.db.num_seqs, len(node_map[max(node_map.keys())]))
        return

    def test_match_taxon_to_internal_nodes(self):
        self.assertEqual([],
                         self.db.match_taxon_to_internal_nodes("c__Methanococcus"))
        self.assertEqual([132],
                         self.db.match_taxon_to_internal_nodes("g__Methanopyrus"))
        self.assertEqual([212, 378, 110, 273, 477, 26, 409, 131, 436, 287, 36, 448, 386,
                          294, 41, 415, 276, 495, 439, 438, 437, 417, 410, 296, 295, 27],
                         self.db.match_taxon_to_internal_nodes("d__Archaea"))
        self.assertEqual([497],
                         self.db.match_taxon_to_internal_nodes("r__Root"))
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
        self.db.f__pkl = self.new_pkl_path
        self.db.pickle_package()
        self.db.slurp()
        self.assertTrue("McrA" == self.db.prefix)
        return

    def test_slurp(self):
        self.db.f__pkl = self.pkl_path
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
        renamed_refpkg.f__pkl = os.path.join(self.disband_path, "McrA2_build.pkl")
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

    def test_write_edited_pkl(self):
        from treesapp.refpkg import write_edited_pkl

        # Test writing to a new directory
        retcode = write_edited_pkl(ref_pkg=self.mutable_ref_pkg, output_dir=self.intermediates_dir + "tmp/", overwrite=False)
        self.assertEqual(0, retcode)
        self.assertTrue(os.path.isfile(os.path.join(self.intermediates_dir, "tmp", os.path.basename(self.db.f__pkl))))

        # Test writing to the same file (should not work)
        retcode = write_edited_pkl(ref_pkg=self.mutable_ref_pkg, output_dir=self.intermediates_dir + "tmp/", overwrite=False)
        self.assertEqual(1, retcode)

        # Test overwriting the existing reference package
        self.mutable_ref_pkg.feature_annotations = {"Something": []}
        retcode = write_edited_pkl(ref_pkg=self.mutable_ref_pkg, output_dir="", overwrite=True)
        self.assertEqual(0, retcode)
        self.assertEqual(os.path.join(self.intermediates_dir, "tmp", os.path.basename(self.db.f__pkl)),
                         self.mutable_ref_pkg.f__pkl)
        # Reset the feature_annotations attribute to ensure the annotation is from the written, edited RefPkg
        self.mutable_ref_pkg.feature_annotations = {}
        self.mutable_ref_pkg.slurp()
        self.assertTrue("Something" in self.mutable_ref_pkg.feature_annotations)
        return

    def test_create_viewable_newick(self):
        newick_str = self.db.create_viewable_newick()
        self.assertEqual(15243, len(newick_str))
        return

    def test_convert_feature_indices_to_inodes(self):
        test_feature_map = {"11_McrA": "Hydrogenotrophic"}
        internal_node_feature_map = self.db.convert_feature_indices_to_inodes(test_feature_map)
        k, _ = internal_node_feature_map.popitem()
        self.assertIsInstance(k, int)
        self.assertTrue(296, k)

        test_feature_map = {'Gom-Arc1-GOS | MSRB0112': "Methanotrophic"}
        internal_node_feature_map = self.db.convert_feature_indices_to_inodes(test_feature_map)
        k, _ = internal_node_feature_map.popitem()
        self.assertTrue(39, k)
        return

    def test_add_feature_annotations(self):
        self.mutable_ref_pkg.add_feature_annotations(feature_name="test", feature_map={"g__Methanosarcina": "Aceticlastic",
                                                                                       "p__Candidatus Helarchaeota": "SCAO",
                                                                                       "c__Methanobacteria": "Hydrogenotrophic",
                                                                                       "116_McrA": "Hydrogenotrophic"})
        self.assertEqual(1, len(self.mutable_ref_pkg.feature_annotations))
        self.assertEqual(3, len(self.mutable_ref_pkg.feature_annotations["test"]))

        self.mutable_ref_pkg.add_feature_annotations(feature_name="test", feature_map={"p__Candidatus Bathyarchaeota": "SCAO",
                                                                                       "g__Candidatus Methanoperedens": "Methanotrophy"}, reset=True)
        self.assertEqual(2, len(self.mutable_ref_pkg.feature_annotations["test"]))

        self.mutable_ref_pkg.add_feature_annotations(feature_name="test", feature_map={}, reset=True)
        self.assertNotIn("test", self.mutable_ref_pkg.feature_annotations)
        return

    def test_edit(self):
        from treesapp.refpkg import edit
        new_value = "Z0002"
        edit(refpkg=self.mutable_ref_pkg,
             output_dir=os.path.dirname(self.mutable_ref_pkg.f__pkl),
             attributes=["refpkg_code", new_value],
             overwrite=True)

        # To ensure the observed attribute value is as we expect, assign the attribute a temporary value
        self.mutable_ref_pkg.refpkg_code = "tmp"

        self.mutable_ref_pkg.slurp()
        self.assertEqual(new_value, self.mutable_ref_pkg.refpkg_code)

        # Test updating a feature_annotation using a taxonomy-phenotype map with internal nodes
        edit(refpkg=self.mutable_ref_pkg,
             output_dir=os.path.dirname(self.mutable_ref_pkg.f__pkl),
             attributes=["feature_annotations", "Paralog"],
             phenotypes=utils.get_test_data("McrA_paralog_map.tsv"),
             overwrite=True,
             reset=False)
        self.assertEqual(1, len(self.mutable_ref_pkg.feature_annotations))
        self.assertEqual(2, len(self.mutable_ref_pkg.feature_annotations["Paralog"]))
        for clade_annot in self.mutable_ref_pkg.feature_annotations["Paralog"]:
            if clade_annot.name == "McrA":
                self.assertEqual(225, len(clade_annot.members))
            elif clade_annot.name == "AcrA":
                self.assertEqual(14, len(clade_annot.members))
            else:
                self.assertTrue(False)

        return


if __name__ == '__main__':
    unittest.main()
