import unittest
import pytest
import os
import shutil

from .testing_utils import get_test_data, get_treesapp_root


class PlacementTrainerTestCase(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp.fasta import FASTA
        from treesapp.refpkg import ReferencePackage
        from treesapp.utilities import fetch_executable_path
        self.output_dir = "pt_test_dir"
        if not os.path.isdir(self.output_dir):
            os.mkdir(self.output_dir)
        # FASTA instance for tests
        self.test_fasta = FASTA(file_name=get_test_data("ENOG4111FIN.txt"))
        self.test_fasta.load_fasta()
        # A different fasta file containing non-homologous sequences
        self.bad_fasta = FASTA(file_name=get_test_data("EggNOG_McrA.faa"))
        self.bad_fasta.load_fasta()
        # ReferencePackage instance for tests
        self.test_refpkg = ReferencePackage("PuhA")
        self.test_refpkg.f__json = get_test_data(os.path.join("refpkgs", "PuhA_build.pkl"))
        self.test_refpkg.slurp()
        # Executables dictionary
        self.exes = {}
        self.treesapp_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.sep
        for dep in ["hmmbuild", "hmmalign", "raxml-ng", "mafft", "BMGE.jar"]:
            self.exes[dep] = fetch_executable_path(dep, get_treesapp_root())
        return

    def tearDown(self) -> None:
        if os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)
        return

    def test_reduce_examples(self):
        from treesapp.placement_trainer import reduce_examples
        test_size = 10  # Maximum is 23 (number of seq names in the dict)
        candidates_dict = {'p': {'bacp1': ['seq1', 'seq2', 'seq3', 'seq4'],
                                 'bacp2': ['seq5', 'seq6', 'seq7'],
                                 'arcp1': ['seq8', 'seq9', 'seq10']},
                           'o': {'bacc1': ['seq1', 'seq2'],
                                 'bacc2': ['seq3', 'seq4'],
                                 'bacc3': ['seq5', 'seq6', 'seq7'],
                                 'arcc1': ['seq8', 'seq9'],
                                 'arcc2': ['seq10']},
                           's': {'bacs1': ['seq1'],
                                 'bacs2': ['seq3'],
                                 'arcs1': ['seq9']}}
        trimmed_candidates = reduce_examples(candidate_seqs=candidates_dict, max_examples=test_size)
        i = 0
        for r, t in trimmed_candidates.items():
            i += sum([len(n) for n in t.values()])
        self.assertEqual(test_size, i)

    def test_clade_exclusion_phylo_placement(self):
        from treesapp.placement_trainer import clade_exclusion_phylo_placement
        # Ensure a rank isn't included in the returned dictionary if there were no PQueries associated with it
        train_seqs = {}
        pqueries = clade_exclusion_phylo_placement(rank_training_seqs=train_seqs,
                                                   test_fasta=self.test_fasta, ref_pkg=self.test_refpkg,
                                                   executables={}, min_seqs=1)
        self.assertEqual({}, pqueries)

        # Return a dictionary with an empty rank for which classification failed due to FASTA discrepancy
        train_seqs = {"species": {'d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Bradyrhizobiaceae; g__Bradyrhizobium; s__Bradyrhizobium sp. ORS 278': ['114'],
                                  'd__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhodospirillales; f__Rhodospirillaceae; g__Rhodospirillum; s__Rhodospirillum centenum': ['4146']}}
        with pytest.raises(KeyError):
            clade_exclusion_phylo_placement(rank_training_seqs=train_seqs,
                                            test_fasta=self.test_fasta, ref_pkg=self.test_refpkg,
                                            executables=self.exes, min_seqs=1, output_dir=self.output_dir)

        # Return a dictionary with an empty rank for which classification failed because its non-homologous
        train_seqs = {"species": {'d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Bradyrhizobiaceae; g__Bradyrhizobium; s__Brady': ['420247.Msm_1015'],
                                  'd__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhodospirillales; f__Rhodospirillaceae; g__Rhodospirillum; s__Rhodo': ['523846.Mfer_0784']},
                      "class": {'d__Bacteria; p__Proteobacteria; c__Betaproteobacteria': ['420247.Msm_1015',
                                                                                          '634498.mru_1924'],
                                'd__Bacteria; p__Proteobacteria; c__Gammaproteobacteria': ['523846.Mfer_0784',
                                                                                           '79929.MTBMA_c15480']}}
        pqueries = clade_exclusion_phylo_placement(rank_training_seqs=train_seqs,
                                                   test_fasta=self.bad_fasta, ref_pkg=self.test_refpkg,
                                                   executables=self.exes, min_seqs=3, output_dir=self.output_dir)
        self.assertEqual({"class": {'d__Bacteria; p__Proteobacteria; c__Betaproteobacteria': [],
                                    'd__Bacteria; p__Proteobacteria; c__Gammaproteobacteria': []},
                          "species": {'d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Bradyrhizobiaceae; g__Bradyrhizobium; s__Brady': [],
                                      'd__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhodospirillales; f__Rhodospirillaceae; g__Rhodospirillum; s__Rhodo': []}},
                         pqueries)
        return

    def test_evo_dists_from_pqueries(self):
        from treesapp.placement_trainer import evo_dists_from_pqueries
        from treesapp.phylo_seq import PQuery, PhyloPlace
        mock_pquery = PQuery()
        mock_pquery.consensus_placement = PhyloPlace()
        mock_pqueries = {"class":  {"tax_class": [mock_pquery]}}
        # Missing a required rank, 'species'
        evo_dists = evo_dists_from_pqueries(mock_pqueries)
        self.assertEqual({}, evo_dists)

        # Contains all training ranks, but it's empty
        mock_pqueries.update({"species": {"tax_species": []}})
        evo_dists = evo_dists_from_pqueries(mock_pqueries)
        self.assertEqual({}, evo_dists)

        # The input pqueries dictionary contain both ranks but incorrect type
        mock_pqueries.update({"species": {"tax_species": mock_pquery}})
        with pytest.raises(TypeError):
            evo_dists_from_pqueries(mock_pqueries)

        # The input pqueries dictionary is as expected
        mock_pqueries = {"class": {"tax_class": [mock_pquery, mock_pquery]},
                         "species": {"tax_species": [mock_pquery]}}
        evo_dists = evo_dists_from_pqueries(mock_pqueries)
        self.assertEqual(2, len(evo_dists["class"]))
        self.assertEqual(1, len(evo_dists["species"]))

        return


if __name__ == '__main__':
    unittest.main()
