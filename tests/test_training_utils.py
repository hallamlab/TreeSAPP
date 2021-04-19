import unittest
import pytest
import os
import shutil

from tqdm import tqdm

from .testing_utils import get_test_data, get_treesapp_root


@pytest.fixture(scope="class")
def refpkg(request):
    from treesapp import refpkg
    request.cls.db = refpkg.ReferencePackage()
    request.cls.db.f__pkl = get_test_data(os.path.join("refpkgs", "PuhA_build.pkl"))
    request.cls.db.slurp()
    request.cls.db.validate()
    return


@pytest.mark.usefixtures("refpkg")
class ClassifierTester(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp.refpkg import ReferencePackage
        self.output_dir = "./training_test/"
        if not os.path.isdir(self.output_dir):
            os.mkdir(self.output_dir)

        self.puha_rp = ReferencePackage()
        self.puha_rp.f__pkl = get_test_data(os.path.join("refpkgs", "PuhA_build.pkl"))
        self.puha_rp.slurp()
        return

    def tearDown(self) -> None:
        if os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)
        return

    def test_augment_training_set(self):
        from treesapp.training_utils import augment_training_set
        import numpy as np
        test_array = np.array([0.33025841, 0.33025841, 0.07989843, 0.12523663, 0.13434813])
        self.assertEqual(3, len(augment_training_set(test_array)))
        self.assertEqual(1, len(augment_training_set(test_array, n_reps=0)))
        self.assertEqual((10, 5), augment_training_set(test_array, n_reps=10).shape)
        return

    def test_generate_pquery_data_for_trainer(self):
        from treesapp.training_utils import generate_pquery_data_for_trainer
        from treesapp.fasta import FASTA
        from treesapp.utilities import fetch_executable_path

        treesapp_dir = get_treesapp_root()
        executables = {}
        for dep in ["hmmbuild", "hmmalign", "hmmsearch", "epa-ng", "raxml-ng", "FastTree", "mafft", "BMGE.jar"]:
            executables[dep] = fetch_executable_path(dep, treesapp_dir)
        pbar = tqdm()
        test_taxon_one = "f__Bradyrhizobiaceae; g__Bradyrhizobium; s__Bradyrhizobium 'sp.' BTAi1"
        fasta_for_test = FASTA(file_name=get_test_data("ENOG4111FIN.txt"))
        fasta_for_test.load_fasta()
        pqueries = generate_pquery_data_for_trainer(ref_pkg=self.db, taxon=test_taxon_one, test_fasta=fasta_for_test,
                                                    training_seqs=['288000.BBta_6411'], rank="species", pbar=pbar,
                                                    output_dir=self.output_dir, executables=executables)
        pbar.close()
        self.assertEqual(1, len(pqueries))
        return

    def test_fetch_executable_path(self):
        from treesapp.utilities import fetch_executable_path
        from re import sub
        treesapp_dir = get_treesapp_root()
        exe_path = fetch_executable_path("BMGE.jar", treesapp_dir)
        self.assertEqual("/sub_binaries/BMGE.jar", sub(treesapp_dir, '', exe_path))
        return

    def test_load_training_data_frame(self):
        from treesapp.training_utils import load_training_data_frame
        from treesapp.phylo_seq import PQuery, PhyloPlace
        test_pqueries = {self.puha_rp.prefix: []}
        positives = {self.puha_rp.prefix: []}
        test_placements = [{'p': [[87, -12537.6198751071, 1.0, 0.0666942114, 0.0361806694]], 'n': ['-1']},
                           {'p': [[63, -12546.3946158455, 1.0, 0.1390379169, 0.2382098565]], 'n': ['-1']},
                           {'p': [[62, -12537.5328155479, 1.0, 0.0359374153, 0.0001]], 'n': ['-1']},
                           {'p': [[33, -12388.89952794, 1.0, 0.2410141, 0.3651239826]], 'n': ['-1']},
                           {'p': [[408, -1, 1.0, 0.0, 0.0]], 'n': ['-1']}]

        i = 0
        for placement in test_placements:
            pq = PQuery()
            pq.evalue = 3E-5
            pq.seq_len = 180
            pq.rank = "species"
            pq.seq_name = "seq" + str(i) + "_0"
            pq.consensus_placement = PhyloPlace(placement)
            test_pqueries[self.puha_rp.prefix].append(pq)
            if i % 2 == 0:  # Ensure some false positives are included by adding even numbered sequence names
                positives[self.puha_rp.prefix].append(pq.seq_name)
            i += 1

        # Fail due to bad internal node in last PQuery
        with pytest.raises(SystemExit):
            load_training_data_frame(pqueries=test_pqueries,
                                     refpkg_map={self.puha_rp.prefix: self.puha_rp},
                                     refpkg_positive_annots=positives)

        # Return empty DataFrame since the refpkg_positive_annots dictionary is empty
        self.assertEqual(0, len(load_training_data_frame(test_pqueries, {}, {})))
        test_pqueries[self.puha_rp.prefix].pop()

        # Test as intended
        training_df = load_training_data_frame(pqueries=test_pqueries,
                                               refpkg_map={self.puha_rp.prefix: self.puha_rp},
                                               refpkg_positive_annots=positives)
        self.assertEqual(4, len(training_df))
        self.assertTrue(set(training_df["tp"]))

        return


if __name__ == '__main__':
    unittest.main()
