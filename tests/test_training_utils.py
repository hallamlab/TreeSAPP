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
    request.cls.db.f__json = get_test_data(os.path.join("refpkgs", "PuhA_build.pkl"))
    request.cls.db.slurp()
    request.cls.db.validate()
    return


@pytest.mark.usefixtures("refpkg")
class ClassifierTester(unittest.TestCase):
    def setUp(self) -> None:
        self.output_dir = "./training_test/"
        if not os.path.isdir(self.output_dir):
            os.mkdir(self.output_dir)
        return

    def tearDown(self) -> None:
        if os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)
        return

    def test_augment_training_set(self):
        from treesapp.training_utils import augment_training_set
        import numpy as np
        test_array = np.array([0.33025841, 0.33025841, 0.07989843, 0.12523663, 0.13434813])
        result = augment_training_set(test_array)
        self.assertTrue(3, len(result))

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


if __name__ == '__main__':
    unittest.main()
