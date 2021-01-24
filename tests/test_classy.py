import unittest
import pytest
import os

from collections import namedtuple
from shutil import rmtree

from .testing_utils import get_test_data


@pytest.fixture(scope="class")
def treesapp_instance(request):
    from treesapp.classy import TreeSAPP, ModuleFunction
    request.cls.db = TreeSAPP("")
    request.cls.db.stages = {0: ModuleFunction("orf-call", 0),
                             1: ModuleFunction("clean", 1),
                             2: ModuleFunction("search", 2),
                             3: ModuleFunction("align", 3),
                             4: ModuleFunction("place", 4),
                             5: ModuleFunction("classify", 5)}
    return


@pytest.mark.usefixtures("treesapp_instance")
class TreeSAPPClassTester(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp.fasta import FASTA
        self.output_dir = "./tests/TreeSAPP_classy_test_dir"
        self.fasta = get_test_data("marker_test_suite.faa")

        self.test_fasta = FASTA(self.fasta)
        self.test_fasta.load_fasta()

        self.db.output_dir = self.output_dir
        return

    def tearDown(self) -> None:
        if os.path.isdir(self.output_dir):
            rmtree(self.output_dir)
        return

    def test_increment_stage_dir(self):
        from treesapp.classy import Purity
        mock_purity_inst = Purity()
        # Set directory paths
        mock_purity_inst.output_dir = self.output_dir
        os.makedirs(os.path.join(mock_purity_inst.output_dir, "intermediates"), exist_ok=True)
        for stage in mock_purity_inst.stages.values():
            stage.dir_path = os.path.join(mock_purity_inst.output_dir, "intermediates", stage.name)

        # Set current stage
        mock_purity_inst.current_stage = mock_purity_inst.stages[0]
        mock_purity_inst.stages[1].run = False  # skip 'assign'
        self.assertEqual('clean', mock_purity_inst.current_stage.name)

        # Begin tests
        mock_purity_inst.increment_stage_dir()
        self.assertEqual('summarize', mock_purity_inst.current_stage.name)
        self.assertTrue(mock_purity_inst.current_stage.run)
        self.assertTrue(os.path.isdir(os.path.join(mock_purity_inst.output_dir, "intermediates", 'summarize')))
        # Increment past last stage
        mock_purity_inst.increment_stage_dir()
        self.assertEqual('summarize', mock_purity_inst.current_stage.name)
        self.assertTrue(mock_purity_inst.current_stage.run)
        return

    def test_furnish_with_arguments(self):
        args = namedtuple("args", ["molecule", "output", "input", "executables"])
        args.input = self.fasta
        args.output = self.output_dir
        args.molecule = "prot"
        args.executables = {'prodigal': '/home/connor/bin/prodigal', 'BMGE.jar': '/usr/local/bin/BMGE.jar',
                            'hmmbuild': '/usr/local/bin/hmmbuild',
                            'hmmalign': '/usr/local/bin/hmmalign',
                            'hmmsearch': '/usr/local/bin/hmmsearch',
                            'epa-ng': '/usr/local/bin/epa-ng', 'raxml-ng': '/usr/local/bin/raxml-ng'}
        self.db.furnish_with_arguments(args)
        self.assertEqual(len(args.executables), len(self.db.executables))
        self.assertEqual(self.fasta, self.db.input_sequences)
        self.assertEqual(args.molecule, self.db.molecule_type)
        return

    def test_check_previous_output(self):
        self.db.final_output_dir = self.output_dir + os.sep + "final_outputs"
        self.db.var_output_dir = self.output_dir + os.sep + "intermediates"
        self.db.check_previous_output(overwrite=True)
        self.assertTrue(os.path.isdir(self.output_dir))
        return

    def test_validate_new_dir(self):
        from treesapp.utilities import validate_new_dir
        output_dir = validate_new_dir(self.output_dir)
        self.assertEqual(os.sep, output_dir[-1])
        with pytest.raises(SystemExit):
            validate_new_dir(os.path.join(output_dir, "nope", "nope"))
        return

    def test_find_stage_dirs(self):
        ts_inst = self.db
        self.assertEqual(0, ts_inst.find_stage_dirs())
        ts_inst.change_stage_status("orf-call", False)
        ts_inst.change_stage_status("clean", False)
        self.assertEqual(2, ts_inst.find_stage_dirs())
        return

    def test_dedup_records(self):
        from treesapp.classy import dedup_records, get_header_info
        test_fasta_headers = get_header_info(self.test_fasta.header_registry)
        self.assertEqual(76, len(self.test_fasta.fasta_dict))
        deduped_records = dedup_records(self.test_fasta, test_fasta_headers)
        self.assertEqual(76, len(deduped_records))
        return

    def test_prep_logging(self):
        from treesapp.classy import prep_logging
        import logging
        # Reset logging for tests
        logging.getLogger().handlers.clear()
        # Should make the file properly
        prep_logging(log_file="./test_log.txt", verbosity=True)
        self.assertTrue(os.path.isfile("./test_log.txt"))
        os.remove("./test_log.txt")

        # Reset logging for tests
        logging.getLogger().handlers.clear()
        # Directory doesn't exist, should fail
        with pytest.raises(SystemExit):
            prep_logging(log_file=os.path.join(os.getcwd(), "must", "fail", "test_log.txt"),
                         verbosity=False)
        return


class EvaluatorTester(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp.classy import Evaluator
        from treesapp.refpkg import ReferencePackage
        self.output_dir = "./TreeSAPP_classy_test_dir"
        self.eval = Evaluator()
        # Set the basic attribtues
        self.eval.var_output_dir = self.output_dir
        # Set the instance's reference package
        self.eval.ref_pkg = ReferencePackage()
        self.eval.ref_pkg.f__json = get_test_data(os.path.join("refpkgs", "PuhA_build.pkl"))
        self.eval.ref_pkg.slurp()

        return

    def tearDown(self) -> None:
        if os.path.isdir(self.output_dir):
            rmtree(self.output_dir)
        return

    def test_new_taxa_test(self):
        rpp = self.eval.ref_pkg.prefix
        taxon = "p__Proteobacteria"
        taxa_test = self.eval.new_taxa_test(lineage="r__Root; d__Bacteria; p__Proteobacteria", tool="treesapp")
        self.assertEqual(taxon, taxa_test.taxon)
        self.assertEqual(os.path.join(self.output_dir, rpp, taxon) + os.sep,
                         taxa_test.intermediates_dir)
        self.assertEqual(os.path.join(self.output_dir, rpp, taxon,
                                      "treesapp_output", "final_outputs", "marker_contig_map.tsv"),
                         taxa_test.classification_table)
        return


if __name__ == '__main__':
    unittest.main()
