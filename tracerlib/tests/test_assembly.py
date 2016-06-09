from __future__ import print_function

import os
import unittest
from mock import patch
import sys
import pandas as pd
from pandas.util.testing import assert_frame_equal

from tracerlib import base_dir
from tracerlib.io import parse_IgBLAST
from tracerlib.tasks import Assembler


class TestAssemble(unittest.TestCase):

    expected_folder = os.path.join(base_dir, 'test_data', 'expected_summary')
    results_folder = os.path.join(base_dir, 'test_data', 'results', 'filtered_TCR_summary')

    def test_args(self):
        config_file = os.path.expanduser('~/.tracerrc')
        assert os.path.isfile(config_file), "Config file ~/.tracerrc not found"
        assemble_args = dict(
                ncores='1', config_file=config_file, resume_with_existing_files=False, species='Mmus',
                seq_method='imgt', fastq1='/not/a/file', fastq2='/not/a/file', cell_name='cell1',
                output_dir=self.results_folder, single_end=False, fragment_length=False, fragment_sd=False)

        # Launcher checks for fastq file existance
        with self.assertRaisesRegex(OSError, 'FASTQ file not found'):
            Assembler(**assemble_args).run()

        # Launcher checks for second fastq if paired end
        assemble_args['fastq2'] = None
        assemble_args['fastq1'] = os.path.join(base_dir, 'test_data', 'cell1_1.fastq')
        with self.assertRaisesRegex(AssertionError, 'Only one fastq'):
            Assembler(**assemble_args).run()

        # Check for fragment length with single end
        assemble_args['single_end'] = True
        with self.assertRaisesRegex(AssertionError, 'fragment length'):
            Assembler(**assemble_args).run()

    def test_parse_igblast(self):
        locus_names = ["TCRA", "TCRB"]
        imgt_seq_location = os.path.join(base_dir, "resources", "mouse", "imgt_sequences")
        const_seq_file = os.path.join(base_dir, "resources", "mouse", "constant_seqs.csv")
        out_dir = os.path.join(base_dir, "test_data", "results", "cell2")

        cell = parse_IgBLAST(locus_names, out_dir, 'cell2', imgt_seq_location, 'Mmus', 'imgt', const_seq_file)
        assert not cell._check_is_empty(), "No Cell results"


if __name__ == '__main__':
    unittest.main()
