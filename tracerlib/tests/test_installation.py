from __future__ import print_function

import os
import unittest
from mock import patch
import sys
import pandas as pd
from pandas.util.testing import assert_frame_equal

from tracerlib.launcher import Launcher
from tracerlib import base_dir


class TestInstall(unittest.TestCase):

    expected_folder = os.path.join(base_dir, 'test_data', 'expected_summary')
    results_folder = os.path.join(base_dir, 'test_data', 'results', 'filtered_TCR_summary')

    def test_installation(self):
        test_args = ['tracer', 'test', '-p', '1', '-c', os.path.expanduser('~/.tracerrc')]
        with patch.object(sys, 'argv', test_args):
            Launcher()

    def test_recombinants(self):

        # Assert the generated recombinants are identical
        def read_recombinants(rec_file):

            recombinants = pd.read_csv(rec_file, sep='\t')
            recombinants.dropna(how='all', inplace=True)
            recombinants.sort_values(by='recombinant_id', inplace=True)
            recombinants.reset_index(inplace=True, drop=True)

            # Reconstructed lengths differ as trinity is stochastic
            recombinants.drop(['reconstructed_length'], axis=1, inplace=True)

            return recombinants

        expected_recombinants = read_recombinants(os.path.join(self.expected_folder, 'recombinants.txt'))
        result_recombinants = read_recombinants(os.path.join(self.results_folder, 'recombinants.txt'))

        assert_frame_equal(expected_recombinants, result_recombinants)

    def test_clonotype_sizes(self):

        # Look at clonotype size files
        expected_clonosize = pd.read_csv(os.path.join(self.expected_folder, 'clonotype_sizes.txt'), sep='\t')
        results_clonosize = pd.read_csv(os.path.join(self.results_folder, 'clonotype_sizes.txt'), sep='\t')
        assert_frame_equal(expected_clonosize, results_clonosize)


if __name__ == '__main__':
    unittest.main()
