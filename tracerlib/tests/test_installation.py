from __future__ import print_function

import os
import unittest
from mock import patch
import sys

from tracer.tracerlib.launcher import Launcher
from tracer import base_dir


class TestInstall(unittest.TestCase):

    def test_installation(self):
        test_args = ['tracer', 'test', '-p', '1', '-c', '/Users/ge2/.tracerrc']
        with patch.object(sys, 'argv', test_args):
            Launcher()
        with open(os.path.join(base_dir, 'test_data', 'expected_summary', 'recombinants.txt')) as fh:
            fh.readline()
            # Get first line of recombinants
            expected_recombinant = fh.readline()
        with open(os.path.join(base_dir, 'test_data', 'results', 'filtered_TCR_summary', 'recombinants.txt')) as fh:
            recombinant_results = fh.read()

        assert expected_recombinant in recombinant_results, 'Expected recombinant not found'


if __name__ == '__main__':
    unittest.main()
