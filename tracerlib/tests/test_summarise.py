from __future__ import print_function

import os
import shutil
import unittest
import pandas as pd
import sys
import pandas as pd
from pandas.util.testing import assert_frame_equal

from tracerlib import base_dir
from tracerlib.tasks import Summariser
from tracerlib.io import parse_IgBLAST, parse_invariant_cells
from tracerlib.tasks import Assembler


class TestSummarise(unittest.TestCase):

    def setUp(self):
        self.results_dir = os.path.join(base_dir, 'test_data', 'results')
        self.summariser = Summariser(
            use_unfiltered=False, keep_invariant=False, graph_format='pdf',
            no_networks=False, root_dir=self.results_dir,
            receptor_name='TCR', loci=['A', 'B'], species='Mmus')

    def test_clonotype_matrix(self):
        self.summariser.run()

        cell_data = pd.read_csv(
            os.path.join(self.results_dir, "filtered_TCRAB_summary",
                         "cell_data.csv"),
            index_col="cell_name")
        assert cell_data.loc["cell1", "A_productive"] is not None
        assert cell_data.loc["cell1", "group_size"] == 2


if __name__ == "__main__":
    unittest.main()
