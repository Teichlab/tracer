from __future__ import print_function

import os
import unittest
import glob

import shutil

from tracerlib.tasks import Builder
from tempfile import gettempdir

tempdir = gettempdir()


class TestBuild(unittest.TestCase):
    species = 'PurpleSnoutedHeffalump'

    def setUp(self):
        # Fake seq files
        seq_files = [('vseqs.fa', 'V01'), ('jseqs.fa', 'J01'),
                     ('cseqs.fa', 'C_01')]
        self.seq_files = [(os.path.join(tempdir, f), h) for f, h in seq_files]
        for filename, header in self.seq_files:
            with open(filename, 'w') as fh:
                fh.write(">%s\naaatttgggccc" % header)

        self.builder = Builder(
            ncores=1, species=self.species, receptor_name='TCR', locus_name='A',
            N_padding=10, V_seqs=self.seq_files[0][0],
            J_seqs=self.seq_files[1][0], C_seq=self.seq_files[2][0],
            force_overwrite=True)
        self.builder.init_dirs()
        self.vdjc_files = self.builder.copy_raw_files()

    def test_copy(self):
        v_fa = os.path.join(self.builder.species_dir, 'raw_seqs', 'TCR_A_V.fa')
        assert os.path.isfile(v_fa)
        with open(v_fa) as fh:
            header = fh.readline()
            assert '>V01' in header

    def test_recombinomes(self):
        self.builder.make_recombinomes(self.vdjc_files)
        # Read in recombinome data
        recombinome_file = os.path.join(
            self.builder.species_dir, 'combinatorial_recombinomes',
            'TCR_A.fa')
        with open(recombinome_file) as fh:
            header = fh.readline()
            assert 'V01_J01' in header
            tcr = fh.readline()
            assert 'aaatttgggcccNNN' in tcr

    def test_bowtie2_index(self):
        # Read in recombinome data
        recombinome_file = os.path.join(
            self.builder.species_dir, 'combinatorial_recombinomes',
            'TCR_A.fa')
        self.builder.make_bowtie2_index(recombinome_file)
        bt_files = glob.glob(os.path.join(
            self.builder.species_dir, 'combinatorial_recombinomes', '*.bt2'))
        assert len(bt_files)

    def test_igblastn_creation(self):
        missing_dbs = self.builder.make_igblast_db(self.vdjc_files)

    def tearDown(self):
        # Remove the species resources
        root_folder = self.builder.get_resources_root(self.species)
        shutil.rmtree(root_folder)
        for filename, _ in self.seq_files:
            if os.path.isfile(filename):
                os.unlink(filename)


if __name__ == '__main__':
    unittest.main()
