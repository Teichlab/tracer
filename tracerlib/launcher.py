from __future__ import print_function
import matplotlib as mpl

mpl.use('pdf')
import argparse
import sys

from tracerlib.tasks import Assembler, Summariser, Tester


def launch():
    parser = argparse.ArgumentParser(
        description='TraCeR: reconstruction of TCR sequences from single-cell RNAseq data',
        usage=''' tracer <mode> [<args>]

              Modes are :

              - assemble: assemble TCR sequences from single-cell RNA-sequencing reads
              - summarise: summarise TCR sequences from set of cells, build clonotype networks
              - test : use a small dataset from three cells to test TraCeR installation

              use tracer <mode> -h for specific help
              ''')
    parser.add_argument('mode', metavar="<MODE>", help='tracer mode (assemble, summarise or get_test_data)',
                        choices=['assemble', 'summarise', 'summarize', 'test'])
    args = parser.parse_args(sys.argv[1:2])

    task_mapper = {
        'assemble': Assembler,
        'summarise': Summariser,
        'summarize': Summariser,
        'test': Tester
    }

    if args.mode not in task_mapper:
        print('Unrecognised mode')
        parser.print_help()
        exit(1)

    Task = task_mapper[args.mode]
    Task().run()
