#!/usr/bin/python
# -*- coding: utf-8 -*-
# Based on Rémy Greinhofer (rgreinho) tutorial on subcommands in docopt
# https://github.com/rgreinho/docopt-subcommands-example

from docopt import docopt

import input_output as io

from transform_gfa import gfa_to_fasta
from solve_ambiguities import solve_ambiguities

# from segment import check_if_all_links_are_sorted

from scipy import sparse
import numpy as np
import argparse
import os.path
import sys
import pickle  # reading and writing files
import time


class AbstractCommand:
    """Base class for the commands"""

    def __init__(self, command_args, global_args):
        """Initialize the commands."""
        self.args = docopt(self.__doc__, argv=command_args)
        self.global_args = global_args

    def execute(self):
        """Execute the commands"""
        raise NotImplementedError
