#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 07:42:14 2020

"""
from graphunzip.arg_helpers import parse_args_command
import graphunzip.input_output as io
from graphunzip.run_submodules import run_submodule

import logging
import time

__version__ = "3.0.0"
__author__ = "Roland Faure"


def main():
    # Set up logger
    logging.basicConfig(
        level=0,
        format="%(asctime)s:%(levelname)s:%(module)s:%(message)s",
    )

    # Fetch sub-module command
    args_command = parse_args_command()
    command = args_command.command

    # Start time
    starttime = time.time()

    run_submodule(command)

    logging.info(f"Finished in {(time.time() - starttime):.2f} seconds")


if __name__ == "__main__":
    main()
