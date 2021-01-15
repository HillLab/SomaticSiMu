#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 28 21:27:05 2020

@author: davidchen
"""

# Include standard modules
import argparse

# Initiate the parser
parser = argparse.ArgumentParser()

# Add long and short argument
parser.add_argument("--generation", "-g", help="generation number")

# Read arguments from the command line
args = parser.parse_args()

# Check for --width
if args.generation:
    print("this is test generation %s" % args.generation)