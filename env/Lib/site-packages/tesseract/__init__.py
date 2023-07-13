#!/usr/bin/python
"""
tesseract
=========

A package for measuring the concentration of halos from Nbody simulations 
non-parametrically using Voronoi tessellation.

Subpackages
-----------
voro 
  Routines for running and manipulating data returned by the Voronoi
  tesselation routine vorovol.
nfw 
  Routines relating to fitting and determining properties of NFW profiles.
io 
  Routines for data input and output.
util
  Misc. utility routines
tests
  Routines for running and plotting different tests for the provided test halos.

"""

# Basic dependencies
import ConfigParser
import os
import shutil

# Initialize config file
_config_file_def = os.path.join(os.path.dirname(__file__),"default_config.ini")
_config_file_usr = os.path.expanduser("~/.tessrc")
if not os.path.isfile(_config_file_usr):
    print 'Creating user config file: {}'.format(_config_file_usr)
    shutil.copyfile(_config_file_def,_config_file_usr)

# Read config file
config_parser = ConfigParser.ConfigParser()
config_parser.optionxform = str
config_parser.read(_config_file_def)
config_parser.read(_config_file_usr) # Overrides defaults with user options

# General options
config = {}
if config_parser.has_option('general','outputdir'):
    config['outputdir'] = os.path.expanduser(config_parser.get('general','outputdir').strip())
else:
    config['outputdir'] = os.getcwd()

# Subpackages
import voro
import nfw
import io
import util
import tests


__all__ = ['voro','nfw','io','util','tests']


