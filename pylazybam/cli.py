#!/usr/bin/env python

"""
Lazy BAM parser for partially decoding reads
Useful for filtering reads into files or extracting
statistics.
Pure python and requires only standard library

Usage:
  pylazybam command [<cmd_arg>]...
  pylazybam -h | --help
  example --version

 Options:
   -h, --help       Show this message.
   -b, --beer       Drink beer.
   -r, --rock       Play AC/DC.
   -p, --pub=<p>    Which pub.
   --version        Print the version.
"""

from docopt import docopt
from pprint import pprint

if __name__ == '__main__':
    arguments = docopt(__doc__, version='FIXME')
    pprint(arguments)