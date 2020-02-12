#!/usr/bin/env python3
# encoding: utf-8
"""
Module      : pylazybam/decoders.py
Description : A basic lazy bam parser.
Copyright   : (c) Matthew Wakefield, 2018-2020
License     : BSD-3-Clause
Maintainer  : matthew.wakefield@unimelb.edu.au
Portability : POSIX
"""

from array import array
from typing import Dict
from pylazybam import bam

FLAGS: Dict[str,int] = {
    "paired": 0x1,
    "aligned": 0x2,
    "unmapped": 0x4,
    "pair_unmapped": 0x8,
    "forward": 0x40,
    "reverse": 0x80,
    "secondary": 0x100,
    "qc_fail": 0x200,
    "duplicate": 0x400,
    "supplementary": 0x800,
}


def decode_sequence(raw_seq: bytes) -> str:
    """Decode raw BAM sequence into ASCII values

    Parameters
    ----------
    raw_seq : bytes
        The sequence section of a BAM alignment record as bytes
        eg the output of pybam.bam.get_raw_sequence()

    Returns
    -------
    str
        The ASCII encoded SAM representation of the query sequence

    """
    bam_bases = "=ACMGRSVTWYHKDBN"
    result = ""
    for x in array("B", raw_seq):
        result += bam_bases[x >> 4]
        result += bam_bases[x & 0b1111]
    return result


def decode_cigar(raw_cigar: bytes) -> str:
    """Decode raw BAM cigar strings into ASCII values

    Parameters
    ----------
    raw_cigar : bytes
        The cigar section of a BAM alignment record as bytes
        eg the output of pylazybam.bam.get_raw_cigar()

    Returns
    -------
    str
        The ASCII encoded SAM representation of the cigar string

    """
    codes = "MIDNSHP=X"
    cigar = [str(x >> 4) + codes[x & 0b1111] for x in array("I", raw_cigar)]
    return "".join(cigar)


def decode_base_qual(raw_base_qual: bytes,
                     offset: int = 33) -> str:
    """Decode raw BAM base quality scores into ASCII values

    Parameters
    ----------
    raw_base_qual : bytes
        The base quality section of a BAM alignment record as bytes
        eg the output from pylazybam.bam.get_raw_base_qual()

    offset : int
        The offset to add to the quality values when converting to ASCII

    Returns
    -------
    str
        The ASCII encoded SAM representation of the quality scores

    """
    return "".join([chr(q + offset) for q in list(raw_base_qual)])


def is_flag(alignment: bytes, flag: int) -> bool:
    """Test BAM flag values against a BAM alignment bytestring

    Parameters
    ----------
    alignment : bytes
        A byte string of a bam alignment entry in raw binary format

    flag :
        An integer representing the bitmask to compare to the

    Returns
    -------
    bool
        Returns true if all bits in the bitmask are set in the flag value

    Notes
    -----
    Common flag values are available from pylazybam.bam.FLAGS

    >>> print(pylazybam.bam.FLAGS)
    {"paired": 0x1,
    "aligned": 0x2,
    "unmapped": 0x4,
    "pair_unmapped": 0x8,
    "forward": 0x40,
    "reverse": 0x80,
    "secondary": 0x100,
    "qc_fail": 0x200,
    "duplicate": 0x400,
    "supplementary": 0x800,}

    See https://samtools.github.io/hts-specs/SAMv1.pdf for details.

    """
    return bool(bam.get_flag(alignment) & flag)
