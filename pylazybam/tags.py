#!/usr/bin/env python3
# encoding: utf-8
"""
Module      : bam.py
Description : A basic lazy bam parser.
Copyright   : (c) Matthew Wakefield, 2018-2020
License     : BSD-3-Clause
Maintainer  : matthew.wakefield@unimelb.edu.au
Portability : POSIX
"""

import struct, re
from typing import Any

# define a very negative int to avoid using MIN32INT and type conversion
MIN32INT: int = -2147483648

def get_AS(tag_bytes: bytes,
           no_tag: Any = MIN32INT) -> int:
    """Extract the high scoring alignment score from an AS tag in a raw BAM
    alignment bytestring

    Parameters
    ----------
        tag_bytes : bytes
            a bytestring containing bam formatted tag elements

        no_tag : Any
            return value for when tag not found (default: MIN32INT)

    Returns
    -------
        AS Tag Value : int
            the integer value of the AS tag
            returns the value of no_tag if tag absent (default:MIN32INT)

    Raises
    ------
        ValueError
            raises a ValueError if more than one tag match

    Notes
    -----
    Recommended try accept for use on raw alignment with fall back
    to calling on only the tag byte string.

    Please test carefully on your BAM output as in complicated output the
    regular expression based extraction of the tag can be error prone

    """
    match = re.findall(b"ASC.", tag_bytes)
    if not match:
        return no_tag
    elif len(match) != 1:
        raise ValueError(
            (
                f"More than one match to b'ASC.' was found in {tag_bytes}"
                "meaning that more than one value of the tag is present"
                "or that another tag contains a match as part of its value"
            )
        )
    else:
        return struct.unpack("<xxxB", match[0])[-1]  # type: int


def get_XS(tag_bytes: bytes,
           no_tag: Any = MIN32INT) -> int:
    """Extract the suboptimal alignment score from an XS tag in a raw BAM
    alignment bytestring

    Parameters
    ----------
        tag_bytes : bytes
            a bytestring containing bam formatted tag elements

        no_tag : Any
            return value for when tag not found (default: MIN32INT)

    Returns
    -------
        XS Tag Value : int
            the integer value of the XS tag
            returns the value of no_tag if tag absent (default:MIN32INT)

    Raises
    ------
        ValueError
            raises a ValueError if more than one tag match

    Notes
    -----
    This function is for the genome aligner definition of XS where XS:i:<int>
    is the alignment score of the suboptimal alignment.
    This is not the same as the spliced aligner XS tag XS:C:<str> that
    represents the strand on which the intron occurs (equiv to TS:C:<str>)

    Recommended try accept for use on raw alignment with fall back
    to calling on only the tag byte string.

    Please test carefully on your BAM output as in complicated output the
    regular expression based extraction of the tag can be error prone

    """
    match = re.findall(b"XSC.", tag_bytes)
    if not match:
        return no_tag
    elif len(match) != 1:
        raise ValueError(
            (
                f"More than one match to b'XSC.' was found in {tag_bytes}"
                "meaning that more than one value of the tag is present"
                "or that another tag contains a match as part of its value"
            )
        )
    else:
        return struct.unpack("<xxxB", match[0])[-1]


def get_ZS(tag_bytes: bytes,
           no_tag: Any = MIN32INT) -> int:
    """Extract the suboptimal alignment score from the ZS tag in a raw BAM
    alignment bytestring

    Parameters
    ----------
        tag_bytes : bytes
            a bytestring containing bam formatted tag elements

        no_tag : Any
            return value for when tag not found (default: MIN32INT)

    Returns
    -------
        ZS Tag Value : int
            the integer value of the ZS tag
            returns the value of no_tag if tag absent (default:MIN32INT)

    Raises
    ------
        ValueError
            raises a ValueError if more than one tag match

    Notes
    -----
    ZS is the equivalent to XS:i:<int> tag in some spliced aligners
    including HISAT2.

    Recommended try accept for use on raw alignment with fall back
    to calling on only the tag byte string.

    Please test carefully on your BAM output as in complicated output the
    regular expression based extraction of the tag can be error prone

    """
    match = re.findall(b"ZSC.", tag_bytes)
    if not match:
        return no_tag
    elif len(match) != 1:
        raise ValueError(
            (
                f"More than one match to b'ZSC.' was found in {tag_bytes}"
                "meaning that more than one value of the tag is present"
                "or that another tag contains a match as part of its value"
            )
        )
    else:
        return struct.unpack("<xxxB", match[0])[-1]


def get_MD(tag_bytes: bytes,
           no_tag: Any = None) -> str:
    """Extract the MD tag from a raw BAM alignment bytestring

    Parameters
    ----------
        tag_bytes : bytes
            a bytestring containing bam formatted tag elements

        no_tag : Any
            return value for when tag not found (default: None)

    Returns
    -------
        MD Tag Value : str
            an ASCII string representing the SAM format value of the MD tag
            returns the value of no_tag if tag absent (default: None)

    Raises
    ------
        ValueError
            raises a ValueError if more than one tag match

    Notes
    -----
    Recommended try accept for use on raw alignment with fall back
    to calling on only the tag byte string.

    Please test carefully on your BAM output as in complicated output the
    regular expression based extraction of the tag can be error prone

    """
    match = re.findall(b"MDZ[0-9ACGTN^]+\x00", tag_bytes)
    if not match:
        return no_tag
    elif len(match) != 1:
        raise ValueError(
            (f"More than one match to b'MDZ[0-9ACGTN^]+\x00' was found in "
             f"{tag_bytes} "
             "meaning that more than one value of the tag is present"
             "or that another tag contains a match as part of its value")
        )
    else:
        return match[0][3:-1].decode()

