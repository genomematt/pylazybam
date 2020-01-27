#!/usr/bin/env python3
# encoding: utf-8
"""
Module      : bam.py
Description : A basic lazy bam parser.
Copyright   : (c) Matthew Wakefield, 2018-2019 
License     : BSD-3-Clause 
Maintainer  : matthew.wakefield@unimelb.edu.au 
Portability : POSIX
"""
from __future__ import annotations
import sys, os
import struct
from array import array
import re
from typing import BinaryIO, Generator, Tuple, Dict
from pylazybam.bgzf import BgzfWriter as FileWriter

FLAGS = {
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


# Parsing functions


def get_ref_index(alignment: bytes) -> int:
    return struct.unpack("<i", alignment[4:8])[0]


def get_pos(alignment: bytes) -> int:
    return struct.unpack("<i", alignment[8:12])[0]


def get_len_read_name(alignment: bytes) -> int:
    return struct.unpack("<B", alignment[12:13])[0]


def get_mapq(alignment: bytes) -> int:
    return struct.unpack("<B", alignment[13:14])[0]


def get_bin(alignment: bytes) -> int:
    return struct.unpack("<H", alignment[14:16])[0]


def get_number_cigar_opperations(alignment: bytes) -> int:
    return struct.unpack("<H", alignment[16:18])[0]


def get_flag(alignment: bytes) -> int:
    return struct.unpack("<H", alignment[18:20])[0]


def get_len_sequence(alignment: bytes) -> int:
    return struct.unpack("<i", alignment[20:24])[0]


def get_pair_ref_index(alignment: bytes) -> int:
    return struct.unpack("<i", alignment[24:28])[0]


def get_pair_pos(alignment: bytes) -> int:
    return struct.unpack("<i", alignment[28:32])[0]


def get_template_len(alignment: bytes) -> int:
    return struct.unpack("<i", alignment[32:36])[0]


def get_read_name(alignment: bytes,
                  read_name_length: int) -> str:
    return alignment[36 : 35 + read_name_length].decode("utf-8")


def get_raw_read_name(alignment: bytes,
                      read_name_length: int) -> bytes:
    return alignment[36 : 36 + read_name_length]


def get_raw_cigar(alignment: bytes,
                  len_read_name: int,
                  number_cigar_opperations: int) -> bytes:
    start = 36 + len_read_name
    end = 36 + len_read_name + (4 * number_cigar_opperations)
    return alignment[start:end]


def get_tag_bytestring(alignment: bytes,
                       len_read_name: int,
                       number_cigar_opperations: int,
                       len_sequence: int, ) -> bytes:
    start = (
        36
        + len_read_name
        + (4 * number_cigar_opperations)
        + (len_sequence + 1) // 2
        + len_sequence
    )
    return alignment[start:]


def get_raw_seq(alignment: bytes,
                len_read_name: int,
                number_cigar_opperations: int,
                len_sequence: int, ) -> bytes:
    start = 36 + len_read_name + (4 * number_cigar_opperations)
    end = start + (len_sequence + 1) // 2
    return alignment[start:end]


def get_raw_base_qual(alignment: bytes,
                      len_read_name: int,
                      number_cigar_opperations: int,
                      len_sequence: int, ) -> bytes:
    start = (
        36
        + len_read_name
        + (4 * number_cigar_opperations)
        + (len_sequence + 1) // 2
    )
    end = start + len_sequence
    return alignment[start:end]


def get_AS(tag_bytes: bytes) -> int:
    """A function to extract an AS tag from a raw bam byte string
    Can be used on a raw alignment entry.
    Raises ValueError if multiple matches.
    Recommended try accept for use on raw alignment with fall back
    to calling on only the tag byte string.
    """
    match = re.findall(b"ASC.", tag_bytes)
    if not match:
        return None
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


def get_XS(tag_bytes: bytes) -> int:
    """A function to extract an XS tag from a raw bam byte string
    for genome aligner definition of XS where XS:i:<int> is the
    alignment score of the next best alignment.
    This is not the same as the spliced aligner XS tag that represents
    the strand on which the intron occurs
    Can be used on a raw alignment entry.
    Raises ValueError if multiple matches.
    Recommended try accept for use on raw alignment with fall back
    to calling on only the tag byte string.
    """
    match = re.findall(b"XSC.", tag_bytes)
    if not match:
        return None
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


def get_ZS(tag_bytes: bytes) -> int:
    """A function to extract an ZS tag from a raw bam byte string
    ZS is the equivalent to XS:i:<int> tag in some spliced aligners
    including HISAT2.
    Can be used on a raw alignment entry.
    Raises ValueError if multiple matches.
    Recommended try accept for use on raw alignment with fall back
    to calling on only the tag byte string.
    """
    match = re.findall(b"ZSC.", tag_bytes)
    if not match:
        return None
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


def get_MD(tag_bytes: bytes) -> str:
    """A function to extract an MD tag from a raw bam byte string
    Can be used on a raw alignment entry.
    Recommended try accept for use on raw alignment with fall back
    to calling on only the tag byte string.

    Parameters
    ----------
        tag_bytes : bytes
            a bytestring containing bam formatted tag elements

    Returns
    -------
        MD Tag Value : str
            an ascii string representing the value of the MD tag

    Raises
    ------
        ValueError
            raises a ValueError if more than one tag match

    """
    match = re.findall(b"MDZ[0-9ACGTN^]+\x00", tag_bytes)
    if not match:
        return None
    elif len(match) != 1:
        raise ValueError(
            (f"More than one match to b'MDZ[0-9ACGTN^]+\x00' was found in "
             f"{tag_bytes} "
             "meaning that more than one value of the tag is present"
             "or that another tag contains a match as part of its value")
        )
    else:
        return match[0][3:-1].decode()


def decode_seq(raw_seq: bytes) -> str:
    bam_bases = "=ACMGRSVTWYHKDBN"
    result = ""
    for x in array("B", raw_seq):
        result += bam_bases[x >> 4]
        result += bam_bases[x & 0b1111]
    return result


def decode_cigar(raw_cigar: bytes) -> str:
    codes = "MIDNSHP=X"
    cigar = [str(x >> 4) + codes[x & 0b1111] for x in array("I", raw_cigar)]
    return "".join(cigar)


def decode_base_qual(raw_base_qual: bytes,
                     offset: int = 33) -> str:
    return "".join([chr(q + offset) for q in list(raw_base_qual)])


def is_flag(alignment: bytes, flag: int) -> bool:
    return bool(bam.get_flag(a) & flag)


class FileReader:
    """A Pure Python Lazy Bam Parser Class

    Parameters
    ----------
        ubam : BinaryIO
            An binary (bytes) file or stream containing a valid uncompressed
            bam file conforming to the specification.

    Yields
    ------
        align : bytes
            A byte string of a bam alignment entry in raw binary format

    Attributes
    ----------
        raw_header : bytes
            The raw bytestring representing the bam header
        raw_refs : bytes
            The raw bytestring representing the reference sequences
        header : str
            The ASCII representation of the header
        refs : Dict[str:int]
            A dictionary of reference_name keys with reference_length
        ref_to_index : Dict[str:int]
            A dictionary mapping reference names to the bam numeric identifier
        index_to_ref : Dict[int:str]
            A dictionary mapping bam reference numeric identifiers to names

    Notes
    -----
        It is advisable not to call the private functions or operate directly
        on the underlying file object.

    Example
    -------
        This class requires an uncompressed bam file as input.
        For this example we will use gzip to decompress the test file included
        as a resource in the package

        >>> import gzip
        >>> from pkg_resources import resource_stream
        >>> ubam = gzip.open('/tests/data/paired_end_testdata_human.bam'))

        A bam filereader object can then be created and headers inspected
        >>> mybam = bam.FileReader(ubam)
        >>> print(mybam.header)

        The filereader object is an iterator and yields alignments in raw format
        >>> align = next(mybam)

        Alignments can be processed using functions from pylazybam.bam
        >>> print(mybam.index_to_ref[get_ref_index(align)], get_AS(align))

    """

    def __init__(self, ubam: BinaryIO):
        """
        """
        self._ubam = ubam
        self.magic = self._ubam.read(4)
        if self.magic != b"BAM\x01":
            raise ValueError(
                (f"Incorrect start to uncompressed bam: {self.magic} "
                "not b'BAM\x01'. Check the file has been decompressed "
                "before being passed to this class")
            )
        self.raw_header, self.header = self._read_header()
        self.raw_refs, self.refs = self._read_refs()
        self.index_to_ref = dict(zip(range(len(self.refs)), self.refs.keys()))
        self.ref_to_index = dict(zip(self.refs.keys(), range(len(self.refs))))
        self.index_to_ref[-1] = "*"
        self.alignments = self._get_alignments()

    def __enter__(self: FileReader) -> FileReader:
        """Return self for use in WITH statement."""
        return self

    def __exit__(self, type, value, traceback):
        """Tidy up at end of WITH statement."""
        self._ubam.close()

    def close(self):
        """Close input file"""
        return self._ubam.close()

    def _read_header(self) -> Tuple[bytes, str]:
        """Utility function for reading bam file headers
        Called by init to create self.raw_header and self.header"""
        raw_header_length = self._ubam.read(4)
        header_length = struct.unpack("<i", raw_header_length)[0]
        raw_header = raw_header_length + self._ubam.read(header_length)
        header = raw_header[4:].decode("utf-8")
        return (raw_header, header)

    def _read_refs(self) -> Tuple[bytes, Dict[str, int]]:
        """Utility function for reading bam file headers
        Called by init to create self.raw_refs and self.refs"""
        ref_buffer = b""
        refs = {}

        raw_nref = self._ubam.read(4)
        ref_buffer += raw_nref
        n_ref = struct.unpack("<i", raw_nref)[0]
        for i in range(n_ref):
            raw_lname = self._ubam.read(4)
            ref_buffer += raw_lname
            l_name = struct.unpack("<i", raw_lname)[0]
            raw_ref_name = self._ubam.read(l_name)
            ref_buffer += raw_ref_name
            ref_name = raw_ref_name[:-1].decode("utf-8")
            raw_ref_length = self._ubam.read(4)
            ref_buffer += raw_ref_length
            ref_length = struct.unpack("<i", raw_ref_length)[0]
            refs[ref_name] = ref_length
        return (ref_buffer, refs)

    def _get_alignments(self) -> Generator[bytes, None, None]:
        while self._ubam:
            try:
                raw_blocksize = self._ubam.read(4)
                block_size = struct.unpack("<i", raw_blocksize)[0]
                alignment = raw_blocksize + self._ubam.read(
                    block_size
                )  # add back size so alignment record intact
                yield alignment
            except ValueError:
                if not self._ubam.read():
                    break
                else:
                    raise

    def __iter__(self):
        return self

    def __next__(self):
        return next(self.alignments)
