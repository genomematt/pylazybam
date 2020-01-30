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
import struct
from typing import BinaryIO, Generator, Tuple, Dict
from pylazybam.bgzf import BgzfWriter as FileWriter
from pylazybam.decoders import *
from pylazybam.tags import *

# Parsing functions


def get_ref_index(alignment: bytes) -> int:
    """
    Extract the reference index from a BAM alignment

    Parameters
    ----------
    alignment : bytes
        A byte string of a bam alignment entry in raw binary format

    Returns
    -------
    int
        The zero based rank of the reference in the BAM header

    Notes
    -----
    The index can be converted to the reference name with
    pylazybam.bam.FileReader().index_to_ref[index]
    """
    return struct.unpack("<i", alignment[4:8])[0]


def get_pos(alignment: bytes) -> int:
    """
    Extract the one based position of this read

    Parameters
    ----------
    alignment : bytes
        A byte string of a bam alignment entry in raw binary format

    Returns
    -------
    int
        The one based left aligned position of this read on the reference

    """
    return struct.unpack("<i", alignment[8:12])[0]


def get_len_read_name(alignment: bytes) -> int:
    """
    Extract the length of the read name from a BAM alignment

    Parameters
    ----------
    alignment : bytes
        A byte string of a bam alignment entry in raw binary format

    Returns
    -------
    int
        The length of the read name
    """
    return struct.unpack("<B", alignment[12:13])[0]


def get_mapq(alignment: bytes) -> int:
    """
    Extract the read mapping quality score from a BAM alignment

    Parameters
    ----------
    alignment : bytes
        A byte string of a bam alignment entry in raw binary format

    Returns
    -------
    int
        The integer mapping quality score
    """
    return struct.unpack("<B", alignment[13:14])[0]


def get_bin(alignment: bytes) -> int:
    """
    Extract the BAI index bin from a BAM alignment

    Parameters
    ----------
    alignment : bytes
        A byte string of a bam alignment entry in raw binary format

    Returns
    -------
    int
        The integer value of the index bin
    """
    return struct.unpack("<H", alignment[14:16])[0]


def get_number_cigar_operations(alignment: bytes) -> int:
    """
    Extract the number of cigar operations from a BAM alignment

    Parameters
    ----------
    alignment : bytes
        A byte string of a bam alignment entry in raw binary format

    Returns
    -------
    int
        The number of operations in the cigar string

    """
    return struct.unpack("<H", alignment[16:18])[0]


def get_flag(alignment: bytes) -> int:
    """
    Extract the alignment flag from a BAM alignment

    Parameters
    ----------
    alignment : bytes
        A byte string of a bam alignment entry in raw binary format

    Returns
    -------
    int
        The alignment flag

    Notes
    -----
    Flag values can be tested on raw BAM alignment with pylazybam.bam.is_flag()
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
    return struct.unpack("<H", alignment[18:20])[0]


def get_len_sequence(alignment: bytes) -> int:
    """
    Extract the sequence length from a BAM alignment

    Parameters
    ----------
    alignment : bytes
        A byte string of a bam alignment entry in raw binary format

    Returns
    -------
    int
        the length of the sequence

    Notes
    -----
    The decoded quality string will be the same length as the sequence

    """
    return struct.unpack("<i", alignment[20:24])[0]


def get_pair_ref_index(alignment: bytes) -> int:
    """Extract the index identifying the reference sequence of this query
    sequences pair from a BAM alignment

    Parameters
    ----------
    alignment : bytes
        A byte string of a bam alignment entry in raw binary format

    Returns
    -------
    int
        The zero based rank of the reference in the BAM header

    Notes
    -----
    The index can be converted to the reference name with
    pylazybam.bam.FileReader().index_to_ref[index]

    """
    return struct.unpack("<i", alignment[24:28])[0]


def get_pair_pos(alignment: bytes) -> int:
    """Extract the one based position of this reads pair from a BAM alignment

    Parameters
    ----------
    alignment : bytes
        A byte string of a bam alignment entry in raw binary format

    Returns
    -------
    int
        The one based left aligned position of this reads pair on the reference

    """
    return struct.unpack("<i", alignment[28:32])[0]


def get_template_len(alignment: bytes) -> int:
    """Extract the template length from a BAM alignment bytestring

    Parameters
    ----------
    alignment : bytes
        A byte string of a bam alignment entry in raw binary format

    Returns
    -------
    int
        The integer length of the template
        (The distance between aligned read pairs)
    """
    return struct.unpack("<i", alignment[32:36])[0]


def get_read_name(alignment: bytes,
                  read_name_length: int) -> str:
    """Extract the read name in ASCII SAM format from a BAM alignment bytestring

    Parameters
    ----------
    alignment : bytes
        A byte string of a bam alignment entry in raw binary format

    len_read_name : int
        The length of the readname string
        eg from pylazybam.bam.get_len_read_name()
    Returns
    -------
    str
        The read name in ASCII SAM format

    """
    return alignment[36 : 35 + read_name_length].decode("utf-8")


def get_raw_read_name(alignment: bytes,
                      read_name_length: int) -> bytes:
    """Extract the raw readname from a BAM alignment bytestring

    Parameters
    ----------
    alignment : bytes
        A byte string of a bam alignment entry in raw binary format

    len_read_name : int
        The length of the readname string
        eg from pylazybam.bam.get_len_read_name()

    Returns
    -------
    bytes
        The raw base readname in BAM format as a binary bytestring

    """
    return alignment[36 : 36 + read_name_length]


def get_raw_cigar(alignment: bytes,
                  len_read_name: int,
                  number_cigar_operations: int) -> bytes:
    """Extract the raw cigar string from a BAM alignment bytestring

    Parameters
    ----------
    alignment : bytes
        A byte string of a bam alignment entry in raw binary format

    len_read_name : int
        The length of the readname string
        eg from pylazybam.bam.get_len_read_name()

    number_cigar_operations : int
        The number of cigar operations
        eg from pylazybam.bam.get_number_cigar_operations()

    Returns
    -------
    bytes
        The raw base cigar string in BAM format as a binary bytestring

    """
    start = 36 + len_read_name
    end = 36 + len_read_name + (4 * number_cigar_operations)
    return alignment[start:end]


def get_tag_bytestring(alignment: bytes,
                       len_read_name: int,
                       number_cigar_operations: int,
                       len_sequence: int, ) -> bytes:
    """Extract the raw tags from a BAM alignment bytestring

    Parameters
    ----------
    alignment : bytes
        A byte string of a bam alignment entry in raw binary format

    len_read_name : int
        The length of the readname string
        eg from pylazybam.bam.get_len_read_name()

    number_cigar_operations : int
        The number of cigar operations
        eg from pylazybam.bam.get_number_cigar_operations()

    len_sequence : int
        The length of the sequence and quality score strings
        eg from pylazybam.bam.get_len_sequence

    Returns
    -------
    bytes
        The raw tags in BAM format as a binary bytestring

    """
    start = (
        36
        + len_read_name
        + (4 * number_cigar_operations)
        + (len_sequence + 1) // 2
        + len_sequence
    )
    return alignment[start:]


def get_raw_sequence(alignment: bytes,
                len_read_name: int,
                number_cigar_operations: int,
                len_sequence: int, ) -> bytes:
    """Extract the raw sequence from a BAM alignment bytestring

    Parameters
    ----------
    alignment : bytes
        A byte string of a bam alignment entry in raw binary format

    len_read_name : int
        The length of the readname string
        eg from pylazybam.bam.get_len_read_name()

    number_cigar_operations : int
        The number of cigar operations
        eg from pylazybam.bam.get_number_cigar_operations()

    len_sequence : int
        The length of the sequence and quality score strings
        eg from pylazybam.bam.get_len_sequence

    Returns
    -------
    bytes
        The raw sequence in BAM format as a binary bytestring

    """
    start = 36 + len_read_name + (4 * number_cigar_operations)
    end = start + (len_sequence + 1) // 2
    return alignment[start:end]


def get_raw_base_qual(alignment: bytes,
                      len_read_name: int,
                      number_cigar_operations: int,
                      len_sequence: int, ) -> bytes:
    """Extract the raw base qualities from a BAM alignment bytestring

    Parameters
    ----------
    alignment : bytes
        A byte string of a bam alignment entry in raw binary format

    len_read_name : int
        The length of the readname string
        eg from pylazybam.bam.get_len_read_name()

    number_cigar_operations : int
        The number of cigar operations
        eg from pylazybam.bam.get_number_cigar_operations()

    len_sequence : int
        The length of the sequence and quality score strings
        eg from pylazybam.bam.get_len_sequence

    Returns
    -------
    bytes
        The raw base qualities in BAM format as a binary bytestring

    """
    start = (
        36
        + len_read_name
        + (4 * number_cigar_operations)
        + (len_sequence + 1) // 2
    )
    end = start + len_sequence
    return alignment[start:end]


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

        The detailed specification for the BAM format can be found at
        https://samtools.github.io/hts-specs/SAMv1.pdf

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
        >>> print(mybam.index_to_ref[get_ref_index(align)],
        >>>       get_pos(align),
        >>>       get_AS(align))

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
        self.index_to_ref: Dict[int,str] = dict(zip(range(len(self.refs)),
                                                    self.refs.keys()))
        self.ref_to_index: Dict[str,int] = dict(zip(self.refs.keys(),
                                                    range(len(self.refs))))
        self.index_to_ref[-1] = "*"
        self.alignments: Generator[bytes, None, None] = self._get_alignments()

    def __enter__(self):
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
                if not raw_blocksize:
                    break
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
