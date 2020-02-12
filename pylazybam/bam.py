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
from pathlib import Path
from typing import BinaryIO, Generator, Tuple, Dict
from pylazybam.bgzf import BgzfWriter
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

class _FileBase:
    def __init__(self):
        pass

    def get_full_raw_header(self):
        """
        Get a complete BAM header with all elements suitable for writing to a
        bam.FileWriter

        Returns
        -------
        bytes
            The complete raw BAM header including the BAM magic and reference
            block
        """
        self.update_header_length()
        return self.magic + self.raw_header + self.raw_refs

    def update_header_length(self, raw_header=None):
        """
        Update the length of the SAM format text component of the header

        Parameters
        ----------
        raw_header : bytes
            optional raw SAM text header as bytes to process
            Default self.raw_header

        Returns
        -------
        bytes
            returns the length corrected raw header if a raw header is provided

        """
        inplace = False
        if not raw_header:
            inplace = True
            raw_header = self.raw_header
        encoded_length = struct.unpack("<i", raw_header[:4])[0]
        actual_length = len(raw_header[4:])
        if not encoded_length == actual_length:
            raw_header = (struct.pack("<i", actual_length)
                               + raw_header[4:])
        if inplace:
            self.raw_header = raw_header
        else:
            return raw_header

    def get_updated_header(self,
                  id: str,
                  program: str,
                  version: str,
                  command: str = None,
                  description: str = None,
                  raw_header = None,
                  ):
        """
        Get a modified version of the header with additional program information
        in a @PG BAM header suitable for writing to a new output file.

        Does not include the BAM magic or the BAM reference block.

        Note that the header should be written after self.magic and before
        self.raw_refs

        Parameters
        ----------
        id : str
            a unique identifier of this action of the program on the BAM

        program : str
            the name of the program used to process the BAM

        version : str
            version number of the program used to process the BAM

        command : str
            the command and arguments to the program used to process the BAM

        description : str
            optional description

        raw_header : bytes
            a raw BAM formatted header without BAM magic or REF block
            used for recursive modification to add multiple

        Notes
        -----
        Parameters will be utf-8 encoded
        If there are existing @PG records the last record will be used as PP

        See Also
        --------
        update_header - calls this method to modify self.raw_header in place
        get_full_raw_header
        """
        if raw_header == None:
            raw_header = self.raw_header
        PG_match = re.search(b'@PG.+',raw_header)
        CO_match = re.search(b'@CO.+',raw_header)
        CO_start = CO_match.span()[0] if CO_match else len(raw_header)
        PG_start = PG_match.span()[0] if PG_match else CO_start
        raw_HD_and_SQ = raw_header[:PG_start]
        raw_PG = raw_header[PG_start:CO_start]
        raw_CO = raw_header[CO_start:]
        assert raw_HD_and_SQ+raw_PG+raw_CO == raw_header

        new_PG_line = f"@PG\tID:{id}\tPN:{program}\tVN:{version}".encode('utf-8')

        if raw_PG:
            new_PG_line += (b"\tPP:"
                            + re.findall(rb'ID:[0-9a-zA-Z]+',raw_PG)[-1][3:])

        if command:
            new_PG_line += b"\tCL:" + command.encode('utf-8')

        if description:
            new_PG_line += b"\tDS:" + description.encode('utf-8')

        new_PG_line += b'\n'

        return self.update_header_length(raw_header= (raw_HD_and_SQ
                                                      + raw_PG
                                                      + new_PG_line
                                                      + raw_CO
                                                      )
                                         )

    def update_header(self, *args, **kwargs):
        """Use get_updated_header to update self.raw_header inplace

        See get_updated_header for Parameters and documentation
        """
        self.raw_header = self.get_updated_header(*args,**kwargs,
                                                  raw_header=self.raw_header)
        self.update_header_length()

class FileReader(_FileBase):
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
        header : str
            The ASCII representation of the header

        index_to_ref : Dict[int:str]
            A dictionary mapping bam reference numeric identifiers to names

        raw_header : bytes
            The raw bytestring representing the bam header

        raw_refs : bytes
            The raw bytestring representing the reference sequences

        refs : Dict[str:int]
            A dictionary of reference_name keys with reference_length

        ref_to_index : Dict[str:int]
            A dictionary mapping reference names to the bam numeric identifier

        sort_order : str
            The value of the SO field indicating sort type. Value is as given in
            the BAM file.
            Should be one of 'unknown', 'unsorted', 'queryname' or 'coordinate'

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
        self.sort_order = re.search(b'SO:[a-zA-Z]+',
                                   self.raw_header)[0][3:].decode()
        if self._ubam.seekable():
            self._start_of_alignments = self._ubam.tell()
        else: #pragma: no cover
            self._start_of_alignments = None
        self.alignments: Generator[bytes, None, None] = self._get_alignments()

    def __enter__(self):
        """Return self for use in WITH statement."""
        return self

    def __exit__(self, type, value, traceback):
        """Tidy up at end of WITH statement."""
        self._ubam.close()

    def close(self):
        """Close input file"""
        self._ubam.close()

    def _read_header(self) -> Tuple[bytes, str]:
        """Utility function for reading bam file headers
        Called by init to create self.raw_header and self.header"""
        raw_header_length = self._ubam.read(4)
        header_length = struct.unpack("<i", raw_header_length)[0]
        raw_header = raw_header_length + self._ubam.read(header_length)
        header = raw_header[4:].decode("latin-1")
        return (raw_header, header)

    def _read_refs(self) -> Tuple[bytes, Dict[str, int]]:
        """Utility function for reading bam file headers
        Called by init to create self.raw_refs and self.refs"""
        ref_buffer = b""
        refs = {}

        raw_nref = self._ubam.read(4)
        ref_buffer += raw_nref
        self.n_ref = struct.unpack("<i", raw_nref)[0]
        for i in range(self.n_ref):
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
        if self.n_ref != len(refs): #pragma: no cover
            raise ValueError(f"Invalid header: should be {self.n_ref} "
                             f"references but only found {len(refs)}"
                             )
        return (ref_buffer, refs)

    def _get_alignments(self) -> Generator[bytes, None, None]:
        """utility function to create an alignment generator"""
        while self._ubam:
            raw_blocksize = self._ubam.read(4)
            if not raw_blocksize:
                break
            block_size = struct.unpack("<i", raw_blocksize)[0]
            alignment = raw_blocksize + self._ubam.read(
                block_size
            )  # add back size so alignment record intact
            yield alignment

    def __iter__(self):
        return self

    def __next__(self):
        return next(self.alignments)

    def reset_alignments(self):
        """Reset the file pointer to the beginning of the alignment block"""
        if self._start_of_alignments:
            self._ubam.seek(self._start_of_alignments)
        else: #pragma: no cover
            raise NotImplementedError('Seek is not implemented for this file')


class FileWriter(_FileBase):
    def __init__(self,
                 file,
                 raw_header = None,
                 raw_refs = None,
                 mode = 'wb',
                 compresslevel = 6,
                 ):
        if not hasattr(file, 'write'):
            self.bgzf_file = BgzfWriter(filename=Path(file),
                                        mode=mode,
                                        fileobj=None,
                                        compresslevel=compresslevel,
                                        )
            self.name = str(Path(file)) #safe for Path or str
        else:
            self.bgzf_file = BgzfWriter(filename=None,
                                        fileobj=file,
                                        compresslevel=compresslevel,
                                        )
            self.name = file.name
        self.header_written = False
        self.magic = b"BAM\x01"
        self.raw_header = raw_header
        self.raw_refs = raw_refs

    def write(self, data):
        """
        Write output to the BAM file.
        Data is buffered by the underlying bgzf method

        Parameters
        ----------
        data : bytes
            The data to be written to the BAM file
        """
        return self.bgzf_file.write(data)

    def close(self, *args, **kwargs):
        """
        Flush and write any data to the BAM file before finalizing and closing
        """
        return self.bgzf_file.close(*args,**kwargs)

    def write_header(self,
                     raw_header = None,
                     raw_refs = None):
        """
        Write the header information to the BAM file

        Parameters
        ----------
        raw_header : bytes
            A raw byte format header containing the standard SAM format header
            Default : self.raw_header
        raw_refs
            A raw bytestring containing the reference sequence header data
            Default : self.raw_refs

        Raises
        -------
        RuntimeError
            Raises a runtime error if header information already written

        """
        if self.header_written:
            raise RuntimeError('Header has already been written to file')

        if not raw_header:
            raw_header = self.raw_header

        if not raw_refs:
            raw_refs = self.raw_refs

        self.update_header_length()
        self.write(self.magic + self.raw_header + self.raw_refs)
        self.header_written = True

    def __enter__(self):
        """Return self for use in WITH statement."""
        return self

    def __exit__(self, type, value, traceback):
        """Tidy up at end of WITH statement."""
        self.close()

    def tell(self):
        """Return the current location of the pointer in the file"""
        return self.bgzf_file.tell()

    def seekable(self):
        """Return the seek state of the file"""
        return False

