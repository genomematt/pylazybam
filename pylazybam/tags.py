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
    The MD field aims to achieve SNP/indel calling independent of the reference.
    Numbers represent matches, letters bases that differ from the reference and
    bases preceeded by ^ are deletions. A ^T0A indicates a base change to an A
    immediately following a deleted T.

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

def get_int_tag(tag_bytes: bytes,
                tag: bytes,
                no_tag: Any = MIN32INT) -> int:
    """
    Extract an integer format tag from a raw BAM alignment bytestring

    Parameters
    ----------
        tag_bytes : bytes
            a bytestring containing bam formatted tag elements

        tag : bytes
            the two byte tag to be returned eg b'AS'
            The tag must represent an integer format tag

        no_tag : Any
            return value for when tag not found (default: MIN32INT)

    Returns
    -------
        int
            the integer value of the tag
            returns the value of no_tag if tag absent (default:MIN32INT)

    Raises
    ------
        ValueError
            raises a ValueError if more than one tag match

    Notes:
    ------
    Potential values for the tag parameter include:

        AM:i:score  The smallest template-independent mapping quality of any
                    segment in the same template as this read. (See also SM.)
        AS:i:score  Alignment score generated by aligner.

        CP:i:pos    Leftmost coordinate of the next hit.

        FI:i:int    The index of segment in the template.

        H0:i:count  Number of perfect hits.

        H1:i:count  Number of 1-difference hits (see also NM).

        H2:i:count  Number of 2-difference hits.

        HI:i:i      Query hit index, indicating the alignment record is the
                    i-th one stored in SAM.

        IH:i:count  Number of alignments stored in the file that contain the
                    query in the current record.

        MQ:i:score  Mapping quality of the mate/next segment.

        NH:i:count  Number of reported alignments that contain the query in
                    the current record.

        NM:i:count  Number of differences (mismatches plus inserted and deleted
                    bases) between the sequence and reference, counting only
                    (case-insensitive) A, C, G and T bases in sequence and
                    reference as potential matches, with everything else being
                    a mismatch.

        PQ:i:score  Phred likelihood of the template, conditional on the mapping
                    locations of both/all segments being correct.

        SM:i:score  Template-independent mapping quality, i.e., the mapping
                    quality if the read were mapped as a single read rather
                    than as part of a read pair or template.

        TC:i:       The number of segments in the template.

        UQ:i:       Phred likelihood of the segment,
                    conditional on the mapping being correct.

        see https://samtools.github.io/hts-specs/SAMtags.pdf

    """
    if len(tag) != 2 or type(tag) != bytes:
        raise ValueError(f"Tags must be two bytes not {tag}")
    match = re.findall(tag + b"C.", tag_bytes)
    if not match:
        return no_tag
    elif len(match) != 1:
        raise ValueError(
            (
                f"More than one match to {tag+b'C.'} was found in {tag_bytes}"
                "meaning that more than one value of the tag is present"
                "or that another tag contains a match as part of its value"
            )
        )
    else:
        return struct.unpack("<xxxB", match[0])[-1]

def get_str_tag(tag_bytes: bytes,
                tag: bytes,
                no_tag: Any = None) -> str:
    """Extract an integer format tag from a raw BAM alignment bytestring

    Parameters
    ----------
        tag_bytes : bytes
            a bytestring containing bam formatted tag elements

        tag : bytes
            the two byte tag to be returned eg b'AS'
            The tag must represent an integer format tag

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
    Potential values for the tag parameter include:

        BQ:Z:qualities  Offset to base alignment quality (BAQ), of the same
                        length as the read sequence. At the i-th read base,
                        BAQi = Qi − (BQi − 64) where Qi is the i-th base quality

        CC:Z:rname      Reference name of the next hit; ‘=’ for same chromosome.

        E2:Z:bases      The 2nd most likely base calls.
                        Same encoding and same length as SEQ.

        FS:Z:str        Segment suffix.

        MC:Z:cigar      CIGAR string for mate/next segment.

        MD:Z:           String for mismatching positions.

        Q2:Z:qualities  Phred quality of the mate/next segment sequence in the
                        R2 tag. Same encoding as QUAL.

        R2:Z:bases      Sequence of the mate/next segment in the template.

        SA:Z:           (rname ,pos ,strand ,CIGAR ,mapQ ,NM ;)+
                        Other canonical alignments in a chimeric alignment,
                        formatted as a semicolon-delimited list.

        U2:Z:           Phred probability of the 2nd call being wrong
                        conditional on the best being wrong.

        RG:Z:readgroup  The read group to which the read belongs.

        LB:Z:library    The library from which the read has been sequenced.

        PG:Z:program id Program. Value matches the header PG-ID tag

        PU:Z:platformunit The platform unit in which the read was sequenced.

        CO:Z:text       Free-text comments.

        BC:Z:sequence   Barcode sequence (Identifying the sample/library),
                        with any quality scores (optionally) stored in  QT tag.
                        The BC tag should match the QT tag in length.

        QT:Z:qualities  Phred quality of the sample barcode sequence in BC tag.
                        Same encoding as QUAL, i.e., Phred score + 33.

        CB:Z:str        Cell identifier, consisting of the optionally-corrected
                        cellular barcode sequence and an optional suffix.
                        The sequence part is similar to the CR tag

        CR:Z:sequence+  Cellular barcode. The uncorrected sequence bases of the
                        cellular barcode as reported by the sequencing machine,
                        with the corresponding base quality scores (optionally)
                        stored in CY.

        CY:Z:qualities+ Phred quality of the cellular barcode sequence in CR tag
                        Same encoding as QUAL, i.e., Phred score + 33.

        MI:Z:str        Molecular Identifier. A unique ID within the SAM file
                        for the source molecule from which this read is derived.

        OX:Z:sequence+  Raw (uncorrected) unique molecular identifier bases,
                        with  quality scores (optionally) stored in the BZ tag.

        BZ:Z:qualities+ Phred quality of the (uncorrected) unique molecular
                        identifier sequence in the OX tag.
                        Same encoding as QUAL, i.e., Phred score + 33.

        RX:Z:sequence+  Sequence bases from the unique molecular identifier.
                        These could be either corrected or uncorrected.
                        Unlike MI, the value may be non-unique in the file.

        QX:Z:qualities+ Phred quality of the unique molecular identifier
                        sequence in the RX tag.
                        Same encoding as QUAL, i.e., Phred score + 33

        OA:Z:(RNAME,POS,strand,CIGAR,MAPQ,NM;)+ The original alignment
                        information of the record prior to realignment or
                        unalignment by a subsequent tool.

        OC:Z:cigar      Original CIGAR, usually before realignment.

        CT:Z:strand ;type (;key (=value )?)* Complete read annotation tag
                        used for consensus annotation dummy features.

        PT:Z:annotag(\\|annotag)* where each annotag matches
                        start;end;strand;type(;key(=value)?)*
                        Read annotations for parts of the padded read sequence.

    see https://samtools.github.io/hts-specs/SAMtags.pdf

    Recommended try accept for use on raw alignment with fall back
    to calling on only the tag byte string.

    Please test carefully on your BAM output as in complicated output the
    regular expression based extraction of the tag can be error prone

    """
    if len(tag) != 2 or type(tag) != bytes:
        raise ValueError(f"Tags must be two bytes not {tag}")
    match = re.findall(tag+b"Z[^\x00]+\x00", tag_bytes)
    if not match:
        return no_tag
    elif len(match) != 1:
        restring = tag+b'Z[+\\x00]+\\x00'
        raise ValueError(
            (f"More than one match to {restring} was found in "
             f"{tag_bytes} "
             "meaning that more than one value of the tag is present"
             "or that another tag contains a match as part of its value")
        )
    else:
        return match[0][3:-1].decode()
