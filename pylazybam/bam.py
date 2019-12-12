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
import sys, os
import struct
from array import array
import re
from pylazybam.bgzf import BgzfWriter as FileWriter

# Parsing functions

def get_flag_as_int(alignment):
    return struct.unpack('<H',alignment[18:20])[0]

def get_read_name_length(alignment):
    return struct.unpack('<B',alignment[12:13])[0]

def get_read_name_text(alignment, read_name_length):
    return alignment[36:36+read_name_length].decode("utf-8")

def get_read_name_raw(alignment, read_name_length):
    return alignment[36:36+read_name_length]

def get_number_of_cigar_operations(alignment):
    return struct.unpack('<H',alignment[16:18])[0]

def get_sequence_length(alignment):
    return struct.unpack('<i',alignment[20:24])[0]

def get_tag_raw(alignment, len_read_name, n_cigar_op, len_sequence):
    start = 36 + len_read_name + (4*n_cigar_op) +  (len_sequence+1)//2 + len_sequence
    return alignment[start:]

def get_sequence_raw(alignment, len_read_name, n_cigar_op, len_sequence):
    start = 36 + len_read_name + (4*n_cigar_op)
    end = start + (len_sequence+1)//2
    return alignment[start:end+1]

def get_ref(alignment):
    return

def get_pos(alignment):
    return

def quick_AS_XS(tag_raw):
    if (tag_raw[4:7] == b'XSC' or tag_raw[4:7] == b'ZSC') and tag_raw[:3] == b'ASC':
        return (struct.unpack('<B',tag_raw[3:4])[0],struct.unpack('<B',tag_raw[8:9])[0])
    else:
        raise TypeError

def decode_sequence(seq_raw):
    bam_bases = '=ACMGRSVTWYHKDBN'
    result = ''
    for x in array('B', seq_raw):
        result += bam_bases[x >> 4]
        result += bam_bases[x & 0b1111]
    return result
    
def get_sequence_text(alignment, len_read_name, n_cigar_op, len_sequence):
    return decode_sequences(get_sequence_raw(alignment,
                                             len_read_name,
                                             n_cigar_op,
                                             len_sequence)

class Reader():
    """A Pure Python Lazy Bam Parser Class
    """
    def __init__(self, ubam):
        self.ubam = ubam
        self.magic = self.ubam.read(4)
        assert self.magic == b'BAM\x01' ###TODO add more informative errors
        self.raw_header, self.ascii_header = self.read_header()
        self.raw_refs, self.refs = self.read_refs()
        self.alignments = self._get_alignments()

    def __enter__(self):
        """Open a file operable with WITH statement."""
        return self
    
    def __exit__(self, type, value, traceback):
        """Close a file operable with WITH statement."""
        self.ubam.close()
    
    def close(self):
        """Close input file"""
        return self.ubam.close()

    def read_header(self):
        raw_header_length = self.ubam.read(4)
        header_length = struct.unpack('<i',raw_header_length)[0]
        raw_header = raw_header_length+self.ubam.read(header_length)
        header = raw_header.decode("utf-8")
        return (raw_header, header)
 
    def read_refs(self):
        ref_buffer = b""
        refs = {}

        raw_nref = self.ubam.read(4)
        ref_buffer += raw_nref
        n_ref = struct.unpack('<i',raw_nref)[0]
        for i in range(n_ref):
            raw_lname = self.ubam.read(4)
            ref_buffer += raw_lname
            l_name = struct.unpack('<i',raw_lname)[0]
            raw_ref_name = self.ubam.read(l_name)
            ref_buffer += raw_ref_name
            ref_name = raw_ref_name.decode("utf-8")
            raw_ref_length = self.ubam.read(4)
            ref_buffer += raw_ref_length
            ref_length = struct.unpack('<i',raw_ref_length)[0]
            refs[ref_name] = ref_length
        return (ref_buffer,refs)
    
    def _get_alignments(self):
        while self.ubam:
            try:
                raw_blocksize = self.ubam.read(4)
                block_size = struct.unpack('<i',raw_blocksize)[0]
                alignment = raw_blocksize+self.ubam.read(block_size) #add back size so alignment record intact
                yield alignment
            except:
                if not self.ubam.read():
                    break
    
    def __iterator__(self):
        return self
    
    def __next__(self):
        return next(self.alignments)