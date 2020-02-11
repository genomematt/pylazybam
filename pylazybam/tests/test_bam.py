#!/usr/bin/env python3
# encoding: utf-8
"""
pylazybam/tests/test_bam.py

Copyright (c) 2018-2020 Matthew Wakefield, The Walter and Eliza Hall Institute and The University of Melbourne. All rights reserved.
"""

import gzip, struct
import unittest

from pkg_resources import resource_stream
from tempfile import NamedTemporaryFile

from pylazybam import bam

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2018-2020 Matthew Wakefield, The Walter and Eliza Hall Institute and The University of Melbourne"
__credits__ = ["Matthew Wakefield", ]
__license__ = "BSD-3-Clause"
__version__ = "0.1.0"
__maintainer__ = "Matthew Wakefield"
__email__ = "wakefield@wehi.edu.au"
__status__ = "Development/Beta"

MIN32INT: int = -2147483648

RAW_HEADER = b'h\x02\x00\x00@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:MT\tLN:16569\n@SQ\tSN:1\tLN:249250621\n@SQ\tSN:2\tLN:243199373\n@SQ\tSN:3\tLN:198022430\n@SQ\tSN:4\tLN:191154276\n@SQ\tSN:5\tLN:180915260\n@SQ\tSN:6\tLN:171115067\n@SQ\tSN:7\tLN:159138663\n@SQ\tSN:8\tLN:146364022\n@SQ\tSN:9\tLN:141213431\n@SQ\tSN:10\tLN:135534747\n@SQ\tSN:11\tLN:135006516\n@SQ\tSN:12\tLN:133851895\n@SQ\tSN:13\tLN:115169878\n@SQ\tSN:14\tLN:107349540\n@SQ\tSN:15\tLN:102531392\n@SQ\tSN:16\tLN:90354753\n@SQ\tSN:17\tLN:81195210\n@SQ\tSN:18\tLN:78077248\n@SQ\tSN:19\tLN:59128983\n@SQ\tSN:20\tLN:63025520\n@SQ\tSN:21\tLN:48129895\n@SQ\tSN:22\tLN:51304566\n@SQ\tSN:X\tLN:155270560\n@SQ\tSN:Y\tLN:59373566\n@PG\tID:bowtie2\tPN:bowtie2\tVN:2.0.0-beta6\n'
RAW_REFS = b'\x19\x00\x00\x00\x03\x00\x00\x00MT\x00\xb9@\x00\x00\x02\x00\x00\x001\x00=C\xdb\x0e\x02\x00\x00\x002\x00\x8d\xed~\x0e\x02\x00\x00\x003\x00\x1e\x95\xcd\x0b\x02\x00\x00\x004\x00d\xc8d\x0b\x02\x00\x00\x005\x00<\x8c\xc8\n\x02\x00\x00\x006\x00;\x023\n\x02\x00\x00\x007\x00gC|\t\x02\x00\x00\x008\x00vV\xb9\x08\x02\x00\x00\x009\x00\xf7\xbej\x08\x03\x00\x00\x0010\x00\x9b\x18\x14\x08\x03\x00\x00\x0011\x004\t\x0c\x08\x03\x00\x00\x0012\x00\xf7j\xfa\x07\x03\x00\x00\x0013\x00VZ\xdd\x06\x03\x00\x00\x0014\x00$\x06f\x06\x03\x00\x00\x0015\x00@\x81\x1c\x06\x03\x00\x00\x0016\x00A\xb4b\x05\x03\x00\x00\x0017\x00\xca\xf0\xd6\x04\x03\x00\x00\x0018\x00@]\xa7\x04\x03\x00\x00\x0019\x00\x97<\x86\x03\x03\x00\x00\x0020\x00p\xb1\xc1\x03\x03\x00\x00\x0021\x00gg\xde\x02\x03\x00\x00\x0022\x00v\xd8\x0e\x03\x02\x00\x00\x00X\x00\xa0=A\t\x02\x00\x00\x00Y\x00\xfe\xf7\x89\x03'
REFS = {'MT': 16569, '1': 249250621, '2': 243199373, '3': 198022430, '4': 191154276, '5': 180915260, '6': 171115067,
        '7': 159138663, '8': 146364022, '9': 141213431, '10': 135534747, '11': 135006516, '12': 133851895,
        '13': 115169878, '14': 107349540, '15': 102531392, '16': 90354753, '17': 81195210, '18': 78077248,
        '19': 59128983, '20': 63025520, '21': 48129895, '22': 51304566, 'X': 155270560, 'Y': 59373566}
HEADER = '@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:MT\tLN:16569\n@SQ\tSN:1\tLN:249250621\n@SQ\tSN:2\tLN:243199373\n@SQ\tSN:3\tLN:198022430\n@SQ\tSN:4\tLN:191154276\n@SQ\tSN:5\tLN:180915260\n@SQ\tSN:6\tLN:171115067\n@SQ\tSN:7\tLN:159138663\n@SQ\tSN:8\tLN:146364022\n@SQ\tSN:9\tLN:141213431\n@SQ\tSN:10\tLN:135534747\n@SQ\tSN:11\tLN:135006516\n@SQ\tSN:12\tLN:133851895\n@SQ\tSN:13\tLN:115169878\n@SQ\tSN:14\tLN:107349540\n@SQ\tSN:15\tLN:102531392\n@SQ\tSN:16\tLN:90354753\n@SQ\tSN:17\tLN:81195210\n@SQ\tSN:18\tLN:78077248\n@SQ\tSN:19\tLN:59128983\n@SQ\tSN:20\tLN:63025520\n@SQ\tSN:21\tLN:48129895\n@SQ\tSN:22\tLN:51304566\n@SQ\tSN:X\tLN:155270560\n@SQ\tSN:Y\tLN:59373566\n@PG\tID:bowtie2\tPN:bowtie2\tVN:2.0.0-beta6\n'
INDEX_TO_REF = {0: 'MT', 1: '1', 2: '2', 3: '3', 4: '4', 5: '5', 6: '6', 7: '7', 8: '8', 9: '9', 10: '10', 11: '11',
                12: '12', 13: '13', 14: '14', 15: '15', 16: '16', 17: '17', 18: '18', 19: '19', 20: '20', 21: '21',
                22: '22', 23: 'X', 24: 'Y', -1: '*'}
REF_TO_INDEX = {'MT': 0, '1': 1, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7, '8': 8, '9': 9, '10': 10, '11': 11,
                '12': 12, '13': 13, '14': 14, '15': 15, '16': 16, '17': 17, '18': 18, '19': 19, '20': 20, '21': 21,
                '22': 22, 'X': 23, 'Y': 24}

ALIGN0 = b'\x12\x01\x00\x00\x0c\x00\x00\x00eB\xf0\x07(&\n2\x02\x00S\x00d\x00\x00\x00\x0c\x00\x00\x00\x1eB\xf0\x07S\xff\xff\xffHWI-ST960:96:COTO3ACXX:3:1101:1220:2089\x000\x06\x00\x00\x14\x00\x00\x00HB\x18"$H\x14\x82"\x14(\x12\x88\x14AD(AD!D\x14\x11\x82B\x88A\x12"\x84AD!\x11D\x88B\x14\x84\x14"AA\x82\x12\x12!(\x12\x1f#"##!\x1f####""!!\x1e##$$$##$"#"%%\'\'\'\'\')((&\')))))))(\'(()((\'#!))))))))((()))(&\'\'))))())()(&))(\'\'%%\'%####\x1c\x10\x02ASC\xc6XSC~XNC\x00XMC\x00XOC\x00XGC\x00NMC\x00MDZ99\x00YSC\xbdYTZCP\x00'
# HWI-ST960:96:COTO3ACXX:3:1101:1220:2089 83      12      133186150       38      99M1S   =       133186079       -173    GTGCATCCCGGTAGTCCCAGCTACTTAGGAGGCTGAGGCAGGAGAATCGCTTGAACCCTGGAGGCAAAGGTTGCAGTGAGCCGAGATCACACCACTACAN    DCDDB@DDDDCCBB?DDEEEDDECDCFFHHHHHJIIGHJJJJJJJIHIIJIIHDBJJJJJJJJIIIJJJIGHHJJJJIJJIJIGJJIHHFFHFDDDD=1#    AS:i:198        XS:i:126        XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:99 YS:i:189        YT:Z:CP
ALIGN42 = b'\xe4\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff(\x00H\x12\x00\x00M\x00d\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff\x00\x00\x00\x00HWI-ST960:96:COTO3ACXX:3:1101:1294:2176\x00\x81\x88\x11((\x14D\x14\x11!\x11\x14\x11\x11\x81\x12\x81A\x14!\x11\x11A\x14\x18!\x14B(\x14B"\x11H\x82\x88A\x12\x88\x88\x14\x82\x82\x18A"\x82\x84\x81\x81"""%%%%%\'\'\'\'\')\'(\'()))))(\'()))))))))()))())))))))))))())()(())()))\'\'$\'&\'%%%$!#$$$$$$#$#$$####"#######YTZUP\x00'
TAGS = b'SAZchr12,82544431,+,102M48S,60,3;\x00XAZchr12,+82544561,122S28M,0;\x00MCZ150M\x00PGZMarkDuplicates\x00ASc\x1eXSc\x1cMDZ0G30\x00NMI\x01\x00\x00\x00RGZNA12778_CTCACCAA-CTAGGCAA_HCKWTDSXX_L001\x00'



if __name__ == "__main__":
    unittest.main()


class test_main(unittest.TestCase):
    def setUp(self):
        pass

    def test_FileReader_init(self):
        test_bam = resource_stream(__name__, 'data/paired_end_testdata_human.bam')
        the_bam = bam.FileReader(gzip.open(test_bam))
        self.assertEqual(the_bam.raw_header, RAW_HEADER)
        self.assertEqual(the_bam.raw_refs, RAW_REFS)
        self.assertEqual(the_bam.refs, REFS)
        self.assertEqual(the_bam.header, HEADER)
        self.assertEqual(the_bam.index_to_ref, INDEX_TO_REF)
        self.assertEqual(the_bam.ref_to_index, REF_TO_INDEX)
        self.assertEqual(the_bam.sort_order, 'unsorted')
        the_bam.close()
        test_bam = resource_stream(__name__, 'data/paired_end_testdata_human.bam')
        with bam.FileReader(gzip.open(test_bam)) as the_bam:
            self.assertEqual(the_bam.header, HEADER)
            first_read = next(the_bam)
            the_bam.reset_alignments()
            self.assertEqual(first_read,next(the_bam))
        self.assertRaises(ValueError, bam.FileReader, test_bam)


    def test_update_header_length(self):
        test_bam = resource_stream(__name__, 'data/paired_end_testdata_human.bam')
        the_bam = bam.FileReader(gzip.open(test_bam))
        the_bam.raw_header = the_bam.raw_header+b'@CO\tThis is an extra comment\n'
        the_bam.update_header_length()
        self.assertEqual(len(the_bam.raw_header[4:]),
                         struct.unpack("<i", the_bam.raw_header[:4])[0])
        header_with_incorrect_len = struct.pack("<i",42)+RAW_HEADER[4:]
        self.assertEqual(the_bam.update_header_length(header_with_incorrect_len)[:4],
                         RAW_HEADER[:4])

    def test_get_updated_header(self):
        test_bam = resource_stream(__name__, 'data/paired_end_testdata_human.bam')
        the_bam = bam.FileReader(gzip.open(test_bam))
        the_bam.raw_header = the_bam.raw_header+b'@CO\tThis is an extra comment\n'
        the_bam.update_header_length()
        updated_header = the_bam.get_updated_header(id='test',
                                         program='test',
                                         version='88.88.8',
                                         command="--argumentative")
        self.assertEqual(updated_header[4:], #skip 4 bytes as len changed
                         (RAW_HEADER[4:]
                         + b'@PG\tID:test\tPN:test\tVN:88.88.8\tPP:bowtie2\t'
                         + b'CL:--argumentative\n'
                         + b'@CO\tThis is an extra comment\n')
                         )

    def test_FileReader_iteration(self):
        test_bam = resource_stream(__name__, 'data/paired_end_testdata_human.bam')
        the_bam = bam.FileReader(gzip.open(test_bam))
        self.assertEqual(next(the_bam), ALIGN0)
        for i in range(41):
            next(the_bam)
        self.assertEqual(next(the_bam), ALIGN42)
        for align in the_bam:
            pass


    def test_get_ref_index(self):
        self.assertEqual(bam.get_ref_index(ALIGN0), 12)

    def test_get_pos(self):
        self.assertEqual(bam.get_pos(ALIGN0), 133186149)

    def test_get_len_read_name(self):
        self.assertEqual(bam.get_len_read_name(ALIGN0), 40)

    def test_get_mapq(self):
        self.assertEqual(bam.get_mapq(ALIGN0), 38)

    def test_get_bin(self):
        self.assertEqual(bam.get_bin(ALIGN0), 12810)

    def test_get_number_cigar_operations(self):
        self.assertEqual(bam.get_number_cigar_operations(ALIGN0), 2)

    def test_get_flag(self):
        self.assertEqual(bam.get_flag(ALIGN0), 83)

    def test_get_len_sequence(self):
        self.assertEqual(bam.get_len_sequence(ALIGN0), 100)

    def test_get_pair_ref_index(self):
        self.assertEqual(bam.get_pair_ref_index(ALIGN0), 12)

    def test_get_pair_pos(self):
        self.assertEqual(bam.get_pair_pos(ALIGN0), 133186078)

    def test_get_template_len(self):
        self.assertEqual(bam.get_template_len(ALIGN0), -173)

    def test_get_read_name(self):
        self.assertEqual(bam.get_read_name(ALIGN0, 40), 'HWI-ST960:96:COTO3ACXX:3:1101:1220:2089')

    def test_get_raw_read_name(self):
        self.assertEqual(bam.get_raw_read_name(ALIGN0, 40), b'HWI-ST960:96:COTO3ACXX:3:1101:1220:2089\x00')

    def test_get_raw_cigar(self):
        self.assertEqual(bam.get_raw_cigar(ALIGN0, 40, 2), b'0\x06\x00\x00\x14\x00\x00\x00')

    def test_get_tag_bytestring(self):
        self.assertEqual(bam.get_tag_bytestring(ALIGN0, 40, 2, 100),
                         b'ASC\xc6XSC~XNC\x00XMC\x00XOC\x00XGC\x00NMC\x00MDZ99\x00YSC\xbdYTZCP\x00')

    def test_get_raw_sequence(self):
        self.assertEqual(bam.get_raw_sequence(ALIGN0, 40, 2, 100),
                         b'HB\x18"$H\x14\x82"\x14(\x12\x88\x14AD(AD!D\x14\x11\x82B\x88A\x12"\x84AD!\x11D\x88B\x14\x84\x14"AA\x82\x12\x12!(\x12\x1f')

    def test_get_raw_base_qual(self):
        self.assertEqual(bam.get_raw_base_qual(ALIGN0, 40, 2, 100),
                         b'#"##!\x1f####""!!\x1e##$$$##$"#"%%\'\'\'\'\')((&\')))))))(\'(()((\'#!))))))))((()))(&\'\'))))())()(&))(\'\'%%\'%####\x1c\x10\x02')

    def test_decode_base_qual(self):
        self.assertEqual(bam.decode_base_qual(
            b'#"##!\x1f####""!!\x1e##$$$##$"#"%%\'\'\'\'\')((&\')))))))(\'(()((\'#!))))))))((()))(&\'\'))))())()(&))(\'\'%%\'%####\x1c\x10\x02'),
                         'DCDDB@DDDDCCBB?DDEEEDDECDCFFHHHHHJIIGHJJJJJJJIHIIJIIHDBJJJJJJJJIIIJJJIGHHJJJJIJJIJIGJJIHHFFHFDDDD=1#')

    def test_decode_sequence(self):
        self.assertEqual(bam.decode_sequence(
            b'HB\x18"$H\x14\x82"\x14(\x12\x88\x14AD(AD!D\x14\x11\x82B\x88A\x12"\x84AD!\x11D\x88B\x14\x84\x14"AA\x82\x12\x12!(\x12\x1f'),
                         'GTGCATCCCGGTAGTCCCAGCTACTTAGGAGGCTGAGGCAGGAGAATCGCTTGAACCCTGGAGGCAAAGGTTGCAGTGAGCCGAGATCACACCACTACAN')

    def test_decode_cigar(self):
        self.assertEqual(bam.decode_cigar(b'0\x06\x00\x00\x14\x00\x00\x00'),
                         '99M1S')

    def test_get_AS(self):
        self.assertEqual(bam.get_AS(ALIGN0), 198)
        self.assertEqual(bam.get_AS(b''), MIN32INT)
        self.assertEqual(bam.get_AS(b'',no_tag=None), None)
        self.assertRaises(ValueError,bam.get_AS,b'ASC\xc6ASC\xc6')

    def test_get_XS(self):
        self.assertEqual(bam.get_XS(ALIGN0), 126)
        self.assertEqual(bam.get_XS(b''), MIN32INT)
        self.assertEqual(bam.get_XS(b'',no_tag=None), None)
        self.assertRaises(ValueError,bam.get_XS,b'XSC\xc6XSC\xc6')

    def test_get_ZS(self):
        self.assertEqual(bam.get_ZS(ALIGN0), MIN32INT)
        self.assertEqual(bam.get_ZS(b''), MIN32INT)
        self.assertEqual(bam.get_ZS(b'',no_tag=None), None)
        self.assertRaises(ValueError,bam.get_ZS,b'ZSC\xc6ZSC\xc6')
        self.assertEqual(bam.get_ZS(b'ASC\xc6ZSC~XNC\x00XMC\x00XOC\x00XGC\x00NMC\x00MDZ99\x00YSC\xbdYTZCP\x00'),
                         126)

    def test_get_MD(self):
        self.assertEqual(bam.get_MD(b''), None)
        self.assertEqual(bam.get_MD(b'',no_tag='SpanishInquisition'),
                         'SpanishInquisition')
        self.assertRaises(ValueError,bam.get_MD,b'MDZ99\x00MDZ99\x00')
        self.assertEqual(bam.get_MD(ALIGN0), '99')

    def test_get_int_tag(self):
        self.assertEqual(bam.get_int_tag(ALIGN0, b"AS"), 198)
        self.assertEqual(bam.get_int_tag(b'', b"AS"), MIN32INT)
        self.assertEqual(bam.get_int_tag(b'', b"AS", no_tag=None), None)
        self.assertRaises(ValueError, bam.get_int_tag,
                          *(b'ASC\xc6ASC\xc6', b"AS"))
        self.assertRaises(ValueError,bam.get_int_tag, *(b'',b''))

    def test_get_str_tag(self):
        self.assertEqual(bam.get_str_tag(b'',b'MD'), None)
        self.assertEqual(bam.get_str_tag(b'',b'MD',no_tag='SpanishInquisition'),
                         'SpanishInquisition')
        self.assertRaises(ValueError,bam.get_str_tag,
                          *(b'MDZ99\x00MDZ99\x00',b'MD'))
        self.assertEqual(bam.get_str_tag(ALIGN0,b'MD'), '99')
        self.assertRaises(ValueError,bam.get_str_tag, *(b'',b''))
        self.assertEqual(bam.get_str_tag(TAGS,b'SA'),
                         'chr12,82544431,+,102M48S,60,3;')
        self.assertEqual(bam.get_str_tag(TAGS,b'XA'),
                         'chr12,+82544561,122S28M,0;')
        self.assertEqual(bam.get_str_tag(TAGS,b'RG'),
                         'NA12778_CTCACCAA-CTAGGCAA_HCKWTDSXX_L001')
        self.assertEqual(bam.get_str_tag(TAGS,b'PG'),
                         'MarkDuplicates')

    def test_is_flag(self):
        self.assertTrue(bam.is_flag(ALIGN0,bam.FLAGS['forward']))

    def test_round_trip(self):
        test_bam = resource_stream(__name__, 'data/minitest.bam')
        # get raw uncompressed content
        content = gzip.open(test_bam).read()
        # go back to start of file and parse with Filereader
        test_bam.seek(0)
        the_bam = bam.FileReader(gzip.open(test_bam))
        out_file = NamedTemporaryFile(delete=False)
        out_file_name = out_file.name
        with bam.FileWriter(out_file_name) as out_bam:
            out_bam.raw_header = the_bam.raw_header
            out_bam.raw_refs = the_bam.raw_refs
            out_bam.write_header()
            self.assertRaises(RuntimeError,out_bam.write_header,None)
            for align in the_bam:
                out_bam.write(align)
        # parse the new temporary file
        new_bam = bam.FileReader(gzip.open(out_file_name))
        self.assertEqual(the_bam.header,new_bam.header)
        try:
            while True:
                align = next(the_bam)
                align2 = next(the_bam)
                self.assertEqual(align,align2)
        except StopIteration:
            pass

    def test_FileBase(selfs):
        the_filebase = bam._FileBase()
        # methods tested in classes that inherit


    def test_FileWriter(self):
        #construct with object
        outfile = NamedTemporaryFile(delete=True)
        out_bam = bam.FileWriter(outfile)
        self.assertEqual(out_bam.name,outfile.name)
        #construct with name
        outfile = NamedTemporaryFile(delete=False)
        outfile.close()
        out_bam = bam.FileWriter(outfile.name)
        self.assertEqual(out_bam.name, outfile.name)
        out_bam.raw_header = RAW_HEADER
        out_bam.raw_refs = RAW_REFS
        self.assertEqual(out_bam.get_full_raw_header(),
                         b'BAM\x01' + RAW_HEADER + RAW_REFS,
                         )
        out_bam.raw_header = RAW_HEADER + b"@CO\tComment\n"
        self.assertEqual(out_bam.get_full_raw_header()[8:],
                 RAW_HEADER[4:] + b"@CO\tComment\n" + RAW_REFS,
                 )
        self.assertEqual(len(out_bam.raw_header[4:]),
                         struct.unpack('<i',out_bam.raw_header[:4])[0])

        out_bam.update_header(id = 'pylazybam',
                              program = 'pylazybam',
                              version = '0.0.0',
                              description= 'testing123')
        self.assertEqual(out_bam.raw_header[4:],
                         (RAW_HEADER[4:]
                          + b'@PG\tID:pylazybam\tPN:pylazybam\tVN:0.0.0\t'
                          + b'PP:bowtie2\tDS:testing123\n@CO\tComment\n'
                          )
                         )
        self.assertFalse(out_bam.seekable())
        point = out_bam.tell()
        out_bam.write(ALIGN0)
        self.assertLess(point, out_bam.tell())


