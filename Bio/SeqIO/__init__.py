# Copyright 2006-2010 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
#Nice link:
# http://www.ebi.ac.uk/help/formats_frame.html

"""Sequence input/output as SeqRecord objects.

Bio.SeqIO is also documented at U{http://biopython.org/wiki/SeqIO} and by
a whole chapter in our tutorial:
 - U{http://biopython.org/DIST/docs/tutorial/Tutorial.html}
 - U{http://biopython.org/DIST/docs/tutorial/Tutorial.pdf}  

Input
=====
The main function is Bio.SeqIO.parse(...) which takes an input file handle
(or in recent versions of Biopython alternatively a filename as a string),
and format string.  This returns an iterator giving SeqRecord objects:

    >>> from Bio import SeqIO
    >>> for record in SeqIO.parse("Fasta/f002", "fasta"):
    ...     print record.id, len(record)
    gi|1348912|gb|G26680|G26680 633
    gi|1348917|gb|G26685|G26685 413
    gi|1592936|gb|G29385|G29385 471

Note that the parse() function will invoke the relevant parser for the
format with its default settings.  You may want more control, in which case
you need to create a format specific sequence iterator directly.

Input - Single Records
======================
If you expect your file to contain one-and-only-one record, then we provide
the following 'helper' function which will return a single SeqRecord, or
raise an exception if there are no records or more than one record:

    >>> from Bio import SeqIO
    >>> record = SeqIO.read("Fasta/f001", "fasta")
    >>> print record.id, len(record)
    gi|3318709|pdb|1A91| 79

This style is useful when you expect a single record only (and would
consider multiple records an error).  For example, when dealing with GenBank
files for bacterial genomes or chromosomes, there is normally only a single
record.  Alternatively, use this with a handle when downloading a single
record from the internet.

However, if you just want the first record from a file containing multiple
record, use the iterator's next() method:

    >>> from Bio import SeqIO
    >>> record = SeqIO.parse("Fasta/f002", "fasta").next()
    >>> print record.id, len(record)
    gi|1348912|gb|G26680|G26680 633

The above code will work as long as the file contains at least one record.
Note that if there is more than one record, the remaining records will be
silently ignored.


Input - Multiple Records
========================
For non-interlaced files (e.g. Fasta, GenBank, EMBL) with multiple records
using a sequence iterator can save you a lot of memory (RAM).  There is
less benefit for interlaced file formats (e.g. most multiple alignment file
formats).  However, an iterator only lets you access the records one by one.

If you want random access to the records by number, turn this into a list:

    >>> from Bio import SeqIO
    >>> records = list(SeqIO.parse("Fasta/f002", "fasta"))
    >>> len(records)
    3
    >>> print records[1].id
    gi|1348917|gb|G26685|G26685

If you want random access to the records by a key such as the record id,
turn the iterator into a dictionary:

    >>> from Bio import SeqIO
    >>> record_dict = SeqIO.to_dict(SeqIO.parse("Fasta/f002", "fasta"))
    >>> len(record_dict)
    3
    >>> print len(record_dict["gi|1348917|gb|G26685|G26685"])
    413

However, using list() or the to_dict() function will load all the records
into memory at once, and therefore is not possible on very large files.
Instead, for *some* file formats Bio.SeqIO provides an indexing approach
providing dictionary like access to any record. For example,

    >>> from Bio import SeqIO
    >>> record_dict = SeqIO.index("Fasta/f002", "fasta")
    >>> len(record_dict)
    3
    >>> print len(record_dict["gi|1348917|gb|G26685|G26685"])
    413

Many but not all of the supported input file formats can be indexed like
this. For example "fasta", "fastq", "qual" and even the binary format "sff"
work, but alignment formats like "phylip", "clustalw" and "nexus" will not.

In most cases you can also use SeqIO.index to get the record from the file
as a raw string (not a SeqRecord). This can be useful for example to extract
a sub-set of records from a file where SeqIO cannot output the file format
(e.g. the plain text SwissProt format, "swiss") or where it is important to
keep the output 100% identical to the input). For example,

    >>> from Bio import SeqIO
    >>> record_dict = SeqIO.index("Fasta/f002", "fasta")
    >>> len(record_dict)
    3
    >>> print record_dict.get_raw("gi|1348917|gb|G26685|G26685")
    >gi|1348917|gb|G26685|G26685 human STS STS_D11734.
    CGGAGCCAGCGAGCATATGCTGCATGAGGACCTTTCTATCTTACATTATGGCTGGGAATCTTACTCTTTC
    ATCTGATACCTTGTTCAGATTTCAAAATAGTTGTAGCCTTATCCTGGTTTTACAGATGTGAAACTTTCAA
    GAGATTTACTGACTTTCCTAGAATAGTTTCTCTACTGGAAACCTGATGCTTTTATAAGCCATTGTGATTA
    GGATGACTGTTACAGGCTTAGCTTTGTGTGAAANCCAGTCACCTTTCTCCTAGGTAATGAGTAGTGCTGT
    TCATATTACTNTAAGTTCTATAGCATACTTGCNATCCTTTANCCATGCTTATCATANGTACCATTTGAGG
    AATTGNTTTGCCCTTTTGGGTTTNTTNTTGGTAAANNNTTCCCGGGTGGGGGNGGTNNNGAAA
    <BLANKLINE>
    >>> print record_dict["gi|1348917|gb|G26685|G26685"].format("fasta")
    >gi|1348917|gb|G26685|G26685 human STS STS_D11734.
    CGGAGCCAGCGAGCATATGCTGCATGAGGACCTTTCTATCTTACATTATGGCTGGGAATC
    TTACTCTTTCATCTGATACCTTGTTCAGATTTCAAAATAGTTGTAGCCTTATCCTGGTTT
    TACAGATGTGAAACTTTCAAGAGATTTACTGACTTTCCTAGAATAGTTTCTCTACTGGAA
    ACCTGATGCTTTTATAAGCCATTGTGATTAGGATGACTGTTACAGGCTTAGCTTTGTGTG
    AAANCCAGTCACCTTTCTCCTAGGTAATGAGTAGTGCTGTTCATATTACTNTAAGTTCTA
    TAGCATACTTGCNATCCTTTANCCATGCTTATCATANGTACCATTTGAGGAATTGNTTTG
    CCCTTTTGGGTTTNTTNTTGGTAAANNNTTCCCGGGTGGGGGNGGTNNNGAAA
    <BLANKLINE>

Here the original file and what Biopython would output differ in the line
wrapping.

Input - Alignments
==================
You can read in alignment files as alignment objects using Bio.AlignIO.
Alternatively, reading in an alignment file format via Bio.SeqIO will give
you a SeqRecord for each row of each alignment:

    >>> from Bio import SeqIO
    >>> for record in SeqIO.parse("Clustalw/hedgehog.aln", "clustal"):
    ...     print record.id, len(record)
    gi|167877390|gb|EDS40773.1| 447
    gi|167234445|ref|NP_001107837. 447
    gi|74100009|gb|AAZ99217.1| 447
    gi|13990994|dbj|BAA33523.2| 447
    gi|56122354|gb|AAV74328.1| 447

Output
======
Use the function Bio.SeqIO.write(...), which takes a complete set of
SeqRecord objects (either as a list, or an iterator), an output file handle
(or in recent versions of Biopython an output filename as a string) and of
course the file format::

    from Bio import SeqIO
    records = ...
    SeqIO.write(records, "example.faa", "fasta")

Or, using a handle::

    from Bio import SeqIO
    records = ...
    handle = open("example.faa", "w")
    SeqIO.write(records, handle, "fasta")
    handle.close()

You are expected to call this function once (with all your records) and if
using a handle, make sure you close it to flush the data to the hard disk.

Output - Advanced
=================
The effect of calling write() multiple times on a single file will vary
depending on the file format, and is best avoided unless you have a strong
reason to do so.

If you give a filename, then each time you call write() the existing file
will be overwriten. For sequential files formats (e.g. fasta, genbank) each
"record block" holds a single sequence.  For these files it would probably
be safe to call write() multiple times by re-using the same handle.


However, trying this for certain alignment formats (e.g. phylip, clustal,
stockholm) would have the effect of concatenating several multiple sequence
alignments together.  Such files are created by the PHYLIP suite of programs
for bootstrap analysis, but it is clearer to do this via Bio.AlignIO instead.


Conversion
==========
The Bio.SeqIO.convert(...) function allows an easy interface for simple
file format conversions. Additionally, it may use file format specific
optimisations so this should be the fastest way too.

In general however, you can combine the Bio.SeqIO.parse(...) function with
the Bio.SeqIO.write(...) function for sequence file conversion. Using
generator expressions or generator functions provides a memory efficient way
to perform filtering or other extra operations as part of the process.

File Formats
============
When specifying the file format, use lowercase strings.  The same format
names are also used in Bio.AlignIO and include the following:

 - ace     - Reads the contig sequences from an ACE assembly file.
 - embl    - The EMBL flat file format. Uses Bio.GenBank internally.
 - fasta   - The generic sequence file format where each record starts with
             an identifer line starting with a ">" character, followed by
             lines of sequence.
 - fastq   - A "FASTA like" format used by Sanger which also stores PHRED
             sequence quality values (with an ASCII offset of 33).
 - fastq-sanger - An alias for "fastq" for consistency with BioPerl and EMBOSS
 - fastq-solexa - Original Solexa/Illumnia variant of the FASTQ format which
             encodes Solexa quality scores (not PHRED quality scores) with an
             ASCII offset of 64.
 - fastq-illumina - Solexa/Illumnia 1.3+ variant of the FASTQ format which
             encodes PHRED quality scores with an ASCII offset of 64 (not 33).
 - genbank - The GenBank or GenPept flat file format.
 - gb      - An alias for "genbank", for consistency with NCBI Entrez Utilities
 - ig      - The IntelliGenetics file format, apparently the same as the
             MASE alignment format.
 - imgt    - An EMBL like format from IMGT where the feature tables are more
             indented to allow for longer feature types.
 - phd     - Output from PHRED, used by PHRAP and CONSED for input.
 - pir     - A "FASTA like" format introduced by the National Biomedical
             Research Foundation (NBRF) for the Protein Information Resource
             (PIR) database, now part of UniProt.
 - sff     - Standard Flowgram Format (SFF), typical output from Roche 454.
 - sff-trim - Standard Flowgram Format (SFF) with given trimming applied.
 - swiss   - Plain text Swiss-Prot aka UniProt format.
 - tab     - Simple two column tab separated sequence files, where each
             line holds a record's identifier and sequence. For example,
             this is used as by Aligent's eArray software when saving
             microarray probes in a minimal tab delimited text file.
 - qual    - A "FASTA like" format holding PHRED quality values from
             sequencing DNA, but no actual sequences (usually provided
             in separate FASTA files).
 - uniprot-xml - The UniProt XML format (replacement for the SwissProt plain
             text format which we call "swiss")

Note that while Bio.SeqIO can read all the above file formats, it cannot
write to all of them.

You can also use any file format supported by Bio.AlignIO, such as "nexus",
"phlip" and "stockholm", which gives you access to the individual sequences
making up each alignment as SeqRecords.
"""
__docformat__ = "epytext en" #not just plaintext

#TODO
# - define policy on reading aligned sequences with gaps in
#   (e.g. - and . characters) including how the alphabet interacts
#
# - How best to handle unique/non unique record.id when writing.
#   For most file formats reading such files is fine; The stockholm
#   parser would fail.
#
# - MSF multiple alignment format, aka GCG, aka PileUp format (*.msf)
#   http://www.bioperl.org/wiki/MSF_multiple_alignment_format 

"""
FAO BioPython Developers
========================
The way I envision this SeqIO system working as that for any sequence file
format we have an iterator that returns SeqRecord objects.

This also applies to interlaced fileformats (like clustal - although that
is now handled via Bio.AlignIO instead) where the file cannot be read record
by record.  You should still return an iterator, even if the implementation
could just as easily return a list.

These file format specific sequence iterators may be implemented as:
* Classes which take a handle for __init__ and provide the __iter__ method
* Functions that take a handle, and return an iterator object
* Generator functions that take a handle, and yield SeqRecord objects

It is then trivial to turn this iterator into a list of SeqRecord objects,
an in memory dictionary, or a multiple sequence alignment object.

For building the dictionary by default the id propery of each SeqRecord is
used as the key.  You should always populate the id property, and it should
be unique in most cases. For some file formats the accession number is a good
choice.  If the file itself contains ambiguous identifiers, don't try and
dis-ambiguate them - return them as is.

When adding a new file format, please use the same lower case format name
as BioPerl, or if they have not defined one, try the names used by EMBOSS.

See also http://biopython.org/wiki/SeqIO_dev

--Peter
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Generic import Alignment
from Bio.Alphabet import Alphabet, AlphabetEncoder, _get_base_alphabet

import FastaIO

#Convention for format names is "mainname-subtype" in lower case.
#Please use the same names as BioPerl or EMBOSS where possible.
#
#Note that this simple system copes with defining
#multiple possible iterators for a given format/extension
#with the -subtype suffix
#
#Most alignment file formats will be handled via Bio.AlignIO

_FormatToIterator = {"fasta" : FastaIO.FastaIterator,
                     }

_FormatToWriter = {"fasta" : FastaIO.FastaWriter,
                   }

_BinaryFormats = ["sff", "sff-trim"]

def write(sequences, handle, format):
    """Write complete set of sequences to a file.

     - sequences - A list (or iterator) of SeqRecord objects, or (if using
                   Biopython 1.54 or later) a single SeqRecord.
     - handle    - File handle object to write to, or filename as string
                   (note older versions of Biopython only took a handle).
     - format    - lower case string describing the file format to write.

    You should close the handle after calling this function.

    Returns the number of records written (as an integer).
    """
    from Bio import AlignIO

    #Try and give helpful error messages:
    if not isinstance(format, basestring):
        raise TypeError("Need a string for the file format (lower case)")
    if not format:
        raise ValueError("Format required (lower case string)")
    if format != format.lower():
        raise ValueError("Format string '%s' should be lower case" % format)

    if isinstance(sequences, SeqRecord):
        #This raised an exception in order version of Biopython
        sequences = [sequences]

    if isinstance(handle, basestring):
        if format in _BinaryFormats :
            handle = open(handle, "wb")
        else :
            handle = open(handle, "w")
        handle_close = True
    else:
        handle_close = False

    #Map the file format to a writer class
    if format in _FormatToWriter:
        writer_class = _FormatToWriter[format]
        count = writer_class(handle).write_file(sequences)
    elif format in AlignIO._FormatToWriter:
        #Try and turn all the records into a single alignment,
        #and write that using Bio.AlignIO
        alignment = MultipleSeqAlignment(sequences)
        alignment_count = AlignIO.write([alignment], handle, format)
        assert alignment_count == 1, "Internal error - the underlying writer " \
           + " should have returned 1, not %s" % repr(alignment_count)
        count = len(alignment)
        del alignment_count, alignment
    elif format in _FormatToIterator or format in AlignIO._FormatToIterator:
        raise ValueError("Reading format '%s' is supported, but not writing" \
                         % format)
    else:
        raise ValueError("Unknown format '%s'" % format)

    assert isinstance(count, int), "Internal error - the underlying %s " \
           "writer should have returned the record count, not %s" \
           % (format, repr(count))

    if handle_close:
        handle.close()
    
    return count
    
def parse(handle, format, alphabet=None):
    r"""Turns a sequence file into an iterator returning SeqRecords.

     - handle   - handle to the file, or the filename as a string
                  (note older verions of Biopython only took a handle).
     - format   - lower case string describing the file format.
     - alphabet - optional Alphabet object, useful when the sequence type
                  cannot be automatically inferred from the file itself
                  (e.g. format="fasta" or "tab")

    Typical usage, opening a file to read in, and looping over the record(s):

    >>> from Bio import SeqIO
    >>> filename = "Fasta/sweetpea.nu"
    >>> for record in SeqIO.parse(filename, "fasta"):
    ...    print "ID", record.id
    ...    print "Sequence length", len(record)
    ...    print "Sequence alphabet", record.seq.alphabet
    ID gi|3176602|gb|U78617.1|LOU78617
    Sequence length 309
    Sequence alphabet SingleLetterAlphabet()

    For file formats like FASTA where the alphabet cannot be determined, it
    may be useful to specify the alphabet explicitly:

    >>> from Bio import SeqIO
    >>> from Bio.Alphabet import generic_dna
    >>> filename = "Fasta/sweetpea.nu"
    >>> for record in SeqIO.parse(filename, "fasta", generic_dna):
    ...    print "ID", record.id
    ...    print "Sequence length", len(record)
    ...    print "Sequence alphabet", record.seq.alphabet
    ID gi|3176602|gb|U78617.1|LOU78617
    Sequence length 309
    Sequence alphabet DNAAlphabet()

    If you have a string 'data' containing the file contents, you must
    first turn this into a handle in order to parse it:

    >>> data = ">Alpha\nACCGGATGTA\n>Beta\nAGGCTCGGTTA\n"
    >>> from Bio import SeqIO
    >>> from StringIO import StringIO
    >>> for record in SeqIO.parse(StringIO(data), "fasta"):
    ...     print record.id, record.seq
    Alpha ACCGGATGTA
    Beta AGGCTCGGTTA

    Use the Bio.SeqIO.read(...) function when you expect a single record
    only.
    """
    #NOTE - The above docstring has some raw \n characters needed
    #for the StringIO example, hense the whole docstring is in raw
    #string mode (see the leading r before the opening quote).
    from Bio import AlignIO

    handle_close = False
    
    if isinstance(handle, basestring):
        #Hack for SFF, will need to make this more general in future
        if format in _BinaryFormats :
            handle = open(handle, "rb")
        else :
            handle = open(handle, "rU")
        #TODO - On Python 2.5+ use with statement to close handle
        handle_close = True

    #Try and give helpful error messages:
    if not isinstance(format, basestring):
        raise TypeError("Need a string for the file format (lower case)")
    if not format:
        raise ValueError("Format required (lower case string)")
    if format != format.lower():
        raise ValueError("Format string '%s' should be lower case" % format)
    if alphabet is not None and not (isinstance(alphabet, Alphabet) or \
                                     isinstance(alphabet, AlphabetEncoder)):
        raise ValueError("Invalid alphabet, %s" % repr(alphabet))

    #Map the file format to a sequence iterator:    
    if format in _FormatToIterator:
        iterator_generator = _FormatToIterator[format]
        if alphabet is None:
            i = iterator_generator(handle)
        else:
            try:
                i = iterator_generator(handle, alphabet=alphabet)
            except TypeError:
                i = _force_alphabet(iterator_generator(handle), alphabet)
    elif format in AlignIO._FormatToIterator:
        #Use Bio.AlignIO to read in the alignments
        #TODO - Can this helper function can be replaced with a generator
        #expression, or something from itertools?
        i = _iterate_via_AlignIO(handle, format, alphabet)
    else:
        raise ValueError("Unknown format '%s'" % format)
    #This imposes some overhead... wait until we drop Python 2.4 to fix it
    for r in i:
        yield r
    if handle_close:
        handle.close()

#This is a generator function
def _iterate_via_AlignIO(handle, format, alphabet):
    """Iterate over all records in several alignments (PRIVATE)."""
    from Bio import AlignIO
    for align in AlignIO.parse(handle, format, alphabet=alphabet):
        for record in align:
            yield record

def _force_alphabet(record_iterator, alphabet):
    """Iterate over records, over-riding the alphabet (PRIVATE)."""
    #Assume the alphabet argument has been pre-validated
    given_base_class = _get_base_alphabet(alphabet).__class__
    for record in record_iterator:
        if isinstance(_get_base_alphabet(record.seq.alphabet),
                      given_base_class):
            record.seq.alphabet = alphabet
            yield record
        else:
            raise ValueError("Specified alphabet %s clashes with "\
                             "that determined from the file, %s" \
                             % (repr(alphabet), repr(record.seq.alphabet)))

def read(handle, format, alphabet=None):
    """Turns a sequence file into a single SeqRecord.

     - handle   - handle to the file, or the filename as a string
                  (note older verions of Biopython only took a handle).
     - format   - string describing the file format.
     - alphabet - optional Alphabet object, useful when the sequence type
                  cannot be automatically inferred from the file itself
                  (e.g. format="fasta" or "tab")

    This function is for use parsing sequence files containing
    exactly one record.  For example, reading a GenBank file:

    >>> from Bio import SeqIO
    >>> record = SeqIO.read("GenBank/arab1.gb", "genbank")
    >>> print "ID", record.id
    ID AC007323.5
    >>> print "Sequence length", len(record)
    Sequence length 86436
    >>> print "Sequence alphabet", record.seq.alphabet
    Sequence alphabet IUPACAmbiguousDNA()

    If the handle contains no records, or more than one record,
    an exception is raised.  For example:

    >>> from Bio import SeqIO
    >>> record = SeqIO.read("GenBank/cor6_6.gb", "genbank")
    Traceback (most recent call last):
        ...
    ValueError: More than one record found in handle

    If however you want the first record from a file containing
    multiple records this function would raise an exception (as
    shown in the example above).  Instead use:

    >>> from Bio import SeqIO
    >>> record = SeqIO.parse("GenBank/cor6_6.gb", "genbank").next()
    >>> print "First record's ID", record.id
    First record's ID X55053.1

    Use the Bio.SeqIO.parse(handle, format) function if you want
    to read multiple records from the handle.
    """
    iterator = parse(handle, format, alphabet)
    try:
        first = iterator.next()
    except StopIteration:
        first = None
    if first is None:
        raise ValueError("No records found in handle")
    try:
        second = iterator.next()
    except StopIteration:
        second = None
    if second is not None:
        raise ValueError("More than one record found in handle")
    return first

def to_dict(sequences, key_function=None):
    """Turns a sequence iterator or list into a dictionary.

     - sequences  - An iterator that returns SeqRecord objects,
                    or simply a list of SeqRecord objects.
     - key_function - Optional callback function which when given a
                    SeqRecord should return a unique key for the dictionary.

    e.g. key_function = lambda rec : rec.name
    or,  key_function = lambda rec : rec.description.split()[0]

    If key_function is ommitted then record.id is used, on the assumption
    that the records objects returned are SeqRecords with a unique id.
    
    If there are duplicate keys, an error is raised.

    Example usage, defaulting to using the record.id as key:

    >>> from Bio import SeqIO
    >>> filename = "GenBank/cor6_6.gb"
    >>> format = "genbank"
    >>> id_dict = SeqIO.to_dict(SeqIO.parse(filename, format))
    >>> print sorted(id_dict)
    ['AF297471.1', 'AJ237582.1', 'L31939.1', 'M81224.1', 'X55053.1', 'X62281.1']
    >>> print id_dict["L31939.1"].description
    Brassica rapa (clone bif72) kin mRNA, complete cds.

    A more complex example, using the key_function argument in order to
    use a sequence checksum as the dictionary key:

    >>> from Bio import SeqIO
    >>> from Bio.SeqUtils.CheckSum import seguid
    >>> filename = "GenBank/cor6_6.gb"
    >>> format = "genbank"
    >>> seguid_dict = SeqIO.to_dict(SeqIO.parse(filename, format),
    ...               key_function = lambda rec : seguid(rec.seq))
    >>> for key, record in sorted(seguid_dict.iteritems()):
    ...     print key, record.id
    /wQvmrl87QWcm9llO4/efg23Vgg AJ237582.1
    BUg6YxXSKWEcFFH0L08JzaLGhQs L31939.1
    SabZaA4V2eLE9/2Fm5FnyYy07J4 X55053.1
    TtWsXo45S3ZclIBy4X/WJc39+CY M81224.1
    l7gjJFE6W/S1jJn5+1ASrUKW/FA X62281.1
    uVEYeAQSV5EDQOnFoeMmVea+Oow AF297471.1

    This approach is not suitable for very large sets of sequences, as all
    the SeqRecord objects are held in memory. Instead, consider using the
    Bio.SeqIO.index() function (if it supports your particular file format).
    """    
    if key_function is None:
        key_function = lambda rec : rec.id

    d = dict()
    for record in sequences:
        key = key_function(record)
        if key in d:
            raise ValueError("Duplicate key '%s'" % key)
        d[key] = record
    return d

def index(filename, format, alphabet=None, key_function=None):
    """Indexes a sequence file and returns a dictionary like object.

     - filename - string giving name of file to be indexed
     - format   - lower case string describing the file format
     - alphabet - optional Alphabet object, useful when the sequence type
                  cannot be automatically inferred from the file itself
                  (e.g. format="fasta" or "tab")
     - key_function - Optional callback function which when given a
                  SeqRecord identifier string should return a unique
                  key for the dictionary.
    
    This indexing function will return a dictionary like object, giving the
    SeqRecord objects as values:

    >>> from Bio import SeqIO
    >>> records = SeqIO.index("Quality/example.fastq", "fastq")
    >>> len(records)
    3
    >>> sorted(records)
    ['EAS54_6_R1_2_1_413_324', 'EAS54_6_R1_2_1_443_348', 'EAS54_6_R1_2_1_540_792']
    >>> print records["EAS54_6_R1_2_1_540_792"].format("fasta")
    >EAS54_6_R1_2_1_540_792
    TTGGCAGGCCAAGGCCGATGGATCA
    <BLANKLINE>
    >>> "EAS54_6_R1_2_1_540_792" in records
    True
    >>> print records.get("Missing", None)
    None

    Note that this psuedo dictionary will not support all the methods of a
    true Python dictionary, for example values() is not defined since this
    would require loading all of the records into memory at once.

    When you call the index function, it will scan through the file, noting
    the location of each record. When you access a particular record via the
    dictionary methods, the code will jump to the appropriate part of the
    file and then parse that section into a SeqRecord.

    Note that not all the input formats supported by Bio.SeqIO can be used
    with this index function. It is designed to work only with sequential
    file formats (e.g. "fasta", "gb", "fastq") and is not suitable for any
    interlaced file format (e.g. alignment formats such as "clustal").

    For small files, it may be more efficient to use an in memory Python
    dictionary, e.g.

    >>> from Bio import SeqIO
    >>> records = SeqIO.to_dict(SeqIO.parse(open("Quality/example.fastq"), "fastq"))
    >>> len(records)
    3
    >>> sorted(records)
    ['EAS54_6_R1_2_1_413_324', 'EAS54_6_R1_2_1_443_348', 'EAS54_6_R1_2_1_540_792']
    >>> print records["EAS54_6_R1_2_1_540_792"].format("fasta")
    >EAS54_6_R1_2_1_540_792
    TTGGCAGGCCAAGGCCGATGGATCA
    <BLANKLINE>

    As with the to_dict() function, by default the id string of each record
    is used as the key. You can specify a callback function to transform
    this (the record identifier string) into your prefered key. For example:

    >>> from Bio import SeqIO
    >>> def make_tuple(identifier):
    ...     parts = identifier.split("_")
    ...     return int(parts[-2]), int(parts[-1])
    >>> records = SeqIO.index("Quality/example.fastq", "fastq",
    ...                       key_function=make_tuple)
    >>> len(records)
    3
    >>> sorted(records)
    [(413, 324), (443, 348), (540, 792)]
    >>> print records[(540, 792)].format("fasta")
    >EAS54_6_R1_2_1_540_792
    TTGGCAGGCCAAGGCCGATGGATCA
    <BLANKLINE>
    >>> (540, 792) in records
    True
    >>> "EAS54_6_R1_2_1_540_792" in records
    False
    >>> print records.get("Missing", None)
    None

    Another common use case would be indexing an NCBI style FASTA file,
    where you might want to extract the GI number from the FASTA identifer
    to use as the dictionary key.

    Notice that unlike the to_dict() function, here the key_function does
    not get given the full SeqRecord to use to generate the key. Doing so
    would impose a severe performance penalty as it would require the file
    to be completely parsed while building the index. Right now this is
    usually avoided.
    """
    #Try and give helpful error messages:
    if not isinstance(filename, basestring):
        raise TypeError("Need a filename (not a handle)")
    if not isinstance(format, basestring):
        raise TypeError("Need a string for the file format (lower case)")
    if not format:
        raise ValueError("Format required (lower case string)")
    if format != format.lower():
        raise ValueError("Format string '%s' should be lower case" % format)
    if alphabet is not None and not (isinstance(alphabet, Alphabet) or \
                                     isinstance(alphabet, AlphabetEncoder)):
        raise ValueError("Invalid alphabet, %s" % repr(alphabet))

    #Map the file format to a sequence iterator:    
    import _index #Lazy import
    try:
        indexer = _index._FormatToIndexedDict[format]
    except KeyError:
        raise ValueError("Unsupported format '%s'" % format)
    return indexer(filename, format, alphabet, key_function)

def to_alignment(sequences, alphabet=None, strict=True):
    """Returns a multiple sequence alignment (DEPRECATED).

     - sequences -An iterator that returns SeqRecord objects,
                  or simply a list of SeqRecord objects.  All
                  the record sequences must be the same length.
     - alphabet - Optional alphabet.  Stongly recommended.
     - strict   - Dummy argument, used to enable strict error
                  checking of sequence lengths and alphabets.
                  This is now always done.

    Using this function is now discouraged. You are now encouraged to use
    Bio.AlignIO instead, e.g.

    >>> from Bio import AlignIO
    >>> filename = "Clustalw/protein.aln"
    >>> alignment = AlignIO.read(filename, "clustal")
    """
    import warnings
    import Bio
    warnings.warn("The Bio.SeqIO.to_alignment(...) function is deprecated. "
                  "Please use the Bio.Align.MultipleSeqAlignment(...) object "
                  "directly instead.", Bio.BiopythonDeprecationWarning)
    return MultipleSeqAlignment(sequences, alphabet)

def convert(in_file, in_format, out_file, out_format, alphabet=None):
    """Convert between two sequence file formats, return number of records.

     - in_file - an input handle or filename
     - in_format - input file format, lower case string
     - out_file - an output handle or filename
     - out_format - output file format, lower case string
     - alphabet - optional alphabet to assume

    NOTE - If you provide an output filename, it will be opened which will
    overwrite any existing file without warning. This may happen if even
    the conversion is aborted (e.g. an invalid out_format name is given).

    For example, going from a filename to a handle:

    >>> from Bio import SeqIO
    >>> from StringIO import StringIO
    >>> handle = StringIO("")
    >>> SeqIO.convert("Quality/example.fastq", "fastq", handle, "fasta")
    3
    >>> print handle.getvalue()
    >EAS54_6_R1_2_1_413_324
    CCCTTCTTGTCTTCAGCGTTTCTCC
    >EAS54_6_R1_2_1_540_792
    TTGGCAGGCCAAGGCCGATGGATCA
    >EAS54_6_R1_2_1_443_348
    GTTGCTTCTGGCGTGGGTGGGGGGG
    <BLANKLINE>
    """
    if isinstance(in_file, basestring):
        #Hack for SFF, will need to make this more general in future
        if in_format in _BinaryFormats :
            in_handle = open(in_file, "rb")
        else :
            in_handle = open(in_file, "rU")
        in_close = True
    else:
        in_handle = in_file
        in_close = False
    #Don't open the output file until we've checked the input is OK?
    if isinstance(out_file, basestring):
        if out_format in ["sff", "sff_trim"] :
            out_handle = open(out_file, "wb")
        else :
            out_handle = open(out_file, "w")
        out_close = True
    else:
        out_handle = out_file
        out_close = False
    #This will check the arguments and issue error messages,
    #after we have opened the file which is a shame.
    from _convert import _handle_convert #Lazy import
    count = _handle_convert(in_handle, in_format,
                            out_handle, out_format,
                            alphabet)
    #Must now close any handles we opened
    if in_close:
        in_handle.close()
    if out_close:
        out_handle.close()
    return count
           
def _test():
    """Run the Bio.SeqIO module's doctests.

    This will try and locate the unit tests directory, and run the doctests
    from there in order that the relative paths used in the examples work.
    """
    import doctest
    import os
    if os.path.isdir(os.path.join("..", "..", "Tests")):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("..", "..", "Tests"))
        doctest.testmod()
        os.chdir(cur_dir)
        del cur_dir
        print "Done"
    elif os.path.isdir(os.path.join("Tests", "Fasta")):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("Tests"))
        doctest.testmod()
        os.chdir(cur_dir)
        del cur_dir
        print "Done"

if __name__ == "__main__":
    #Run the doctests
    _test()