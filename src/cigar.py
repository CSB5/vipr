"""CIGAR operations.

There are two ways to represent CIGAR strings: Either as regular strings,
such as "17M1D5M4S", or in an encoded, numeric form as a list of pairs:
[ (17, 0), (2, 1), (5, 0), (4, 0) ]

The mapping of CIGAR operator to numbers is: MIDNSHP => 0123456

The encoded form is used by pysam.


Andreas Wilm: Downloaded from
https://sqt.googlecode.com/git/sqt/cigar.py on 2012-10-11. Part of
sqt: Command-line tools for second generation sequencing data See
https://code.google.com/p/sqt/ Slightly modified original.
"""
from itertools import repeat, chain, groupby
__author__ = 'Marcel Martin'


# cigar table from http://samtools.sourceforge.net/SAM1.pdf
# M 0 alignment match (can be a sequence match or mismatch)
# I 1 insertion to the reference
# D 2 deletion from the reference
# N 3 skipped region from the reference
# S 4 soft clipping (clipped sequences present in SEQ)
# H 5 hard clipping (clipped sequences NOT present in SEQ)
# P 6 padding (silent deletion from padded reference)
# = 7 sequence match
# X 8 sequence mismatch


# constants
M = 0
I = 1
D = 2
N = 3
S = 4
H = 5
P = 6

# use this as a sequence to map an encoded operation to the appropriate
# character
DECODE = 'MIDNSHP'

# this dictionary maps operations to their integer encodings
_ENCODE = dict( (c,i) for (i, c) in enumerate(DECODE) )

def parse(cigar):
    """
	Parse CIGAR string and return a list of pysam-encoded tuples.
    
	MIDNSHP => 0123456
    
	>>> parse("3S17M8D4M9I3H")
	[(4, 3), (0, 17), (2, 8), (0, 4), (1, 9), (5, 3)]
	"""
    result = []
    n = ''
    for c in cigar:
        if c.isdigit():
            n += c
        elif c in _ENCODE:
            if n == '':
                raise ValueError("end of CIGAR string reached, but an operator was expected")
            result.append( (_ENCODE[c], int(n)) )
            n = ''
    return result


def concat(left, right):
    """
    Concatenate two CIGARs that are given as pysam tuples.

    MIDNSHP => 0123456

    >>> concat(parse("1M"), parse("3M"))
    [(0, 4)]
    """
    if left and right:
        left_last = left[-1]
        right_first = right[0]
        # same operation?
        if left_last[0] == right_first[0]:
            right[0] = ( right[0][0], right[0][1] + left[-1][1] )
            left = left[:-1]
    return left + right


def aligned_length(cigar):
    """
    Return aligned length of CIGAR string given as a pysam-encoded
    list of tuples.
    """
    length = 0
    for op, l in cigar:
        if op in [0, 1]: # MI
            length += l
        elif op == 2: # D
            pass
        else:
            raise ValueError("CIGAR operation %s not supported" % cigar)
    return length


def ops(encoded_cigar):
    """
    Yield all operations (as numbers, not characters) one by one.

    >>> list(ops(parse("3S2I3M")))
    [4, 4, 4, 1, 1, 0, 0, 0]
    """
    return chain.from_iterable(repeat(op, l) for (op, l) in encoded_cigar)


def decoded_ops(encoded_cigar):
    """
    Yield all operations (as characters) one by one.

    >>> ''.join(decoded_ops(parse("3S2I3M")))
    'SSSIIMMM'
    """
    return chain.from_iterable(repeat(DECODE[op], l) for (op, l) in encoded_cigar)


def cigar_from_decoded_ops(decoded_ops):
    """Convert decoded operations into cigar string. For some reason
    missing in cigar.py. Added to original code (AW)

    >>> cigar_from_decoded_ops(decoded_ops(parse("3S2I3M")))
    '3S2I3M'
    """

    # alternative:
    #from collections import deque
    #dq = deque(new_cigar_decoded)
    #while len(dq)>1:
    #    conc = dq.popleft() + dq.popleft()
    #    dq.append(conc)
    #new_cigar = dq[0]

    cigar = ""
    for (k, g) in groupby(decoded_ops, lambda x: x):
        cigar = "%s%d%c"  % (cigar, len(list(g)), k)
    return cigar


if __name__ == "__main__":
    import sys
    sys.stderr.write("Go away. Nothing to see here")

