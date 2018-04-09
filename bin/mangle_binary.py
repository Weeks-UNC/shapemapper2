"""
Insert random bytes into a binary stream

"""

# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2017 Steven Busan.                     #
# --------------------------------------------------------------------- #

import os, sys, random, string

p = 0.1
leading_bytes = 100 # initial bytes to ignore (avoid corrupting most headers)

possible_chars = string.printable

try:
    i = 0
    byte = sys.stdin.buffer.read(1)
    while byte:
        i += 1
        if ( i > leading_bytes and
             random.random() < p ):
            sys.stdout.buffer.write(random.choice(possible_chars).encode('utf-8'))
        sys.stdout.buffer.write(byte)
        byte = sys.stdin.buffer.read(1)
except BrokenPipeError:
    # hack attempting to suppress broken pipe warning messages
    # - only partially successful
    # see http://stackoverflow.com/questions/16314321/suppressing-printout-of-exception-ignored-message-in-python-3
    sys.stdout = os.fdopen(1)
    sys.stderr = os.fdopen(1)