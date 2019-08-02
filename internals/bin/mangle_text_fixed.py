#!/usr/bin/env python3
"""
Replace input stream with bad chars. 
Try to mostly preserve lines and whitespace, but mess with
non-whitespace characters and the lengths of each field.

"""

# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2018 Steven Busan.                     #
# --------------------------------------------------------------------- #

import sys
import os

bad_char = '$'

try:
    for line in sys.stdin:
        i = 0
        for c in line:
            if c not in ' \n\t':
                c = bad_char
                i += 1
            if i % 2 == 0:
                sys.stdout.write(c)
except BrokenPipeError:
    pass
    # hack attempting to suppress error messages
    # see http://stackoverflow.com/questions/16314321/suppressing-printout-of-exception-ignored-message-in-python-3
    #sys.stdout = os.fdopen(1)
