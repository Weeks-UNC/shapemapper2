"""
Replace input stream with bad chars

"""

# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2017 Steven Busan.                     #
# --------------------------------------------------------------------- #

import sys
import os

bad_char = '$'

try:
    for line in sys.stdin:
        for c in line:
            if c!='\n':
                c = bad_char
            sys.stdout.write(c)
except BrokenPipeError:
    pass
    # hack attempting to suppress error messages
    # see http://stackoverflow.com/questions/16314321/suppressing-printout-of-exception-ignored-message-in-python-3
    #sys.stdout = os.fdopen(1)
