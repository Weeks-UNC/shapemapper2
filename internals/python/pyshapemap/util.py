# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2018 Steven Busan.                     #
# --------------------------------------------------------------------- #

import os
import random
import string
import datetime
import shutil
import fcntl
from copy import deepcopy

def non_block_read(output):
    """
    Non-blocking read. Used to periodically check Pipe Viewer
    output without hanging.

    from http://stackoverflow.com/questions/375427/non-blocking-read-on-a-subprocess-pipe-in-python

    """
    fd = output.fileno()
    fl = fcntl.fcntl(fd, fcntl.F_GETFL)
    fcntl.fcntl(fd, fcntl.F_SETFL, fl | os.O_NONBLOCK)
    try:
        return output.read()
    except:
        return ""


def get_extension(filepath):
    assert isinstance(filepath, str)
    extension = ".".join(os.path.basename(filepath).split(".")[1:])
    return extension


def rand_id():
    return "".join([random.choice(string.ascii_lowercase) for x in range(8)])


def timestamp():
    return '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())


def version():
    this_dir = os.path.dirname(os.path.realpath(__file__))
    release_dir = os.path.join(this_dir, "../../release")
    f = open(os.path.join(release_dir, "version.txt"), "rU")
    return f.readline().strip()


def require_explicit_kwargs(locs):
    """
    Given the dict returned by locals() inside a function,
    will check that no objects are None
    """
    for x in locs:
        if locs[x] is None:
            raise RuntimeError('Missing required kwarg "{}"'.format(x))


def format_message(name,
                   level,
                   message):
    """
    Indent a message to a particular level and wrap inside
    terminal-width separators.

    """
    # TODO: wrap long lines to fit indent level
    if message is None or len(message)<1:
        return ''
    msg = ''
    lev = ''.join([' ' for x in range(level)])
    s = "{}/".format(lev)
    term_width, term_height = shutil.get_terminal_size()
    pad = ''.join(['`' for x in range(term_width - len(s))])+'\n'
    msg += s+pad
    s = "{} output message:".format(name)
    underline = ''.join(['-' for x in range(len(s))])
    s = "{}|{} \n".format(lev, s)
    s += "{}|{} \n".format(lev, underline)
    msg += s
    for line in ['']+message.splitlines():
        msg += "{}| {}\n".format(lev, line)
    s = "{}\\".format(lev)
    term_width, term_height = shutil.get_terminal_size()
    pad = ''.join(['_' for x in range(term_width - len(s))])+'\n'
    msg += s+pad
    return msg

class Logger:
    """
    Duplicate stream writes to a text file.

    """

    def __init__(self,
                 filepath=None,
                 stream=None):
        self.file=None
        if filepath is not None:
            folder, filename = os.path.split(filepath)
            if len(folder) > 0:
                os.makedirs(folder, exist_ok=True)
            self.file = open(filepath, "a")
        self.stream = stream

    def write(self, s):
        if self.stream is not None:
            self.stream.write(s)
        if self.file is not None:
            # bit of a hack to avoid writing progress bar
            # updates to the logfile
            # - skip all text between a carriage return and a newline or carriage return
            while(True):
                i = s.find("\r")
                if i==-1:
                    break
                j1 = s.find("\r", i)
                j2 = s.find("\n", i)
                j = min(j1, j2)
                if j==-1:
                    j = len(s)
                s = s[:i]+s[j:]
            self.file.write(s)

        self.flush()

    def flush(self):
        if self.stream is not None:
            self.stream.flush()
        if self.file is not None:
            self.file.flush()


def sanitize(s,
             replace_whitespace=True,
             allow_slash=False,
             check_directory_traversal=False,
             check_root_path=False):
    """
    Replace disallowed chars in a string.

    Optionally check for filesystem traversals.
    """
    assert isinstance(s, str)

    allowed_chars = string.ascii_letters + string.digits + ".-_+:"

    if allow_slash:
        allowed_chars += '/'
    if not replace_whitespace:
        allowed_chars += string.whitespace
    new_str = ""
    for c in s:
        if c not in allowed_chars:
            new_str += '_'
        else:
            new_str += c
            # raise RuntimeError(
            #    'Error: "' + s + '" contains disallowed characters.\nAllowed characters are: "' + allowed_chars + '"')
    if check_directory_traversal:
        if "../" in s:
            raise RuntimeError(
                'Error: "' + s + '" contains a directory traversal sequence.\nFile paths are not allowed to point outside of the current directory.')
    if check_root_path:
        if s.startswith('/'):
            raise RuntimeError(
                'Error: "' + s + '" references the root filesystem path.\nFile paths must be relative to the current directory.')
    return new_str


def read_fasta_names_lengths(fastas):
    # Get fasta sequence names from fasta file(s) and return list.
    # Check for duplicated names
    assert isinstance(fastas, list)
    def check_no_dups(l):
        if len(list(set(l))) != len(l):
            # TODO: list fasta filenames in error message
            raise RuntimeError("Error: fasta file(s) contain duplicated sequence names")
        sanitized_names = [sanitize(s) for s in l]
        if len(list(set(sanitized_names))) != len(sanitized_names):
            msg = "Error: after replacing disallowed characters in sequence names, "
            msg += " some sequences have the same name."
            raise RuntimeError(msg)

    names = []
    lengths = []
    lowercase_error = "One or more target sequences have all lowercase characters. "
    lowercase_error += "Please provide at least some unmasked sequence (uppercase characters) for each sequence."

    for filename in fastas:
        f = open(filename, "rU")
        is_all_lower = False
        for line in f:
            if line[0] == '>':
                names.append(line[1:].rstrip())
                lengths.append(0)
                if is_all_lower:
                    raise RuntimeError(lowercase_error)
                is_all_lower = True
            else:
                s = ''.join(line.strip().split())
                if len(s) > 0:
                    if len(lengths)==0:
                        raise RuntimeError("Error: fasta file missing sequence name (should look like '>name'")
                    lengths[-1] += len(s)
                    islow = s.islower()
                    if not islow:
                        is_all_lower = False
        if is_all_lower:
            raise RuntimeError(lowercase_error)
    check_no_dups(names)
    return names, lengths




