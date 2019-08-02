#!/usr/bin/env python3
# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2018 Steven Busan.                     #
# --------------------------------------------------------------------- #


import sys, os


def check_fasta(fasta_path,
                corrected_fasta_path):
    """
    Check a fasta file for format requirements, raise errors if needed.
    Attempt to correct errors and write to a new file. Return True if no
    errors, False otherwise.

    """

    fa = open(fasta_path,"r")
    read = fa.read()
    fa.seek(0)
    lines = fa.readlines()
    fa.close()

    corrected_msgs = []
    uncorrected_msgs = []

    if "\r" in read:
        lines = [line.strip()+"\n" for line in lines]
        corrected_msgs.append("Line endings not in unix format.")

    if len(lines) < 2:
        uncorrected_msgs.append("Not enough lines (need at least a header and a sequence).")

    found_U = False
    found_whitespace = False
    found_nonalpha = False
    missing_header = False
    missing_seq = False
    extra_headers = False
    dup_headers = False

    fixed_lines = [str(line) for line in lines]
    break_flag = False
    headers = []
    header_count = 0
    seq_count = 0
    for i in range(len(fixed_lines)):
        if fixed_lines[i][0] == '>':
            if header_count != seq_count:
                extra_headers = True
                break
            headers.append(lines[i][1:].rstrip())
            header_count += 1
            continue
        elif len(lines[i].strip()) != 0:
            if header_count == seq_count+1:
                seq_count += 1
            elif header_count == 0:
                missing_header = True
        for j in range(len(fixed_lines[i])):
            c = fixed_lines[i][j]
            if c == "\n":
                continue
            elif c.upper() == "U":
                found_U = True
            elif c == ' ':
                found_whitespace = True
            elif not c.isalpha():
                if header_count == 0:
                    missing_header = True
                    break_flag = True
                    break
                else:
                    found_nonalpha = True
                    break_flag = True
                    break
        if found_whitespace or found_U:
            fixed_lines[i] = fixed_lines[i].replace('u', 't').replace('U','T').replace(' ','')
        if break_flag:
            break
    if header_count < seq_count:
        missing_header = True
    elif header_count > seq_count:
        missing_seq = True
    if len(list(set(headers))) != len(headers):
        dup_headers = True

    if found_U:
        msg = "One or more 'U' found in sequence (should be 'T')."
        corrected_msgs.append(msg)
    if found_whitespace:
        msg = "Space(s) found within sequence (STAR aligner interprets these as Ns)."
        corrected_msgs.append(msg)
    if found_nonalpha:
        msg = "One or more non-alphabetic characters found in sequence"
        msg += " or sequence name is misformatted (must start with '>')."
        uncorrected_msgs.append(msg)
    if missing_header:
        msg = "Missing header(s). Preceding each sequence, each sequence"
        msg += " name should be on its own line and start with the '>' character."
        uncorrected_msgs.append(msg)
    if missing_seq:
        msg = "Missing one or more sequences."
        uncorrected_msgs.append(msg)
    if extra_headers:
        msg = "Extra header(s), missing sequence(s), or headers in"
        msg += " wrong location. Preceding each sequence, each sequence"
        msg += " name should be on its own line and start with the '>' character."
        uncorrected_msgs.append(msg)
    if dup_headers:
        msg = "Duplicated sequence names."
        uncorrected_msgs.append(msg)
    msg = ''

    if len(uncorrected_msgs) == 0:
        fa = open(corrected_fasta_path,"w")
        fa.write("".join(fixed_lines))
        fa.close()
    if len(corrected_msgs) > 0 and len(uncorrected_msgs) == 0:
        msg += "ERROR: FASTA file {} contains errors:\n".format(fasta_path)
        for m in corrected_msgs:
            msg += " - {}\n".format(m)
        msg += "These errors were corrected and written to {}\n".format(corrected_fasta_path)
        msg += "Replace the original FASTA file with the new one\n"
        msg += "before rerunning shapemapper.\n"
    if len(uncorrected_msgs) > 0:
        msg += "ERROR: FASTA file {} contains errors:\n".format(fasta_path)
        for m in uncorrected_msgs:
            msg += " - {}\n".format(m)
        for m in corrected_msgs:
            msg += " - {}\n".format(m)
    if len(msg) > 0:
        sys.stderr.write(msg)
        return False
    else:
        return True


if __name__=="__main__":
    if len(sys.argv) != 3:
        print("Usage: python3.5 check_fasta_format.py <example.fa> <corrected.fa>")
        sys.exit(0)
    success = check_fasta(sys.argv[1],
                          sys.argv[2])
    if success:
        print("No formatting errors detected.")
    else:
        sys.exit(1)
