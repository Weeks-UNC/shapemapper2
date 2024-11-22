#!/usr/bin/env python3
import sys, os
from argparse import ArgumentParser

def iterate_fastq(file_obj):
    lines = []
    for line in file_obj:
        # Suppress jdb socket message
        if line.startswith("Listening"):
            continue
        lines.append(line)
        if len(lines) == 4:
            #l = [lines[0].split()[0][1:], lines[1].strip(), "+", lines[3].strip()]
            l = [lines[0].split()[0], lines[1].strip(), "+", lines[3].strip()]
            yield l
            lines = []
        

def unpaired_gen(inp_file_obj):
    iterate_fastq_gen = iterate_fastq(inp_file_obj)
    reads = []
    try:
        reads.append(next(iterate_fastq_gen))
    except StopIteration: 
        raise EOFError("Error! Empty merger output, check shapemapper log for other error messsages.")

    else: # No stop iteration raised, meaning file is not empty
        for r in iterate_fastq(inp_file_obj):
            reads.append(r)

            if len(reads) == 2:
                r1_id = reads[0][0]
                r2_id = reads[1][0]

                if r2_id == r1_id:    # If paired
                    reads = []        # Pass

                else:                 # If unpaired
                    for line in reads[0]:
                        yield line
                    reads = [reads[1]]
        
        if len(reads) == 1: # If there is a read leftover it must be unpaired
            for line in reads[0]:
                yield line


def paired_gen(inp_file_obj, which_read: str): # which_read should correspond to either R1 or R2

    if which_read not in ["R1", "R2"]:
        raise ValueError('Error! which_read must be either "R1" or "R2"')

    iterate_fastq_gen = iterate_fastq(inp_file_obj)
    reads = []
    try:
        reads.append(next(iterate_fastq_gen))
    except StopIteration: 
        raise EOFError("Error! Empty merger output, check shapemapper log for other error messsages.")

    else: # No stop iteration raised, meaning file is not empty
        for r in iterate_fastq(inp_file_obj):
            reads.append(r)

            if len(reads) == 2:
                r1_id = reads[0][0]
                r2_id = reads[1][0]

                if r2_id == r1_id:    # If paired yield R1 or R2 depending on input
                    if which_read == "R1":
                        for line in reads[0]:
                            yield line
                    elif which_read == "R2":
                        for line in reads[1]:
                            yield line
                    reads = []        

                else:                 # If unpaired, pass
                    reads = [reads[1]] 
        
        # If there is a read leftover it must be unpaired, so we don't need to deal with cases
        # where there is one read left

# Simple wrapper that buffers individual lines into reads then formats accordingly
def format_generator_lines(generator):
    lines = []
    for l in generator:
        lines.append(l)

        if len(lines) == 4:
            #yield "@" + lines[0] + "\n"
            #for i in range(1, 4):
            for i in range(0, 4):
                yield lines[i] + "\n"
            lines = []

# Pick generator based on input
def strategy(inp_obj, R1 = False, R2 = False, unpaired = False):
    if R1:
        return format_generator_lines(paired_gen(inp_obj, "R1"))
    elif R2:
        return format_generator_lines(paired_gen(inp_obj, "R2"))
    elif unpaired:
        return format_generator_lines(unpaired_gen(inp_obj))


if __name__ == "__main__":
    ap = ArgumentParser()
    ap.add_argument("--input", type=str)
    ap.add_argument("--output", type=str)
    ap.add_argument("--split_output", default = None, type=str) # Sometimes need to split output
                                                                # into two files.

    group = ap.add_mutually_exclusive_group()
    group.add_argument("--R1", default = False, action="store_true")
    group.add_argument("--R2", default = False, action="store_true")
    group.add_argument("--unpaired", default = False, action="store_true")

    args = ap.parse_args()

    
    if sum([args.R1, args.R2, args.unpaired]) == 0:
        raise ValueError("Must pass --R1, --R2, or --unpaired")
    
    if args.split_output:
        with open(args.input, "r") as inp_obj, open(args.output, "w") as out_obj, \
             open(args.split_output, "w") as spl_out_obj:
            gen = strategy(inp_obj, args.R1, args.R2, args.unpaired)
            for line in gen:
                out_obj.write(line)
                spl_out_obj.write(line)

    else:
        with open(args.input, "r") as inp_obj, open(args.output, "w") as out_obj:
            gen = strategy(inp_obj, args.R1, args.R2, args.unpaired)
            for line in gen:
                out_obj.write(line)
