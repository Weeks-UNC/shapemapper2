## ShapeMapper output file format descriptions

Primary outputs
---------------

Executable: `bin/make_reactivity_profiles.py`

Filename: `<name_<RNA>_profile.txt`
Format:
    Tab-delimited text columns. First line is column names.
Columns:
    "Nucleotide"                Nucleotide number (1-based)
    "Sequence"                  Nucleotide (AUGCaugc)
    "<Sample>_mutations"        Mutation counts
    "<Sample>_read_depth"       Read depth
    "<Sample>_effective_depth"  Description below
    "<Sample>_rate"             Effective mutation rate calculated as
                                (count/(depth-occluded_depth))
    "Reactivity_profile"        Calculated reactivity profile
                                ((Modified_rate-Untreated_rate)/Denatured_rate)
    "Std_err"                   Standard error
    "HQ_profile"                Reactivity profile with high-background and
                                low-depth positions excluded (set to nan)
    "HQ_stderr"                 Standard error with high-background and
                                low-depth positions excluded
    "Norm_profile"              Reactivity profile after normalization
                                (scaling)
    "Norm_stderr"               Standard error after normalization

    Effective depth explanation:
    Simple sequencing read depths are computed as the number of reads crossing
    a given position in the target sequence. For computing mutation rates,
    however, these depths are somewhat inadequate. The mutation rate at a 
    given position should be computed as the number of observed mutations
    divided by the total number of observations (i.e. the number of
    "no mutation" observations plus the number of "mutation" observations). 
    Mutations are excluded from being counted if they occur within 
    `--random-primer-len` of the 3-prime end of a read or if they contain or
    are neighbored by basecalls not meeting the quality threshold set by
    `--min-qual-to-count`. These same criteria are also used to compute the
    effective read depth. 

    Additional note:
    A given position within the region covered by a multinucleotide mutation is
    not actually an observation of "no mutation" - we are effectively blind
    over these regions in individual reads. Because of this, the region over
    each mutation is excluded from contributing to the effective read depth, 
    with the exception of the inferred adduct site. In practice, this
    correction is usually a small change to a large denominator, and has a 
    negligable effect on reactivity profile accuracy, at least with mutation 
    rates in the range produced by SHAPE adducts. 


Filename: `<name>_<RNA>.shape` 
Description:
    Reactivity profile in format expected by RNAstructure software. 
    Tab-delimited text. Two columns. First column is 1-based nucleotide
    position. Second is reactivity, with excluded positions set to -999.

Filename: `<name>_<RNA>.map`
Description:
    Same as .shape file, but with an additional two columns. Third column is
    stderr, fourth is nucleotide sequence.

Filename: `<name>_<RNA>_profiles.pdf`
Description:
    Figures showing read depths, mutation rates, and reactivity profile.
    
Filename: `<name>_<RNA>_histograms.pdf`
Description:
    Figures showing read depth, mutation rate and reactivity histograms.


Optional intermediate outputs
-----------------------------

# Aligned reads

Commandline option: `--output-aligned`
Executable: `bowtie2` or `STAR`
Filename: `*_aligned.sam`
Format: SAM

# Parsed mutations

Commandline option: `--output-parsed`
Executable: `bin/shapemapper_mutation_parser`
Filename: `<name>_<sample>_<RNA>_parsed.mut`
Format:
    Text, one line per mapped read (excluding reads with aligner-reported
    mapping qualities below 10). Each line is space-delimited.
    Fields:
    1) read name
    2) 0-based leftmost mapping position (inclusive)
    3) 0-based rightmost mapping position (inclusive)
    4) Target sequence (*not* read sequence) over mapped region
    5) Basecall quality scores for each target position in mapped region
       in ASCII encoding (phred score = ASCII value - 33)
  ...) Remaining fields (if any) are in groups of four. Each group represents
       a mutation from the target sequence present in the read sequence.
        1. 0-based nearest unchanged target sequence position on the left
        2. 0-based nearest unchanged target sequence position on the right
        3. Double-quoted read sequence that replaces the target sequence
           between the two positions.
        4. Double-quoted basecall quality scores for each read position
           in the previous field.

# Classified mutations

Commandline option: `--output-classified`
Executable: `bin/shapemapper_mutation_counter`
Filename: `<name>_<sample>_<RNA>_classified_mutations.txt`
Description:
    Text, one line per mapped sequence read, excluding reads with
    aligner-reported mapping qualities below 10. Each line is a list of
    mutations (if any) from the target sequence in a single read, not including
    any soft-clipped regions or any mutations overlapping 
    `--random-primer-len`+1 nucleotides from the 3-prime read end.
    
    Fields:
    1) read name
    2) 0-based leftmost mapping position (inclusive)
    3) 0-based rightmost mapping position (inclusive)
    4) Effective read depth over mapped region
    5) Effective mutation count over mapped region
  ...) Remaining fields (if any) are in groups of five. Each group represents
       a mutation from the target sequence present in the read sequence.
        1. 0-based nearest unchanged target sequence position on the left
        2. 0-based nearest unchanged target sequence position on the right
        3. Double-quoted read sequence that replaces the target sequence
           between the two positions.
        4. Double-quoted basecall quality scores for each read position
           in the previous field.
        5. Double-quoted mutation classification
    
Mutation classifications:
    "A-","T-","G-","C-" (single-nucleotide deletions)
    "-A","-T","-G","-C","-N" (single-nucleotide insertions)
    "AT", "AG", "AC",
    "TA", "TG", "TC",
    "GA", "GT", "GC",
    "CA", "CT", "CG" (single-nucleotide mismatches)
    "multinuc_deletion"
    "multinuc_insertion"
    "multinuc_mismatch"
    "complex_deletion"
    "complex_insertion"
    "N_match" (not a real mutation, just an ambiguous basecall)
    
    All classifications may also have "_ambig" appended if a mutation involves 
    any ambiguously aligned nucleotides.


# Mutation counts

Commandline option: `--output-counted`
Executable: `bin/shapemapper_mutation_counter`
Filename: `<name>_<sample>_<RNA>_mutation_counts.txt`
Format: 
    Tab-delimited text columns. First line is column headers. One column for
    each mutation classification listed in the previous section, with the
    exception of "N_match". These columns contain mutation counts listed
    5-prime to 3-prime. Multi-nucleotide mutations are counted at their 3-prime
    end. Columns "read_depth" and "effective_depth" contain sequencing depths
    calculated as previously described above.




Last updated Sep. 2017, Steven Busan.