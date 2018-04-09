ShapeMapper quality control checks
==================================

ShapeMapper performs several basic quality-control checks for each 
reactivity profile. These are necessarily heuristic, since different
downstream analyses require different levels of data quality, and 
since RNAs of differing flexibility will show different overall
signal levels above background. In general, more sequencing read 
depth is always helpful, as are higher modification/mutation rates 
(to a point).


Read depth check
----------------

We require that at least 80% of non-masked nucleotides meet a minimum
sequencing depth of 5000 in all provided samples. Note that with
uneven sequencing coverage, some regions of an RNA may have 
higher-quality reactivity data than other regions. For example, 
coverage is often lower near transcript ends.

If read depths are low, check alignment stats to see the amount of 
target sequence present in each sample. Better target enrichment or 
recovery and/or additional sequencing can often help.


Positive mutation rates above background check
----------------------------------------------

If an untreated control sample is provided (highly recommended), we
require that at least 50% of non-masked nucleotides with depths above
5000 have a higher mutation rate in the modified sample than in the 
untreated sample.


High background mutation rates check
------------------------------------

If an untreated control sample is provided, we require that no more
than 5% of non-masked nucleotide with depths above 5000 have an
untreated mutation rate above 0.05.

An unusual number of high-background nucleotides can result from the
presence of native modifications or sequence variants.


Number of highly reactive nucleotides check
-------------------------------------------

We require that at least 8% of non-masked nucleotides with depths 
above 5000 have a modified mutation rate above 0.006 after background
subtraction.

Possible causes for failure: 
    - DNA contamination. Unusually low background mutation rates can 
      be a secondary indication that this is the problem (since 
      reverse transcription under MaP conditions usually generates
      some errors).

    - Poor mixing of chemical reagents and RNA and/or poor reagent 
      diffusion (if modifying in cells), resulting in low 
      modification rates

    - Expired reagents, resulting in low modification rates

    - Poor reverse transcription conditions, resulting in low adduct 
      read-through

    - Extremely highly structured RNA. A molecule that genuinely
      contains no flexible nucleotides is indistinguishable from a
      highly flexible one that was mistakenly unmodified (they will
      both have low mutation rates above background). In this case, 
      additional control experiments or complementary techniques may 
      be needed.




Last updated Mar. 2017, Steven Busan.