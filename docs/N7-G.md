<!---
NOTE:
If you're reading this, instead try opening README.html in a web browser 
or view this file from within the github repository website.

This is a github-flavored markdown file not meant to be easily readable.
-->

N7 mode
========

Overview
--------
The --N7 option was added in v2.3 to allow processing and analysis of N7-G
modifications induced under the experimental conditions specified for dms
treatment - [DMSmode](docs/dmsmode.md). Traditionally, measurement and analysis
of N7-G modifications on RNA has involved harsh biochemical processing separate
from conventional N1/3 modification analysis. We have developed msDMS-MaP to
allow simultaneous detection and analysis of both N1/3 and N7-G modifications.
We have shown that N7-G reactivity is informative about RNA tertiary and quaternary
structure. The N7-G reactivity analysis complements the traditional N1/3 reactivity
analysis which is conventionally interpreted in the context of secondary structure. 

Under the aformentioned experimental conditions, these N7-G modifications manifest
as G>A mutations. Shapemapper isolates and processes these modifications in a 
channel separate from N1/3 modifications. When the --N7 flag is used, N7-G related
information will be written to a profile.txtga file (and one or more .mutga files
if the --output-parsed-mutations flag is used). Additionally, the N7-G data will
be visualized alongside N1/3 data in the profile.pdf output file.


Experimental conditions
------------------
see [DMSmode](docs/dmsmode.md)


Normalization
-------------
Please see [place publication here] for a detailed description of N7-G normalization.

Due to the way these sites are normalized, the raw rate has as inverse correlation
with normalized reactivity. In other words, an N7-G position with a low raw rate will 
have a high normalized reactivity. We term sites with high normalized reactivity "protected".
Additionally, N7-G reactivity normalization incorporates a log2 transformation. Thus,
successive increments of 1 correspond to a "doubling" of the N7-G normalized reactivity.
For example a normalized reactivity of 2 is twice as protected as a normalized reactivity of 1
due to the preceding log2 transformation.

Based on prior experiments, we have set thresholds to determine how protected each N7-G
position is. Cutoffs have been set at 1.6 and 2.3 corresponding to bases which
are protected and highly protected respectively. In the profiles.pdf data visualization
unprotected bases (N reactivity < 1.6) are colored black, protected bases 
(1.6 <= N reactivity < 2.3) are colored pink, and highly protected bases (N reactivity >= 2.3) 
are colored purple.


Further Analysis
------------------
Additional analysis of N7-G data may be performed in [RingMapper](https://github.com/Weeks-UNC/RingMapper), [DanceMapper](https://github.com/MustoeLab/DanceMapper), and [ArcPlot](https://github.com/MustoeLab/StructureAnalysisTools).

Each of these packages has functionality specific to N7-G data processing.


Citation and reference
----------------------
Please cite Saleem and Miller et al, Journal To Be Determined, 202X, for publications using the --N7 option

Please cite Mitchell et al, Nucleic Acids Research, 2023, for publications using the --dms option

Bicine buffering conditions were first described in Mustoe et al, PNAS 2019

