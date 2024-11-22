#!/bin/bash
# Hack to suppress skipped mate pair warnings from bowtie2-align

# (This works, but puts stderr into stdout.
# Tried redirecting stderr using process substitution, but that results in 
# bowtie2-align getting called with a final arg 2/dev/fd/63 and crashing)

#28 April 2023 - appened sed expressions to remove odd bowtie2 seed mismatch glitch stemming from 
bowtie2 "$@" 2>&1 | sed '/^Warning: skipping mate #/d' | sed '/ because it was < 2 characters long$/d' | sed '/ because length (1) <= # seed mismatches (0)$/d'
exit ${PIPESTATUS[0]}


