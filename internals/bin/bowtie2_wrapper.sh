#!/bin/bash
# Hack to suppress skipped mate pair warnings from bowtie2-align

# (This works, but puts stderr into stdout.
# Tried redirecting stderr using process substitution, but that results in 
# bowtie2-align getting called with a final arg 2/dev/fd/63 and crashing)
bowtie2 "$@" 2>&1 | sed '/^Warning: skipping mate #/d'
exit ${PIPESTATUS[0]}


