#!/bin/bash

# see discussion at https://ma.juii.net/blog/scale-page-content-of-pdf-files
# and https://stackoverflow.com/questions/18343813/scale-pdf-to-add-border-for-printing-full-size-pages

xx="842" # width of A4 landscape in points
yy="595" # height of A4 landscape
char_width=3
margin=75
max_insert="${1}"
#bc not available on some platforms
#page_width=$(bc <<< "scale=4; ${char_width}*${max_insert}+2*${margin}")
page_width=$(python -c "print(${char_width}*${max_insert}+2*${margin})")
#scale_factor=$(bc <<< "scale=4; ${xx}/${page_width}")
scale_factor=$(python -c "print(\"{:.5f}\".format(${xx}/float(${page_width})))")

gs -o "${3}" \
   -q \
   -sDEVICE=pdfwrite \
   -dDEVICEWIDTHPOINTS=${xx} \
   -dDEVICEHEIGHTPOINTS=${yy} \
   -dFIXEDMEDIA \
   -dPDFFitPage \
   -dNOPAUSE \
   -dBATCH \
   -dSAFER \
   -dEmbedAllFonts=true \
   -dCompatibilityLevel=1.4 \
   -c "<</BeginPage{${scale_factor} ${scale_factor} scale}>> setpagedevice" \
   -f "${2}"


#   -sPAPERSIZE=a4 \
