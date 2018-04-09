#!/bin/bash
# Attempt to render markdown readme to html

BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. && pwd )"

VERSION="$(<${BASE_DIR}/release/version.txt)"

reformat () {
    sed -i "s/README.md - Grip/ShapeMapper ${VERSION} README/g" "$1"
    sed -i 's/border:\([0-9]*\)px solid #ddd/border:\1px solid #FF0000/g' "$1"
    sed -i 's/#FF0000}\.markdown-body/#ddd}\.markdown-body/g' "$1"
    sed -i 's/FF0000/FFF/g' "$1"
    # want to find the <div id="readme" . . . > opening tag, then remove
    # the first <h3> . . . </h3> element after that
    div_line=$(grep -n '<div id="readme"' "$1" | cut -d ':' -f 1)
    # just assume <h3> and </h3> are on separate lines and grep for their locations too
    start=$(( $(tail -n +${div_line} "$1" | grep -n -m1 '<h3>' | cut -d ':' -f 1) + ${div_line} - 1 ))
    end=$(( $(tail -n +${div_line} "$1" | grep -n -m1 '</h3>' | cut -d ':' -f 1) + ${div_line} - 1))
    sed -i "${start},${end}d" "$1"
}

grip "${BASE_DIR}/README.md" --export "${BASE_DIR}/README.html"
reformat "${BASE_DIR}/README.html"