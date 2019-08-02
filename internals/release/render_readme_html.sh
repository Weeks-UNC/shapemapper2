#!/bin/bash
# Attempt to render markdown readme to html

BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd ../.. && pwd )"

VERSION="$(<${BASE_DIR}/internals/release/version.txt)"

replace_links () {
    # swap out relative links to .md files for .html
    # this may be a bit much, doesn't check for href tags or anything
    sed -i 's/\.md/\.html/g' "$1"
}

fix_anchors () {
    # remove user-content- prefixes in anchors. not sure why they're added - they break links
    sed -i 's/<a name="user-content-/<a name="/g' "$1"
}


test_replace () {
    echo '<>>>,:{}>>>heyHeyHey<a href="docs/analysis_steps.md#quality-control-checks">>' > tmp
    expected='<>>>,:{}>>>heyHeyHey<a href="docs/analysis_steps.html#quality-control-checks">>'
    replace_links tmp
    actual=$(cat tmp)
    if [ "$expected" != "$actual" ]; then
       echo "$expected"
       echo "$actual"
       echo "replace_links() FAILED"
       exit 1
    fi
    rm tmp
}

test_fix_anchors() {
    echo '<p><a name="user-content-warning"></a></p>' > tmp
    expected='<p><a name="warning"></a></p>'

    fix_anchors tmp
    actual=$(cat tmp)
    if [ "$expected" != "$actual" ]; then
        echo "$expected"
        echo "$actual"
        echo "fix_anchors() FAILED"
        rm tmp
        exit 1
    fi
    rm tmp
}


reformat_main () {
    replace_links "$1"
    fix_anchors "$1"

    # remove border elements
    sed -i "s/README.md - Grip/ShapeMapper ${VERSION} README/g" "$1"
    sed -i 's/border:\([0-9]*\)px solid #ddd/border:\1px solid #FF0000/g' "$1"
    sed -i 's/#FF0000}\.markdown-body/#ddd}\.markdown-body/g' "$1"
    sed -i 's/FF0000/FFF/g' "$1"

    # make the page a bit narrower
    sed -i "s/width:980px/width:720px/g" "$1"

    # want to find the <div id="readme" . . . > opening tag, then remove
    # the first <h3> . . . </h3> element after that
    div_line=$(grep -n '<div id="readme"' "$1" | cut -d ':' -f 1)
    # just assume <h3> and </h3> are on separate lines and grep for their locations too
    start=$(( $(tail -n +${div_line} "$1" | grep -n -m1 '<h3>' | cut -d ':' -f 1) + ${div_line} - 1 ))
    end=$(( $(tail -n +${div_line} "$1" | grep -n -m1 '</h3>' | cut -d ':' -f 1) + ${div_line} - 1))
    sed -i "${start},${end}d" "$1"
}

reformat () {
    replace_links "$1"
    fix_anchors "$1"

    # remove border elements
    sed -i 's/border:\([0-9]*\)px solid #ddd/border:\1px solid #FF0000/g' "$1"
    sed -i 's/#FF0000}\.markdown-body/#ddd}\.markdown-body/g' "$1"
    sed -i 's/FF0000/FFF/g' "$1"

    # make the page a bit narrower
    sed -i "s/width:980px/width:720px/g" "$1"

    # want to find the <div id="readme" . . . > opening tag, then remove
    # the first <h3> . . . </h3> element after that
    div_line=$(grep -n '<div id="readme"' "$1" | cut -d ':' -f 1)
    # just assume <h3> and </h3> are on separate lines and grep for their locations too
    start=$(( $(tail -n +${div_line} "$1" | grep -n -m1 '<h3>' | cut -d ':' -f 1) + ${div_line} - 1 ))
    end=$(( $(tail -n +${div_line} "$1" | grep -n -m1 '</h3>' | cut -d ':' -f 1) + ${div_line} - 1))
    sed -i "${start},${end}d" "$1"
}

grip "${BASE_DIR}/README.md" --export "${BASE_DIR}/README.html"
reformat_main "${BASE_DIR}/README.html"

# also render all .md files in docs/
for f in ${BASE_DIR}/docs/*.md; do
    name=$(basename $f .md)
    opath="${BASE_DIR}/docs/${name}.html"
    grip "$f" --export "${opath}"
    reformat "${opath}"
done
