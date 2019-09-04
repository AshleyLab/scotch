#!/bin/bash

fm="$1"
out="$2"

toSlide="15,16,17,18,19,20,21,22,23,24"
nToSlide="10"
useToFill="0"
delimiter="\t"

fillLine=$(yes "$useToFill" | head -n "$nToSlide" | tr '\n' "$delimiter" | sed "s|$delimiter\$||")

paste "$fm" <(cut -f"$toSlide" "$fm" | tail -n +2; echo "$fillLine") | gzip > "$out"
