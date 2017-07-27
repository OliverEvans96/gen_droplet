#!/bin/bash
awk '{if (NR==5) print > "'"$2"'"; if(NR>5) print $7,$2,$3,$4,$9 > "'"$2"'"}' $1 
