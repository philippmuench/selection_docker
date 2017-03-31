#/bin/sh

# runs codeml on example files and displays the dN/dS ratio
echo "running codeml"
codeml codeml.ctl > /dev/null 2> /dev/null
grep "omega (dN/dS) = " test.mlc

# pairwise codeml
echo "running codeml in a pairwise fahsion"
codeml codeml_pairwise.ctl > /dev/null 2> /dev/null
# get mean of all pairwise comparisons
grep "dN/dS=" test_pair.mlc | awk '{print $8}' | awk -F : '{sum+=$1} END {print "average dN/dS=",sum/NR}'

# fastcodeml
echo "running fastcodeml"
fast -m 22 -bf -hy 1 tree.phy codon.phy > fastcodeml_output.txt

