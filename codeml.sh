#/bin/sh
# runs codeml on files in folder and displays the dN/dS ratio
echo "file;method;dnds;dN;dS"
# cleanup
rm -rf output
mkdir -p output
rm codon/*.reformatted.*
rm codon/*.reformatted
for codon_file in codon/*.msa.*; do
  basename=${codon_file##*/}
  # preprocess input files
  perl -p -e 's/\t/  /g' $codon_file > $codon_file.reformatted
  tree_file="tree/${basename%.*.*}.ffn.phy"
  output_file="output/${basename%.*.*}.txt"
  # add file location to codeml config file
  sed -i "1s|.*|seqfile = $codon_file.reformatted * sequence data file name|" codeml.ctl
  sed -i "2s|.*|treefile = $tree_file * tree structure file name|" codeml.ctl
  sed -i "3s|.*|outfile = $output_file * main result file name|" codeml.ctl
  codeml codeml.ctl > /dev/null 2> /dev/null
  dnds_codeml=$(grep "omega (dN/dS) = " $output_file | awk '{print $4}')
  # sum up over all branches in tree
  dN=$(cat $output_file |\
    sed -n -e '/dN & dS for each branch/,/tree length for dN/ p' \
    | sed 1,4d \
    | head -n -2 \
    | awk '{ SUM += $6} END { print SUM }')
  dS=$(cat $output_file |\
    sed -n -e '/dN & dS for each branch/,/tree length for dN/ p' \
    | sed 1,4d \
    | head -n -2 \
    | awk '{ SUM += $7} END { print SUM }')

  # output values
  echo "${basename%.*.*};codeml;$dnds_codeml;$dN;$dS"
done
