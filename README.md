# selection_docker
calculation of global dN/dS based on codeml, pairwise codeml. Inside the docker you can also run fastcodeml with the `fast` command. 

### implemented methods

- PAML (codeml)
- slimcodeml
- fastcodeml
- yn00 (should be installed but not inlcluded in `codeml_example.sh`)

### input file specification

see 'tree.phy' and 'codon.phy' for example input format
trees have to be unrooted. You may want to use 'unroot.R'
for each run a config file have to be provided, see 'codeml.ctl'(for codeml) and 'codeml_pairwise.ctl' (for pairwise codeml)

### usage

this script will run dN/dS analysis on `codon.phy` and the corresponding tree `tree.phy` and output the dN/dS

##### run the scripts
```
sudo docker run hello-world # test if docker is installed
git clone https://github.com/hzi-bifo/selection_docker.git #clone
cd selection_docker 
sudo docker build -t selection_docker . # Should be finished after 3-5 minutes with [...] Successfully built ...
sudo docker run -i -v /absolute_path/to/selection_docker/folder:/data -t --entrypoint /bin/bash selection_docker # enter the docker image
cd data/ 
./codeml_example.sh
```

### output

- dN/dS will be printed to stdout
- detailed output will be written to `test.mlc` and `test_pair.mlc`
