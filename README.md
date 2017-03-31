# selection_docker
calculation of global dN/dS based on codeml, pairwise codeml. Inside the docker you can also run fastcodeml with the `fast` command. 

# usage

this script will run dN/dS analysis on `codon.phy` and the corresponding tree `tree.phy` and output the dN/dS

```
sudo docker bild -t selection_docker . 
sudo docker run -i -v /home/pmuench/github.com/philippmuench/codeml_docker:/data -t --entrypoint /bin/bash selection_docker
./codeml_example.sh
```
