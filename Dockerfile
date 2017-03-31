FROM philippmuench/fastcodeml:1.1.0

# install octave
RUN add-apt-repository ppa:octave/stable
RUN apt-get update
RUN apt-get install -y octave

# install PAML
RUN mkdir -p /usr/src/paml \
  && curl -SL "http://abacus.gene.ucl.ac.uk/software/paml4.9c.tgz" \
  | tar zxC /usr/src/paml \
  && cd /usr/src/paml/paml4.9c/src \
  && make -j"$(nproc)" \
  && mv codeml /usr/bin/ \
  && mv baseml /usr/bin/ \
  && mv basemlg /usr/bin/ \
  && mv chi2 /usr/bin/ \
  && mv evolver /usr/bin/ \
  && mv infinitesites /usr/bin/ \
  && mv mcmctree /usr/bin/ \
  && mv pamp /usr/bin/ \
  && mv yn00 /usr/bin/ \
  && rm -rf /usr/src/paml

