FROM neuralensemble/base
MAINTAINER andrew.davison@unic.cnrs-gif.fr

ENV NRN_VER=7.4
ENV NRN=nrn-$NRN_VER
ENV PATH=$PATH:$VENV/bin
RUN ln -s /usr/bin/2to3-3.4 $VENV/bin/2to3

WORKDIR $HOME/packages
RUN wget http://www.neuron.yale.edu/ftp/neuron/versions/v$NRN_VER/$NRN.tar.gz
RUN tar xzf $NRN.tar.gz; rm $NRN.tar.gz
RUN git clone --depth 1 https://github.com/INCF/libneurosim.git
RUN cd libneurosim; ./autogen.sh

RUN mkdir $VENV/build
WORKDIR $VENV/build
RUN mkdir libneurosim; \
    cd libneurosim; \
    PYTHON=$VENV/bin/python $HOME/packages/libneurosim/configure --prefix=$VENV; \
    make; make install; ls $VENV/lib $VENV/include
RUN mkdir $NRN; \
    cd $NRN; \
    $HOME/packages/$NRN/configure --with-paranrn --with-nrnpython=$VENV/bin/python --disable-rx3d --without-iv --prefix=$VENV; \
    make; make install; \
    cd src/nrnpython; $VENV/bin/python setup.py install; \
    cd $VENV/bin; ln -s ../x86_64/bin/nrnivmodl

RUN $VENV/bin/pip3 install lazyarray nrnutils PyNN==0.9.1 jupyter

WORKDIR /home/docker/
RUN echo "source $VENV/bin/activate" >> .bashrc

WORKDIR /home/docker/packages
RUN git clone https://github.com/DavidTingley/papers

WORKDIR /home/docker/packages/papers/LS_phasecoding/modeling
