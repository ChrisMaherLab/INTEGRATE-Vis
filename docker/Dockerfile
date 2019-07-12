FROM ubuntu:trusty

# get linux fundamentals
RUN apt-get update \
   && apt-get -y install build-essential libboost-all-dev \
   libncurses5-dev libncursesw5-dev libbz2-dev liblzma-dev cmake libfreetype6-dev pkg-config \
   git python2.7 python-pip wget vim \
   samtools bedtools

# install required python modules
RUN pip install matplotlib==2.0.2

# install gtfToGenePred
RUN wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred -O /usr/local/bin/gtfToGenePred \
   && chmod +x /usr/local/bin/gtfToGenePred

# setup integrate vis
RUN cd / \
   && git clone https://github.com/ChrisMaherLab/INTEGRATE-Vis.git \
   && cd INTEGRATE-Vis/INTEGRATE-Vis.1.0.0 \
   && chmod +x install.sh \
   && ./install.sh -o /usr/local/bin
