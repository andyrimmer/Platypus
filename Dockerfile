FROM ubuntu:14.04

RUN apt-get update && apt-get install -y build-essential \
			python-dev \
			python-pip \
			git \
			wget \
			autoconf \
			zlib1g-dev

# Make a working directory
RUN mkdir /code

ADD . /code/Platypus

# Install htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.3/htslib-1.3.tar.bz2
RUN tar xjf htslib-1.3.tar.bz2
RUN cd htslib-1.3 && autoconf && ./configure && make && make install

# Install Cython
RUN pip install cython

ENV C_INCLUDE_PATH /usr/local/include
ENV LIBRARY_PATH /usr/local/lib
ENV LD_LIBRARY_PATH /usr/local/lib

# Compile
RUN cd /code/Platypus && make
