# Use Ubuntu 22.04 LTS as the base image
FROM ubuntu:22.04

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# 1. Install necessary system dependencies
# These packages are equivalent to the required conda dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    python3 \
    python3-dev \
    python3-pip \
    wget \
    autoconf \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libncurses5-dev \
    && rm -rf /var/lib/apt/lists/*

# 2. Install Python packages from requirements.txt
COPY requirements.txt .
RUN pip3 install --no-cache-dir -r requirements.txt

# 3. Build HTSlib 1.23.1 from source
# Matches the htslib version used in your conda environment
RUN wget https://github.com/samtools/htslib/releases/download/1.23.1/htslib-1.23.1.tar.bz2 && \
    tar xjf htslib-1.23.1.tar.bz2 && \
    cd htslib-1.23.1 && \
    ./configure && \
    make && \
    make install && \
    ldconfig && \
    cd .. && rm -rf htslib-1.23.1*

# Set working directory and copy project files
WORKDIR /app
COPY . .

# Compile the project
RUN make PYTHON=python3

# Set the default command to run Platypus
CMD ["python3", "Platypus.py"]