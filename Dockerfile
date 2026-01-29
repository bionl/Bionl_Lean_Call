# Start with Python 3.11 base image
FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Install system dependencies including procps for Nextflow
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    make \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    zlib1g-dev \
    libssl-dev \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip
RUN pip install --upgrade pip setuptools wheel

# Install Python packages
RUN pip install --no-cache-dir \
    numpy \
    pandas \
    cython \
    cyvcf2 \
    openpyxl
# Install samtools and bcftools
RUN apt-get update && apt-get install -y samtools bcftools && rm -rf /var/lib/apt/lists/*
# Verify installations
RUN python -c "import pandas; import cyvcf2; print('All packages installed successfully')"

# Verify ps command exists
RUN ps --version

# Set the default command
CMD ["python"]