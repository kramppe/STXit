FROM mambaorg/micromamba:1.4.9

# Set working directory
WORKDIR /app

# Copy environment file
COPY environment.yml /tmp/environment.yml

# Install conda environment
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Set PATH to include conda environment
ENV PATH /opt/conda/bin:$PATH

# Copy STXit source code
COPY stxit/ ./stxit/
COPY databases/ ./databases/
COPY setup.py ./
COPY README.md ./

# Install STXit
RUN pip install -e .

# Create output directory
RUN mkdir -p /data/output

# Set environment variables
ENV STXIT_DB_PATH=/app/databases/stx_references.fasta
ENV PYTHONPATH=/app:$PYTHONPATH

# Default command
CMD ["stxit", "--help"]

# Volume for data
VOLUME ["/data"]

# Expose port for potential web interface
EXPOSE 8080

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
  CMD stxit --version || exit 1

LABEL maintainer="STXit Team <stxit@tool.dev>"
LABEL version="1.0.3"
LABEL description="STXit - Shiga Toxin Detection and Analysis Pipeline"
