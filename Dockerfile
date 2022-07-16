FROM continuumio/miniconda3:4.10.3
LABEL author="Paul Cantalupo" \
      description="Docker image containing R and R packages necessary for running TILscore"


# Install procps so that Nextflow can poll CPU usage and
# deep clean the apt cache to reduce image/layer size
RUN apt-get update \
      && apt-get install -y procps \
      && apt-get clean -y && rm -rf /var/lib/apt/lists/*

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/tilscore/bin:$PATH

ENV R_LIBS_USER=/opt/conda/envs/tilscore/lib/R/library

