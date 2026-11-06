FROM rocker/r-ver:4.4.1

RUN apt-get update && apt-get install -y --no-install-recommends \
      libcurl4-openssl-dev libssl-dev libxml2-dev \
      libnlopt-dev gfortran libblas-dev liblapack-dev \
      git \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY renv.lock renv.lock
COPY .Rprofile .Rprofile
COPY renv/ renv/
RUN R -e "install.packages('renv', repos='https://cloud.r-project.org'); renv::restore(prompt=FALSE)"

COPY api/ api/
COPY R/ R/
COPY config.yml config.yml

# latest.rds is mounted at runtime from the pipeline artifact, NOT baked in.
# If you want a self-contained image with the current model, add:
#   COPY outputs/latest.rds outputs/latest.rds
# But that couples the image to a point-in-time model, which is what the
# artifact pattern exists to avoid.

EXPOSE 8000
CMD ["Rscript", "api/run.R"]