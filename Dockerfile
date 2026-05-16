# syntax=docker/dockerfile:1.4
# zzcollab manuscript-package profile
#======================================================================
# Purpose: R package development + Rmd manuscript PDF rendering
# Base:    rocker/verse:4.4.2 (R + tidyverse + pandoc + tinytex)
# Profile: manuscript-package
#
# Build: DOCKER_BUILDKIT=1 docker build -t compendium-env .
# Run:   docker run --rm -v $(pwd):/project -w /project compendium-env Rscript -e '...'
#======================================================================

FROM rocker/verse:4.4.2

ARG DEBIAN_FRONTEND=noninteractive

ENV LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8 \
    TZ=UTC \
    OMP_NUM_THREADS=1 \
    RENV_PATHS_CACHE=/root/.cache/R/renv \
    RENV_CONFIG_REPOS_OVERRIDE="https://packagemanager.posit.co/cran/__linux__/noble/latest" \
    ZZCOLLAB_CONTAINER=true

# jq: used by zzcollab CI scripts
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt/lists,sharing=locked \
    apt-get update && \
    apt-get install -y --no-install-recommends jq

RUN R -e "install.packages('renv')"

# LaTeX packages required by manuscripts.
# Must be installed as root; tlmgr auto-install fails for non-root processes.
RUN R -e "tinytex::tlmgr_install(c('unicode-math', 'fontspec', 'euenc', 'l3packages', 'tipa', 'xunicode', 'lineno'))"

# Pre-install project package library into the image layer.
COPY renv.lock renv.lock
RUN Rscript -e "renv::restore(prompt = FALSE)"

WORKDIR /project

CMD ["R", "--quiet"]
