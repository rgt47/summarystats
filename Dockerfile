# syntax=docker/dockerfile:1.4
# zzcollab Dockerfile v2.8.0

ARG BASE_IMAGE=rocker/tidyverse
ARG R_VERSION=4.6.0
ARG USERNAME=analyst

FROM rocker/tidyverse:4.6.0

ARG USERNAME=analyst
ARG DEBIAN_FRONTEND=noninteractive

# RENV_PATHS_LIBRARY is outside the project bind-mount so the baked library
# is not shadowed at runtime. ZZCOLLAB_AUTO_RESTORE=false disables the
# startup restore so the image library is authoritative.
ENV LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 TZ=UTC \
    RENV_PATHS_LIBRARY=/opt/renv/library \
    RENV_PATHS_CACHE=/opt/renv/cache \
    RENV_CONFIG_REPOS_OVERRIDE="https://packagemanager.posit.co/cran/__linux__/noble/2026-06-09" \
    ZZCOLLAB_CONTAINER=true \
    ZZCOLLAB_AUTO_RESTORE=false

# No additional system dependencies required

# Configure R to use Posit Package Manager for pre-compiled binaries
RUN echo 'options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/noble/2026-06-09"))' \
        >> /usr/local/lib/R/etc/Rprofile.site && \
    echo 'options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"])))' \
        >> /usr/local/lib/R/etc/Rprofile.site

# Install languageserver for IDE support and yaml for R Markdown dependencies
RUN R -e "install.packages(c('languageserver', 'yaml'))"

# Install renv and restore packages from lockfile (using PPM binaries).
RUN R -e "install.packages('renv')"
RUN mkdir -p /opt/renv/library /opt/renv/cache && chmod 755 /opt/renv/library /opt/renv/cache
COPY renv.lock renv.lock
# renv::init creates the platform-specific library directory structure that
# renv::restore() requires to link packages from the cache.
RUN R -e "renv::init(bare=TRUE, force=TRUE, restart=FALSE); renv::restore()"

# Create non-root user
RUN useradd --create-home --shell /bin/bash ${USERNAME} && \
    chown -R ${USERNAME}:${USERNAME} /usr/local/lib/R/site-library

USER ${USERNAME}
WORKDIR /home/${USERNAME}/project

CMD ["R", "--quiet"]
