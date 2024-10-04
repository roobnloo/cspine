FROM docker.io/rocker/verse:4
RUN /rocker_scripts/setup_R.sh
RUN /rocker_scripts/default_user.sh

USER root

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    file \
    curl \
    vim \
    python3 \
    python3-pip \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/* /tmp/library-scripts

RUN pip3 install -U radian

ENV GITHUB_CLI_VERSION 2.58.0
RUN wget https://github.com/cli/cli/releases/download/v${GITHUB_CLI_VERSION}/gh_${GITHUB_CLI_VERSION}_linux_amd64.tar.gz -O - | \
    tar xvzf - -C /usr/bin gh_${GITHUB_CLI_VERSION}_linux_amd64/bin/gh --strip-components=2

USER rstudio
RUN R -e "dotR <- file.path(Sys.getenv('HOME'), '.R'); dir.create(dotR); Makevars <- file.path(dotR, 'Makevars'); file.create(Makevars); cat('\nCXXFLAGS=-O3 -fPIC -Wno-unused-variable -Wno-unused-function -Wno-ignored-attributes', 'CXX = g++', 'CXX11 = g++', file = Makevars, sep = '\n', append = TRUE)"
RUN R -e "dotRprofile <- file.path(Sys.getenv('HOME'), '.Rprofile'); file.create(dotRprofile); cat('local({r <- getOption(\"repos\")', 'r[\"CRAN\"] <- \"https://cloud.r-project.org\"', 'options(repos=r, pkg.build_extra_flags=FALSE)',  'library(devtools)', '})', file = dotRprofile, sep = '\n', append = TRUE)"
RUN install2.r sparsegl devtools languageserver httpgd tictoc testthat igraph
