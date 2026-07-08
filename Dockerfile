# Environment for the downstream RNA-seq analysis of the S. canicula dental lamina.
#
# The base image provides R 4.3 and the system libraries needed to build the
# Bioconductor 3.18 packages; the exact package set (R 4.3.3 / Bioconductor 3.18,
# 169 packages) is then restored from renv.lock so the image matches the
# environment used to generate the results in suppl_files/dea_output/.
FROM bioconductor/bioconductor_docker:RELEASE_3_18

LABEL org.opencontainers.image.description="Downstream RNA-seq analysis environment for the Scyliorhinus canicula dental lamina project (R 4.3 / Bioconductor 3.18)"

# The analysis directory is mounted at run time and carries an renv project
# (.Rprofile + renv/) from the host. Disable renv's autoloader so the container
# uses the packages installed in the image library rather than trying to activate
# the host's (wrong-platform, empty) project library.
ENV RENV_CONFIG_AUTOLOADER_ENABLED=FALSE

WORKDIR /project

# Restore the pinned package set from renv.lock into the image's default library,
# so the packages are available without renv project activation at run time.
RUN Rscript -e "install.packages('renv', repos = 'https://packagemanager.posit.co/cran/latest')"
COPY renv.lock renv.lock
RUN Rscript -e "renv::restore(lockfile = 'renv.lock', library = .libPaths()[1], prompt = FALSE)"
