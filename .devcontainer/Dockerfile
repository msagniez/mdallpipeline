# Dev container for nf-core-mdallpipeline
FROM mcr.microsoft.com/devcontainers/base:ubuntu-22.04

ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=$USER_UID

ENV DEBIAN_FRONTEND=noninteractive
ENV JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64
ENV PATH=$JAVA_HOME/bin:$PATH

RUN apt-get update \
    && apt-get -y install --no-install-recommends \
        openjdk-17-jdk \
        curl \
        ca-certificates \
        git \
        wget \
        jq \
        unzip \
        build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install Nextflow
RUN curl -fsSL https://get.nextflow.io | bash -s -- -q \
    && mv nextflow /usr/local/bin/nextflow \
    && chmod +x /usr/local/bin/nextflow

# Quick verification
RUN java -version && nextflow -version

# Use non-root developer user from base image
USER ${USERNAME}
