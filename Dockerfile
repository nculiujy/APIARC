FROM mambaorg/micromamba:1.5.8

LABEL org.opencontainers.image.title="APIARC" \
      org.opencontainers.image.description="ChIP-seq / RNA-seq integration pipeline with explicit statistical inference" \
      org.opencontainers.image.source="https://github.com/nculiujy/APIARC"

USER root

RUN apt-get update \
 && apt-get install -y --no-install-recommends git procps file ca-certificates \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/apiarc

COPY environment.yml /opt/apiarc/environment.yml
RUN micromamba install -y -n base -f /opt/apiarc/environment.yml \
 && micromamba clean --all --yes

ENV PATH=/opt/conda/bin:$PATH \
    PYTHONUNBUFFERED=1 \
    MAMBA_ROOT_PREFIX=/opt/conda

COPY --chown=micromamba:micromamba Snakefile /opt/apiarc/Snakefile
COPY --chown=micromamba:micromamba config /opt/apiarc/config
COPY --chown=micromamba:micromamba workflow /opt/apiarc/workflow

VOLUME ["/opt/apiarc/workflow/resources", "/opt/apiarc/result"]

USER micromamba

RUN cp -a /opt/apiarc /tmp/smoke \
 && cd /tmp/smoke \
 && python workflow/scripts/common/stub_resources.py \
 && snakemake -n --quiet \
 && rm -rf /tmp/smoke

CMD ["snakemake", "-n"]
