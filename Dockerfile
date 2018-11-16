FROM nfcore/base
LABEL description="Docker image containing all requirements for nf-core/plant-rnaseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-plant-rnaseq-1.0dev/bin:$PATH
