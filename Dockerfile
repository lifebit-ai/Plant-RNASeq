FROM nfcore/base

LABEL authors="phil@lifebit.ai" \
      description="Docker image containing requirements for lifebit-ai/plant-rnaseq pipeline fastqc, REAL & multiqc processes"

## Get those environment.yml goodies for REAL
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-plant-rnaseq-1.0dev/bin:$PATH

COPY ./bin/noReads.py /usr/local/bin/
COPY ./bin/geneExp.py /usr/local/bin/

RUN chmod +x //usr/local/bin/noReads.py
RUN chmod +x //usr/local/bin/geneExp.py
