# Optional alternative to Conda on the host; SLURM does NOT invoke Docker directly.
FROM condaforge/miniforge3:latest
COPY envs/wgs.yaml /tmp/wgs.yaml
RUN conda env create -f /tmp/wgs.yaml && conda clean -afy
ENV PATH="/opt/conda/envs/wgs-gatk/bin:${PATH}"
WORKDIR /pipeline
COPY . /pipeline
CMD ["bash", "run_pipeline.sh", "--help"]
