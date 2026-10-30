FROM mambaorg/micromamba:latest
COPY --chown=$MAMBA_USER:$MAMBA_USER envs/wgs-gatk-mem3.yaml /tmp/environment.yaml
RUN micromamba create -y -f /tmp/environment.yaml && micromamba clean --all --yes
COPY --chown=$MAMBA_USER:$MAMBA_USER . /app
WORKDIR /app
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN micromamba run -n wgs-gatk-mem3 python -m pip install --no-deps -e .
ENTRYPOINT ["micromamba", "run", "-n", "wgs-gatk-mem3", "./run_pipeline.sh"]
