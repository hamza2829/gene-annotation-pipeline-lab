FROM python:3.11-slim
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 MPLCONFIGDIR=/tmp/matplotlib HOME=/tmp
WORKDIR /app
COPY pyproject.toml ./
COPY src ./src
RUN pip install --no-cache-dir . && useradd --uid 10001 --create-home lab
COPY --chown=10001:10001 config.yaml ./
COPY --chown=10001:10001 data ./data
RUN mkdir -p /app/artifacts && chown 10001:10001 /app/artifacts
USER 10001
CMD ["gene-lab", "all", "--config", "config.yaml", "--genes", "data/genes.txt", "--offline"]
