# Agent environment for the Core MBK-2 pollen influx task.
# Contains only what a solver legitimately gets: the core dataset, the two
# atmospheric calibration curves, and a Python toolchain. It must never
# reference anything under solution/ or tests/.

FROM python:3.11-slim-bookworm

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir numpy==2.2.6

WORKDIR /app
COPY data/ /app/data/
