# Verifier image. Built and run after the agent's container is torn down
# (environment_mode = "separate"), so everything the grading code needs must be
# present here. test.sh installs nothing at verify time.

FROM python:3.13-slim-bookworm

# numpy is needed because the verifier re-derives the answer from the sealed
# raw data rather than trusting a stored key.
RUN pip install --no-cache-dir pytest==9.1.1 pytest-json-ctrf==0.5.2 numpy==2.2.6

# Bake the grading code in; Harbor does not upload tests/ at verify time.
COPY . /tests/

# Parent directory of every path in task.toml's top-level artifacts list.
RUN mkdir -p /app
