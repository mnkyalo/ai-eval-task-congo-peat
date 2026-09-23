#!/bin/bash
# Verifier entry point. Runs from /tests, which is baked into the verifier
# image; the agent's container is torn down before this runs.
#
# Anything the agent may have written under /logs is destroyed first. A reward
# pre-seeded by the agent must never survive to be read, whether or not the
# test run that follows succeeds.
#
# Deliberately no `set -e`: a failing or crashing pytest run must still fall
# through to the else branch and write a reward.

rm -rf /logs/verifier
mkdir -p /logs/verifier

# `python3 -m pytest`, never the bare console script -- if that resolved to a
# different interpreter than pytest-json-ctrf was installed into, pytest dies
# with "unrecognized arguments: --ctrf" and a correct submission scores 0.
python3 -m pytest --ctrf /logs/verifier/ctrf.json /tests/test_outputs.py -rA
EXIT_CODE=$?

if [ "$EXIT_CODE" -eq 0 ]; then
  echo 1 > /logs/verifier/reward.txt
else
  echo 0 > /logs/verifier/reward.txt
fi

exit 0
