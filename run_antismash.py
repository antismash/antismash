#!/usr/bin/env python3
"""
Convenience wrapper to run antiSMASH directly from the source tree without installing.
"""
import sys
# wrap the import in a try/except as the import is lengthy and SIGINT is noisy
# when it happens during the import process
try:
    from antismash.__main__ import entrypoint
except KeyboardInterrupt:
    sys.exit(2)

if __name__ == '__main__':
    entrypoint()  # already wrapped in a graceful SIGINT try/except
