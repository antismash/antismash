#!/usr/bin/env python3
"""
Convenience wrapper to run antiSMASH directly from the source tree without installing.
"""

from antismash.__main__ import entrypoint

if __name__ == '__main__':
    entrypoint()
