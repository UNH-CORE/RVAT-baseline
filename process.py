#!/usr/bin/env python
"""This script processes data from the RVAT baseline experiment."""

import pyrvatbl.processing as pr

if __name__ == "__main__":
    pr.batchperf()
    pr.batchwake()
