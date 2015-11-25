#!/usr/bin/env python3

import snakemake
import sys
import os
import argparse
import datetime
import time

ts = time.time()
now = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H:%M:%S')

parser=argparse.ArgumentParser("Script to generate and run portcullis, and other similar tools against a genome.")
parser.add_argument("config",
                    help="Configuration file to use for running pp.")
parser.add_argument("--force_incomplete",
                    help="Force snakemake to rerun incomplete steps")
args=parser.parse_args()


this_dir = os.path.dirname(os.path.realpath(__file__))


snakemake.snakemake(this_dir + "/read_gen.snakefile",
                    cores=100,
                    nodes=10,
                    configfile=args.config,
                    workdir=".",
                    cluster_config=this_dir + "/hpc.json",
                    drmaa=" -R rusage[mem={cluster.memory}]span[ptile={threads}] -n {threads} -q {cluster.queue} -oo /dev/null",
                    printshellcmds=True,
                    snakemakepath="/tgac/software/testing/python/3.4.2/x86_64/bin/snakemake",
                    stats="pp_" + now + ".stats",
                    force_incomplete=args.force_incomplete,
                    latency_wait=30
                    )
