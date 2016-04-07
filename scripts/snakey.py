#!/usr/bin/env python3

import snakemake
import sys
import os
import argparse
import datetime
import time
import yaml

from snakemake.utils import min_version

min_version("3.5")


ts = time.time()
now = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H:%M:%S')

snakey_dir = os.path.dirname(os.path.realpath(__file__))

parser=argparse.ArgumentParser("A general purpose script for running snakemake in a scheduled environment.")
parser.add_argument("snakefile",
                    help="Snakefile to use.")
parser.add_argument("config",
                    help="Configuration file to use.")
parser.add_argument("-c", "--hpc_conf", default=snakey_dir+"/hpc.json",
                    help="Configuration file that allows the user to override resource requests for each rule when running under a scheduler in a HPC environment.")
parser.add_argument("-N", "--max_nodes", type=int, default="10",
                    help="Maximum number of nodes to use concurrently")
parser.add_argument("-n", "--max_cores", type=int, default="1000",
                    help="Maximum number of cores to use concurrently")
parser.add_argument("-d", "--no_drmaa", action='store_true', default=False,
                    help="Use this flag if running on a HPC and DRMAA is not available")
parser.add_argument("--force_incomplete", action='store_true', default=False,
                    help="Force snakemake to rerun incomplete steps")
parser.add_argument("--detailed_summary", action='store_true', default=False,
                    help="Print detailed summary of all input and output files")
parser.add_argument("--list_resources", action='store_true', default=False,
                    help="List resources used in the workflow")
parser.add_argument("--make_dag", action='store_true', default=False,
                    help="Produce a DAG rather than execute the workflow")
args=parser.parse_args()


with open(args.config, 'r') as f:
        doc = yaml.load(f)


SCHEDULER=doc["scheduler"] if doc["scheduler"] else ""
CWD=os.path.abspath(".")

QUEUE=doc["queue"]

res_cmd = ""
sub_cmd = ""

if not os.path.exists("logs"):
    os.makedirs("logs")

if SCHEDULER == "LSF":
    sub_cmd = "bsub"
    res_cmd = " -R rusage[mem={cluster.memory}]span[ptile={threads}] -n {threads} -q " + QUEUE + " -J {rule} -oo logs/" + now + "_{rule}_%j.out"
elif SCHEDULER == "PBS":
    sub_cmd = "qsub"
    res_cmd = " -lselect=1:mem={cluster.memory}MB:ncpus={threads} -q " + QUEUE + " -N {rule} -o logs/{rule}_%j.out -e logs/{rule}_%j.out"
elif SCHEDULER == "SLURM":
    sub_cmd = "sbatch"
    cores = "-c {threads} "
    res_cmd = " " + cores + "-p " + QUEUE + " --mem={cluster.memory} -J {rule} -o logs/" + now + "_{rule}_%j.out -e logs/" + now + "_{rule}_%j.out"


sf = args.snakefile if os.path.exists(args.snakefile) else snakey_dir + "/" + args.snakefile


snakemake.snakemake(sf,
                    cores=args.max_nodes,
                    nodes=args.max_cores,
                    configfile=args.config,
                    workdir=".",
                    cluster_config=args.hpc_conf,
                    cluster=sub_cmd + res_cmd if args.no_drmaa else None,
                    drmaa=res_cmd if not args.no_drmaa else None,
                    printshellcmds=True,
                    snakemakepath="/tgac/software/testing/python_anaconda/2.4.1_dm/x86_64/bin/snakemake",
                    stats="snakemake_" + now + ".stats",
                    force_incomplete=args.force_incomplete,
                    detailed_summary=args.detailed_summary,
                    list_resources=args.list_resources,
                    latency_wait=60 if not SCHEDULER == "" else 1,
                    printdag=args.make_dag,
                    forceall=args.make_dag
                    )
