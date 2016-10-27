#!/usr/bin/env python3

import re
import sys
import json
import argparse
import os
import csv
from copy import copy


def to_csv(string):
	if not os.path.exists(string) or not os.path.isfile(string):
		raise OSError("The input file does not exist!")
	return csv.DictReader(open(string), delimiter="\t")


def to_csv_out(handle, fieldnames):
	writer = csv.DictWriter(handle, fieldnames, delimiter="\t")
	writer.writeheader()
	return writer


def to_json(handle, fieldnames, debug=False):
	'''This function will load the JSON configuration, check that
	both parameters and expression are defined, and prepare the configuration for the analysis.
	Specifically, it will transform the operators/values into predefined functions
	that will be used for the relative value in the row.
	Thus, a configuration of:
	"param": {
	   "operator": "eq",
	   "value": 1
	   }

	will yield a function that checks x==1.
	All functions are in the boolean space.

	The expression instead is a boolean expression which will evaluate the truth status
	of the row given the results of the boolean expressions for each parameter.

	'''

	json_dict = json.load(handle)
	if "parameters" not in json_dict or "expression" not in json_dict:
		print(
			"Configuration is faulty - please ensure that the JSON has valid \"parameters\" and \"expression\" fields.",
			sys.stderr)
		sys.exit(1)

	for param in set(json_dict["parameters"]):
		parameter_name = param.split(".")[0]
		if json_dict["parameters"][param]["operator"] not in ("gt", "ge", "eq", "lt", "le", "in", "not in"):
			print("Unrecognized operator for {1}: {0}".format(json_dict["parameters"][param]["operator"], param),
				  sys.stderr)
			sys.exit(1)

		json_dict["parameters"][param]["name"] = parameter_name

	parameter_names = set(json_dict["parameters"][x]["name"] for x in json_dict["parameters"])

	diff = set.difference(parameter_names, set(fieldnames))
	if len(diff) > 0:
		print("Unrecognized parameters: {0}".format(",".join(list(diff))),
			  sys.stderr)
		print("Fieldnames:\n\t{0}".format("\n\t".join(fieldnames)), file=sys.stderr)
		print("Parameter names:\n\t{0}".format("\n\t".join(parameter_names)), file=sys.stderr)
		sys.exit(1)

	keys = list(filter(lambda x: x not in ("and", "or"), re.findall("([^ ()]+)", json_dict["expression"])))
	diff_params = set.difference(set(keys), set(json_dict["parameters"].keys()))
	if len(diff_params) > 0:
		print("Expression and required parameters mismatch:\n\t{0}".format("\n\t".join(list(diff_params))),
			  sys.stderr)
		sys.exit(1)

	newexpr = json_dict["expression"][:]
	for key in keys:
		newexpr = re.sub(key, "evaluated[\"{0}\"]".format(key), newexpr)
	json_dict["expression"] = newexpr
	newexpr = compile(newexpr, "<string>", "eval")
	if debug is True:
		print(json.dumps(json_dict), file=sys.stderr)

	json_dict["compiled"] = newexpr

	return json_dict


def evaluate(param, conf):
	'''This function evaluates the truth status of a row parameter given the configuration present in the JSON file.'''

	if conf["operator"] == "eq":
		return float(param) == float(conf["value"])
	elif conf["operator"] == "gt":
		return float(param) > float(conf["value"])
	elif conf["operator"] == "lt":
		return float(param) < float(conf["value"])
	elif conf["operator"] == "ge":
		return float(param) >= float(conf["value"])
	elif conf["operator"] == "le":
		return float(param) <= float(conf["value"])
	elif conf["operator"] == "in":
		return param in conf["value"]
	elif conf["operator"] == "not in":
		return param not in conf["value"]


def main():
	parser = argparse.ArgumentParser("Script to automate CSV filtering based on a JSON configuration.")
	parser.add_argument("--json", type=argparse.FileType('r'), required=True,
						help="JSON configuration file.")
	parser.add_argument("--input", type=to_csv, required=True,
						help="CSV junction file.")  #
	parser.add_argument("--debug", action='store_true', default=False,
						help="Debug flag")
	parser.add_argument("--out", type=argparse.FileType('w'), required=True,
						help="CSV output junction file.")
	args = parser.parse_args()

	args.out = to_csv_out(args.out, args.input.fieldnames)

	json_dict = to_json(args.json, args.input.fieldnames, debug=args.debug)

	for row in args.input:
		evaluated = {}
		for param in json_dict["parameters"]:
			name = json_dict["parameters"][param]["name"]
			evaluated[param] = evaluate(row[name], json_dict["parameters"][param])

			# json_dict["parameters"][param]["function"](row[param])
		if eval(json_dict['compiled']) is True:
			args.out.writerow(row)

	sys.exit(0)


if __name__ == '__main__': main()
