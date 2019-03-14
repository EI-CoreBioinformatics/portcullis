#!/usr/bin/env python3

import argparse
import json
import re
import os
import sys
from pandas import DataFrame
import pandas as pd

try:
	from .performance import Performance
except:
	from performance import Performance

def replace_op(op):
    if op == "eq":
        return "=="
    elif op == "gt":
        return ">"
    elif op == "lt":
        return "<"
    elif op == "gte":
        return ">="
    elif op == "lte":
        return "<="
    elif op == "in":
        return ".isin("
    elif op == "not in":
        return ".isin("

def load_genuine(genuine_file):
	glist = []
	with open(genuine_file) as gin:
		for line in gin:
			glist.append(bool(line))
	return glist


def json2pandas(handle, fieldnames, data_frame):
	'''This function will load the JSON configuration, check that
	both parameters and expression are defined, and convert to a pandas string that will operate on a dataframe
	call df
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
		raise ValueError(
			"Configuration is faulty - please ensure that the JSON has valid \"parameters\" and \"expression\" fields.")

	for param in set(json_dict["parameters"]):
		parameter_name = param.split(".")[0]
		if json_dict["parameters"][param]["operator"] not in ("gt", "gte", "eq", "lt", "lte", "in", "not in"):
			raise ValueError("Unrecognized operator for {1}: {0}".format(json_dict["parameters"][param]["operator"], param))

		json_dict["parameters"][param]["name"] = parameter_name

	parameter_names = set(json_dict["parameters"][x]["name"] for x in json_dict["parameters"])

	diff = set.difference(parameter_names, set(fieldnames))
	if len(diff) > 0:
		raise ValueError("Unrecognized parameters: {0}".format(",".join(list(diff))) + "\n" +
						"Fieldnames:\n\t{0}".format("\n\t".join(fieldnames)) + "\n" +
						"Parameter names:\n\t{0}".format("\n\t".join(parameter_names)))

	keys = list(filter(lambda x: x not in ("&", "|"), re.findall("([^ ()]+)", json_dict["expression"])))
	diff_params = set.difference(set(keys), set(json_dict["parameters"].keys()))
	if len(diff_params) > 0:
		raise ValueError("Expression and required parameters mismatch:\n\t{0}".format("\n\t".join(list(diff_params))))

	newexpr = json_dict["expression"][:]
	for key in keys:
		k = key[:-2] if key[-2] == '.' and key[-1].isdigit() else key
		if json_dict['parameters'][key]['operator'] == "in":
			newexpr = re.sub(key, "({3}[\"{0}\"].isin({2}))".format(k, replace_op(json_dict['parameters'][key]['operator']),
															   json_dict['parameters'][key]['value'], data_frame), newexpr)
		elif json_dict['parameters'][key]['operator'] == "not in":
			newexpr = re.sub(key,
							 "(~{3}[\"{0}\"].isin({2}))".format(k, replace_op(json_dict['parameters'][key]['operator']),
															  json_dict['parameters'][key]['value'], data_frame), newexpr)
		else:
			newexpr = re.sub(key, "({3}[\"{0}\"] {1} {2})".format(k, replace_op(json_dict['parameters'][key]['operator']),
															 json_dict['parameters'][key]['value'], data_frame), newexpr)

	return data_frame + ".loc[" + newexpr + "]", data_frame + ".loc[~" + newexpr + "]"


def calcPerformance(passed, failed, invert=False):
	tp = 0
	tn = 0
	fp = 0
	fn = 0
	passed_res = pd.value_counts(passed['genuine'].values, sort=True)
	failed_res = pd.value_counts(failed['genuine'].values, sort=True)
	if invert:
		tn = passed_res['False']
		fn = passed_res['True']
		tp = failed_res['True']
		fp = failed_res['False']
	else :
		tp = passed_res['True']
		fp = passed_res['False']
		tn = failed_res['False']
		fn = failed_res['True']

	return Performance(tp=tp, tn=tn, fp=fp, fn=fn)


def create_training_sets(args):

	# Load portcullis junctions into dataframe
	print("Loading input junctions ... ", end="", flush=True)
	original = DataFrame.read_csv(args.input, sep='\t', header=0)
	fieldnames = [key for key in dict(original.dtypes)]
	print("done.", len(original), "junctions loaded.")

	# Before we go further make sure we have a sufficent number of junctions to work with.  Minimum 1000.
	if len(original) < 500:
		raise ValueError("Not enough junctions to create training set")

	if args.genuine:
		glist = load_genuine(args.genuine)
		if len(glist) != len(original):
			raise ValueError("Genuine list and input junctions do not contain the same number of elements.  Genuine:" + len(glist) + "; input:" + len(original))

		original["genuine"] = pd.Series(glist, index=original.index)


	print()
	print("Creating initial positive set for training")
	print("------------------------------------------")
	print()
	print("Applying the following set of rule-based filters to create initial positive set.")
	for i, pjson in enumerate(args.pos_json, start=1):
		print("\t".join([str(i), pjson]))
	print()
	print("LAYER\t", end="")
	if args.genuine:
		print(Performance.longHeader())
	else:
		print("PASS\tFAIL")

	df = original.copy()	# Required for pandas eval
	pos_juncs = None

	# Run through layers of logic to get the positive set
	for i, json_file in enumerate(args.pos_json, start=1):

		print(str(i) + "\t", end="", flush=True)

		# Create pandas command
		pandas_cmd_in, pandas_cmd_out = json2pandas(open(json_file), fieldnames, "df")

		# print(pandas_cmd)

		# Execute the pandas command, result should be a filtered dataframe
		pos_juncs = eval(pandas_cmd_in)

		nb_not_pos = len(original) - len(pos_juncs)

		if args.genuine:
			not_pos_juncs = original.reset_index().merge(pos_juncs, indicator=True, how='outer').set_index('index')
			not_pos_juncs = not_pos_juncs.loc[not_pos_juncs['_merge'] == 'left_only']
			del not_pos_juncs['_merge']
			print(calcPerformance(pos_juncs, not_pos_juncs).longStr())
		else:
			print(str(len(pos_juncs)) + "\t" + str(nb_not_pos))

		if args.save_layers:
			pos_juncs.to_csv(args.prefix + ".pos_layer_" + str(i) + ".tab", sep='\t')

		# Check we have enough junctions left in positive se (100), if not then stop here
		if len(pos_juncs) <= 100:
			print("WARNING: We recommend at least 100 junctions in the positive set and this set of rules lowered " + \
				  "the positive set to", len(pos_juncs), ".  Will not filter positive set further.", file=sys.stderr)
			pos_juncs = df.copy()	# Override previous filter
			break

		df = pos_juncs.copy()

	# Get L95 for intron sizes
	if len(pos_juncs) == 0:
		raise ValueError("Can't build training sets, positive set filter left no junctions remaining.")

	pos_intron_sizes = pos_juncs["size"].tolist()
	pos_intron_sizes.sort(key=int)
	L95_pos = int(len(pos_intron_sizes) * 0.95)
	L95 = pos_intron_sizes[L95_pos]
	pos_length_limit = int(L95 * 1.2)

	print("Intron size at L95 =", L95, " positive set maximum intron size limit set to L95 x 1.2:", pos_length_limit)

	# Also save this to file as we'll need it back in the C program
	with open(args.prefix + ".L95_intron_size.txt", 'w') as l95out:
		print("Length of intron at 95th percentile", file=l95out)
		print(L95, file=l95out)

	if len(pos_juncs) > 100:
		pos_juncs = pos_juncs.loc[pos_juncs["size"] <= pos_length_limit]
		print("\t".join([str(x) for x in [i+1, len(pos_juncs), len(original) - len(pos_juncs)]]))

		if args.save_layers:
			pos_juncs.to_csv(args.prefix + ".pos_layer_intronsize.tab", sep='\t')

	print()
	print("Positive set contains:", len(pos_juncs), "junctions")
	print()
	print("Saving positive set to disk ... ", end="", flush=True)
	if args.genuine:
		del pos_juncs["genuine"]
	pos_file = args.prefix + ".pos.junctions.tab"
	pos_juncs.to_csv(pos_file, sep='\t')
	print("done. File saved to:", pos_file)

	not_pos_juncs = original.reset_index().merge(pos_juncs, indicator=True, how='outer').set_index('index')
	not_pos_juncs = not_pos_juncs.loc[not_pos_juncs['_merge'] == 'left_only']
	del not_pos_juncs['_merge']

	print(len(not_pos_juncs), "remaining for consideration as negative set")

	print()
	print("Creating initial negative set for training")
	print("------------------------------------------")
	print("Applying a set of rule-based filters to create initial negative set.")
	for i, njson in enumerate(args.neg_json, start=1):
		print("\t".join([str(i), njson]))
	print()
	print("LAYER\t", end="")
	if args.genuine:
		print(Performance.longHeader())
	else:
		print("PASS\tFAIL")

	neg_set = None
	other_juncs = not_pos_juncs

	# Run through layers of logic to get the positive set
	for i, json_file in enumerate(args.neg_json, start=1):

		print(str(i) + "\t", end="", flush=True)

		# Create pandas command
		pandas_cmd_in, pandas_cmd_out = json2pandas(open(json_file), fieldnames, "other_juncs")

		#print(pandas_cmd)

		# Execute the pandas command, result should be a filtered dataframe
		neg_juncs = eval(pandas_cmd_in)
		other_juncs = eval(pandas_cmd_out)

		if i > 1:
			neg_set = pd.concat([neg_set, neg_juncs])
		else:
			neg_set = neg_juncs

		if args.genuine:
			print(calcPerformance(neg_juncs, df).longStr())
		else:
			print(str(len(neg_juncs)) + "\t" + str(len(other_juncs)))

		if args.save_layers:
			neg_juncs.to_csv(args.prefix + ".neg_layer_" + str(i) + ".tab", sep='\t')


	neg_length_limit = int(L95 * 8)
	print("Intron size L95 =", L95, "negative set will use junctions with intron size over L95 x 8:", neg_length_limit, "and with maxmmes < 12")
	neg_juncs = other_juncs.loc[other_juncs["size"] > neg_length_limit]
	neg_juncs = neg_juncs.loc[neg_juncs["maxmmes"] < 12]
	neg_set = pd.concat([neg_set, neg_juncs])
	if args.genuine:
		print(str(i+1) + "\t" + calcPerformance(neg_juncs, df).longStr())
	else:
		print(str(i+1) + "\t" + str(len(neg_juncs)) + "\t" + str(len(other_juncs)))

	if args.save_layers:
		neg_juncs.to_csv(args.prefix + ".neg_layer_intronsize.tab", sep='\t')


	print()
	print("Negative set contains:", len(neg_set), "junctions")

	if args.genuine:
		del neg_set["genuine"]

	print("Saving negative set to disk ... ", end="", flush=True)
	neg_set.sort_index(inplace=True)
	neg_file = args.prefix + ".neg.junctions.tab"
	neg_set.to_csv(neg_file, sep='\t')
	print("done. File saved to:", neg_file)

	if args.save_failed:
		print("Creating file containing junctions not in positive or negative set ... ", end="", flush=True)
		training = pd.concat([pos_juncs, neg_set])
		remaining = original.reset_index().merge(training, indicator=True, how='outer').set_index('index')
		others = remaining.loc[remaining['_merge'] == 'left_only']
		del others['_merge']
		other_file = args.prefix + ".others.tab"
		others.to_csv(other_file, sep='\t')
		print("done.  File saved to:", other_file)


	print()
	print("Final train set stats:")
	print(" - Positive set:", len(pos_juncs), "junctions.")
	print(" - Negative set:", len(neg_set), "junctions.")
	print(" - Others:", len(original) - len(pos_juncs) - len(neg_set), "junctions.")
	print()

def filter_one(args):
	# Load portcullis junctions into dataframe
	if args.verbose:
		print("Loading input junctions ... ", end="", flush=True)
	original = DataFrame.read_csv(args.input, sep='\t', header=0)
	if args.verbose:
		print("done.")

	df = original.copy()

	print("Input junction list contains", len(original), "junctions.")

	fieldnames = [key for key in dict(original.dtypes)]

	# Create pandas command
	pandas_cmd_in, pandas_cmd_out = json2pandas(open(args.json), fieldnames, "df")

	# print(pandas_cmd)

	# Execute the pandas command, result should be a filtered dataframe
	passed = eval(pandas_cmd_in)

	print("After filtering", len(passed), "junctions remain.")

	passed.to_csv(args.prefix + ".passed.junctions.tab", sep='\t')

	if args.save_failed:
		if args.verbose:
			print("Saving junctions failing filter ... ", end="", flush=True)
		failed = original.reset_index().merge(passed, indicator=True, how='outer').set_index('index')
		failed.loc[failed['_merge'] == 'left_only']
		del failed['_merge']
		failed.to_csv(args.prefix + ".failed.junctions.tab", sep='\t')
		if args.verbose:
			print("done.")


def main():
	parser = argparse.ArgumentParser("Script to automate CSV filtering based on a JSON configuration.")
	parser.add_argument("--json", help="Rules for filtering")
	parser.add_argument("--pos_json", nargs="*", help="File containing rules for positive set filtering.  Multiple positive rule sets allowed.  Intersection of all files taken as positive set.")
	parser.add_argument("--neg_json", nargs="*", help="File containing rules for negative set filtering.  Multiple negative rule sets allowed.  Union of all files taken as negative set")
	parser.add_argument("--genuine", help="A simple line separated list file indicating whether each junction in the input file is genuine or not 0 means not genuine, 1 means genuine.  This is used to evaulate the performance of the rule filtering.")
	parser.add_argument("--prefix", default="portcullis_filtered",
						help="The prefix to apply to all portcullis junction output files.")
	parser.add_argument("--save_layers", action='store_true', help="Whether to output the junctions at each layer")
	parser.add_argument("--save_failed", action='store_true', help="Whether to output the junctions not passing the filter")
	parser.add_argument("--verbose", "-v", action='store_true',
						help="Output additional information")

	parser.add_argument("input", help="Portcullis junction file (tab separated).")

	args = parser.parse_args()

	if args.json and (args.pos_json or args.neg_json):
		raise ValueError("Cannot use --json with --pos_json or --neg_json options ")


	if args.json:
		filter_one(args)

	elif args.pos_json and args.neg_json:
		create_training_sets(args)

	else:
		raise ValueError("Invalid configuration")



if __name__=='__main__': main()
