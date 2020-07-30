#!/usr/bin/env python
import warnings
import MYUTILS
import sys
import argparse
parser = argparse.ArgumentParser(description='Provided with a set of barcodes, find all those that collide within N edit distance.')
parser.add_argument('-i',dest='inFP',	metavar='<inFile>',help='Input file of barcodes (barcode\tname)', required=True);
parser.add_argument('-o',dest='outFP', metavar='<outFile>',help='Where to output results [default=stdout]', required=False);
parser.add_argument('-d',dest='distance', metavar='<distance>',help='Max edit distance [default=1]', default=1, required=False);
parser.add_argument('-l',dest='logFP', metavar='<logFile>',help='Where to output errors/warnings [default=stderr]', required=False);
parser.add_argument('-v',dest='verbose', action='count',help='Verbose output?', required=False, default=0);

args = parser.parse_args();
args.distance = int(args.distance)

inFile=MYUTILS.smartGZOpen(args.inFP,'r');


if (args.logFP is not None):
	logFile=MYUTILS.smartGZOpen(args.logFP,'w');
	sys.stderr=logFile;

if (args.outFP is None):
	outFile= sys.stdout;
else:
	if args.verbose>0: warnings.warn("Outputting to file "+args.outFP);
	outFile = MYUTILS.smartGZOpen(args.outFP,'w');

bases = ["A","T","G","C"];

def addWithinDistanceRecur(edHash, curOligo, nEdits, curName): h
	global bases;
	if nEdits>0:
		for i in range(0, len(curOligo)):
			for b in bases:
				oligoI = curOligo[0:i]+b+curOligo[i+1:len(curOligo)];
				addWithinDistanceRecur(edHash,oligoI,nEdits-1, curName);
	else:
		if curOligo not in edHash:
			edHash[curOligo] = {};
		edHash[curOligo][curName]=curName;


withinEditDistance = {};
bc2names = {};
#raise Exception("Reached bad state=%d for '%s.%d' '%s' at line '%s'" %(state,mid,ver,tfid,line));
for line in inFile:
	if line is None or line == "" or line[0]=="#": continue
	data=line.rstrip().split("\t");
	bc = data[0];
	name = data[1];
	addWithinDistanceRecur(withinEditDistance, bc, args.distance, bc);
	bc2names[bc]=name;

inFile.close();

notedCollisions = {};
for k in withinEditDistance:
	if len(withinEditDistance[k])>1:
		message = "Collision found for %s between: "%(k);
		curSet = ""
		for bc in sorted(withinEditDistance[k], key = withinEditDistance[k].get):
			message = "%s %s (%s)"%(message, bc, bc2names[bc]);
			curSet = curSet + bc;
		if curSet not in notedCollisions:
			outFile.write(message+"\n");
			notedCollisions[curSet]=1;


outFile.close();
if (args.logFP is not None):
	logFile.close();
