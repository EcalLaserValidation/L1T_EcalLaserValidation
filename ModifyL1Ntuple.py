#!/usr/bin/env python
# encoding: utf-8

# File        : change.py
# Author      : Zhenbin Wu
# Contact     : zhenbin.wu@gmail.com
# Date        : 2018 Mar 01
#
# Description :

import argparse
import re
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--globalTag', dest='GT', help='input globalTag')
    parser.add_argument('--sqlite', dest='SQ', help='input sqlite for Ecal laser corretion')
    args = parser.parse_args()
    inputfile = open("l1Ntuple_%s.py" % args.GT, 'r')
    outfile = open("l1Ntuple_%s_%s.py" % (args.GT, args.SQ), 'w')
    for line  in inputfile.readlines():
        if re.match("^process\s=\scms.Process\('RAW2DIGI'.*\)$", line) is not None:
            outfile.write(line)
            outfile.write('\n')
            outfile.write("import FWCore.ParameterSet.VarParsing as VarParsing\n")
            outfile.write("# setup 'analysis'  options\n")
            outfile.write("options = VarParsing.VarParsing ('analysis')\n")
            outfile.write("# get and parse the command line arguments\n")
            outfile.write("options.parseArguments()\n")
            outfile.write('\n')
        elif re.match("\s*fileNames\s=\scms.untracked.vstring\('inputFiles'\),", line) is not None:
            outfile.write("\tfileNames = cms.untracked.vstring(options.inputFiles),\n")
        elif re.match("^process.GlobalTag\s=\sGlobalTag.*%s.*" % args.GT, line) is not None:
            outfile.write(line)
            outfile.write('process.GlobalTag.toGet = cms.VPSet(\n')
            outfile.write('   cms.PSet(record = cms.string("EcalTPGLinearizationConstRcd"),\n')
            outfile.write('            tag = cms.string("EcalTPGLinearizationConst_IOV_%s_beginning_at_1"),\n' % args.SQ)
            outfile.write('            connect =cms.string("sqlite_file:%s/EcalTPG_%s_moved_to_1.db"),\n' % (os.getcwd(), args.SQ))
            outfile.write('            ),\n')
            outfile.write('   cms.PSet(record = cms.string("EcalTPGPedestalsRcd"),\n')
            outfile.write('            tag = cms.string("EcalTPGPedestals_IOV_%s_beginning_at_1"),\n' % args.SQ)
            outfile.write('            connect =cms.string("sqlite_file:%s/EcalTPG_%s_moved_to_1.db"),\n' % (os.getcwd(), args.SQ))
            outfile.write('            ),\n')
            outfile.write(')\n')
            outfile.write('process.MessageLogger.cerr.FwkReport.reportEvery = 10000\n')
        else:
            outfile.write(line)

    outfile.write("process.TFileService.fileName= cms.string(options.outputFile)\n")
