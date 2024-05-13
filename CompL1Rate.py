#!/usr/bin/env python
# encoding: utf-8

# File        : CompL1Rate.py
# Author      : Zhenbin Wu
# Contact     : zhenbin.wu@gmail.com
# Date        : 2018 Mar 02
#
# Description :


import pandas
import argparse
import re
import os

def ParseCSV(GT, SQ):
    df_menu = pandas.read_csv("results/L1Menu_%s_%s_emu.csv" % ("130X_dataRun3_Prompt_v3", SQ))
    df_seed = pandas.read_csv("results/L1Seed_%s_%s_emu.csv" % (GT, SQ))
    df_comp = df_seed[df_seed.L1Bit < 1000][["L1SeedName","rate0", "error_rate0"]]
    menu = df_menu.loc[df_menu.L1Bit == 9999]
    df_comp.loc[len(df_comp)] =["L1Menu", menu["rate0"].values[0], menu["error_rate0"].values[0]]
    df_comp[['rate_%s' %SQ ,'error_rate_%s'% SQ]] = df_comp[['rate0','error_rate0']].apply(pandas.to_numeric)
    df_comp.loc[:, SQ] = df_comp[['rate0','error_rate0']].apply(lambda x : '{:=8.2f}+-{:=6.2f}'.format(x[0],x[1]), axis=1)
    return df_comp.copy()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--globalTag', dest='GT', help='input globalTag')
    parser.add_argument('--sqlite1', dest='SQ1', help='input sqlite1 for Ecal laser corretion')
    parser.add_argument('--sqlite2', dest='SQ2', help='input sqlite2 for Ecal laser corretion')
    args = parser.parse_args()

    df_comp1 = ParseCSV(args.GT, args.SQ1)
    df_comp2 = ParseCSV(args.GT, args.SQ2)
    comp = df_comp1.merge(df_comp2, left_on="L1SeedName", right_on="L1SeedName", how="outer")
    comp.loc[:, "diff_value"]= (comp["rate_%s" % args.SQ2] / comp["rate_%s" % args.SQ1]) -1
    comp.loc[:, "diff_error"]=((comp["error_rate_%s" % args.SQ2]/ comp["rate_%s" % args.SQ2])**2 +
                                (comp["error_rate_%s" % args.SQ1]/ comp["rate_%s" % args.SQ1])**2).pow(1./2)
    comp.loc[:, "diff"] = comp[['diff_value','diff_error']].apply(lambda x : '{:=4.2f}+-{:=3.2f}'.format(x[0],x[1]), axis=1)
    comp = comp[["L1SeedName", args.SQ1, args.SQ2, "diff"]]
    comp.to_csv("compRate.csv")
    pandas.set_option('display.expand_frame_repr', False)
    pandas.options.display.float_format = '{:.4f}'.format
    with pandas.option_context('display.max_rows', 15, 'display.max_columns', 5):
        print(comp)

