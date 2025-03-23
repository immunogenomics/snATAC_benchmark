import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import scvi

from os.path import isdir, isfile
import time
import argparse
from argparse import ArgumentParser

def print_time(msg):
    print("\n%s - %s\n"%(time.strftime("%c"),msg))

def getOneShellOutput(cmd):
    print(cmd)
    from subprocess import call, PIPE, Popen
    job_process = Popen(cmd, shell=True, stdout = PIPE, stderr = PIPE)
    out_stream, err_stream = job_process.communicate()
    return [out_stream.strip().decode("UTF-8"),err_stream.decode("UTF-8")]

def mkdir(x):
    if not isdir(x):
        cmd = "mkdir %s"%(x)
        ret = getOneShellOutput(cmd)
        print(ret)

        if not isdir(x):
            return False
        else:
            return True
    else:
        return True

arg_obj = ArgumentParser()
arg_obj.add_argument('inDir', metavar = 'inDir', help='input directory in 10x file format')
arg_obj.add_argument('meta_file', metavar = 'meta_file', help='meta file')
arg_obj.add_argument('outDir', metavar = 'outDir', help='output directory')
arg_obj.add_argument('prefix', metavar = 'prefix', help='prefix for markers overlay files')
args = arg_obj.parse_args()

if not isdir(args.inDir) or not isfile(args.meta_file):
    raise Exception('Input directory does not exist. Please try again.')

if args.outDir[-1]!='/':
    args.outDir+='/'
mkdir(args.outDir)
outPrefix = args.outDir+args.prefix

"""Printing Arguments"""
print("Argument List:")
args_print = vars(args)
for k in args_print:
    print("%s = %s"%(k,args_print[k]))


print_time('Load Files')
meta = pd.read_csv(args.meta_file,sep='\t',header=0,index_col=0)
print(meta.head())

adata = scvi.data.read_10x_atac(args.inDir)
print(adata)
print(adata.obs.head())


print_time('Add metadata')
if not meta.index.equals(adata.obs.index):
    raise Exception('meta and AnnData indices do not match.')
    
for xx in meta.dtypes.index:
    if meta.dtypes[xx]=='object':
        adata.obs[xx] = pd.Categorical(meta[xx])
    else:
        adata.obs[xx] = meta[xx]
print(adata.obs.head())


print_time('Make variables unique')
adata.obs_names_make_unique()
adata.var_names_make_unique()


print_time('Save AnnData')
adata.write("%s_AnnData.h5ad"%(outPrefix))


print_time('Done.')
