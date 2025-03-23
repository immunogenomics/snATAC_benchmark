import sys
import scvi
import anndata
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

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
arg_obj.add_argument('adataFile', metavar = 'adataFile', help='input adata file')
arg_obj.add_argument('outDir', metavar = 'outDir', help='output directory')
arg_obj.add_argument('prefix', metavar = 'prefix', help='prefix for markers overlay files')
arg_obj.add_argument('-seed', metavar = 'seed', type=int, default=420, help='randomization seed; default 420')
arg_obj.add_argument('-filtPerc', metavar = 'filtPerc', type=float, default=0.05, help='filter percentage for features in cells; default 0.05')
arg_obj.add_argument('-batch_key', metavar = 'batch_key', help='batch_key for adata')
arg_obj.add_argument('-cat_key', metavar = 'cat_key', action='append', help='categorical_covariate_keys for adata')
arg_obj.add_argument('-max_dim', metavar = 'max_dim', type=int, help='Number of latent dimensions to generate; default None')
args = arg_obj.parse_args()

if not isfile(args.adataFile):
    raise Exception('Input file does not exist. Please try again.')

if args.outDir[-1]!='/':
    args.outDir+='/'
mkdir(args.outDir)
outPrefix = args.outDir+args.prefix

"""Printing Arguments"""
print("Argument List:")
args_print = vars(args)
for k in args_print:
    print("%s = %s"%(k,args_print[k]))

print_time("Set seed: %s"%(args.seed))
scvi.settings.seed = args.seed

print_time("Load Data")
adata_pvi = anndata.read_h5ad(args.adataFile)
adata_pvi.var_names_make_unique()
print(adata_pvi)

print_time("Filter features to keep those in at least X% of cells")
print(adata_pvi.shape)
sc.pp.filter_genes(adata_pvi, min_cells=int(adata_pvi.shape[0] * args.filtPerc))
print(adata_pvi.shape)

print_time("Setup anndata")
scvi.model.PEAKVI.setup_anndata(adata_pvi, batch_key=args.batch_key, categorical_covariate_keys=args.cat_key)

print_time("Train model")
pvi = scvi.model.PEAKVI(adata_pvi,n_latent=args.max_dim)
pvi.train()

print_time("Save model")
pvi.save("%strained_peakvi"%(args.outDir),save_anndata=True)

print_time("Extract latent space")
adata_pvi.obsm["PeakVI_latent"] = pvi.get_latent_representation()
adata_pvi.write("%s_pVI_object.h5ad"%(outPrefix))

print_time("Save latent space")
latent_mat = pd.DataFrame(adata_pvi.obsm["PeakVI_latent"],
                          index=adata_pvi.obs.index,
                          columns=['latent%s'%(i) for i in range(1,adata_pvi.obsm["PeakVI_latent"].shape[1]+1)])
latent_mat.to_csv('%s_pVI_embeddings.txt'%(outPrefix), index=True, header=True, sep='\t')

print_time("Done.")
