import snapatac2 as snap
import scanpy as sc
import numpy as np
import pandas as pd
from os.path import isdir, isfile
from argparse import ArgumentParser
import random

def print_time(msg):
    import time
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
arg_obj.add_argument('inDir', metavar = 'inDir', help='input directory')
arg_obj.add_argument('-colFile', metavar = 'colFile', help='if given, metadata column names')
arg_obj.add_argument('-GAS', metavar = 'GAS', choices=['counts','norm'], help='if given, gene or gene activity score matrix; choices: counts or norm')
arg_obj.add_argument('outDir', metavar = 'outDir', help='output directory')
arg_obj.add_argument('prefix', metavar = 'prefix', help='prefix for markers overlay files')
arg_obj.add_argument('-seed', metavar = 'seed', type=int, default=0, help='randomization seed; default 0')
arg_obj.add_argument('-n_feat', metavar = 'n_feat', type=int, default=50000, help='number of features; default 50000')
arg_obj.add_argument('-max_dim', metavar = 'max_dim', type=int, default=30, help='Number of spectral dimensions to generate; default 30')
args = arg_obj.parse_args()

if not isdir(args.inDir):
    raise Exception('Input directory does not exist. Please try again.')

if not isfile(args.inDir+'matrix.mtx') and not isfile(args.inDir+'matrix.mtx.gz'):
    raise Exception('Matrix file must exist.')
if not isfile(args.inDir+'features.tsv') and not isfile(args.inDir+'features.tsv.gz'):
    raise Exception('Feature file must exist.')
if not isfile(args.inDir+'barcodes.tsv') and not isfile(args.inDir+'barcodes.tsv.gz'):
    raise Exception('Barcode file must exist.')

if args.colFile is not None and not isfile(args.colFile):
    raise Exception('If given, meta column name file must exist.')

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
random.seed(args.seed)

print_time("Load Data")
data = snap.read_10x_mtx(args.inDir)
print(data)

if args.colFile is not None:
    colNames = np.loadtxt(args.colFile, delimiter="\t",dtype='str')
    if colNames.shape[0]!=data.obs.shape[1]:
        raise Exception('meta column names must be the same length as data.obs')
    data.obs.columns = colNames
print(data.obs)
print(data.var)

assert data.n_obs == np.unique(data.obs_names).size
assert data.var_names.size == np.unique(data.var_names).size

if args.GAS=='counts':
    print_time("HVG, Norm, Log")
    sc.pp.highly_variable_genes(data, flavor='seurat_v3', n_top_genes=args.n_feat)
    data = data[:, data.var.highly_variable]
    sc.pp.normalize_total(data)
    sc.pp.log1p(data)

    print_time("Spectral Embedding")
    snap.tl.spectral(data,n_comps=args.max_dim,random_state=args.seed, features=None)
elif args.GAS=='norm':
    print_time("Log, HVG")
    sc.pp.log1p(data)
    sc.pp.highly_variable_genes(data, flavor='seurat', n_top_genes=args.n_feat)
    data = data[:, data.var.highly_variable]

    print_time("Spectral Embedding")
    snap.tl.spectral(data,n_comps=args.max_dim,random_state=args.seed, features=None)
else:
    print_time("Filter features")
    snap.pp.select_features(data, n_features=args.n_feat)
    
    print_time("Spectral Embedding")
    snap.tl.spectral(data,n_comps=args.max_dim,random_state=args.seed)


spectral_mat = pd.DataFrame(data.obsm["X_spectral"],
                            index=data.obs_names,
                            columns=['spectral%s'%(i) for i in range(1,data.obsm["X_spectral"].shape[1]+1)])
spectral_mat.to_csv('%s_SA2_embeddings.txt'%(outPrefix), index=True, header=True, sep='\t')

print_time("Save")
data.write(outPrefix+'_SA2_object.h5ad')

print_time("Done.")


