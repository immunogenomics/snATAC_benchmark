#Imports
import argparse
from argparse import ArgumentParser
from os.path import isfile, isdir, basename, dirname, realpath
from os import listdir
import re
from submitERIStwo import SubmitJobCompact, SubmitJobsAndWait, WhichJobsAreRunning, WaitForJobs, CheckQueues
import time

#Functions
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

def print_time(msg):
    print("\n%s - %s\n"%(time.strftime("%c"),msg))

def exists(f):
    if not isfile(f):
        raise Exception("File doesn't exist: %s"%(f))

def iter_exists(l):
    for x in l:
        exists(x)


#Current script directory
homeDir = dirname(realpath(__file__))+"/"

#Scripts
varGenes_script = homeDir+'PCA/variableGenes.R'
PCA_script = homeDir+'PCA/PCA_embedding.R'
Harmony_script = homeDir+'Harmony/Harmony_correction.R'
NNgraph_script = homeDir+'metrics/makeNNgraph.R'
UMAP_script = homeDir+'metrics/makeUMAP.R'
iter_exists([varGenes_script,PCA_script,Harmony_script,NNgraph_script,UMAP_script])


#Input arguments
arg_obj = ArgumentParser(description='Make RNA PCA + Harmony NN graph(s) + UMAP')
arg_obj.add_argument('gxc_file',metavar = 'gxc_file', help='Genes x cells file - raw counts')
arg_obj.add_argument('meta_file',metavar = 'meta_file', help='Cell metadata file - rownames as cellnames')
arg_obj.add_argument('--sample_col', metavar='sample_col', default='sample', help='sample column in metadata file; default sample')
arg_obj.add_argument('--norm_scale_factor', metavar='norm_scale_factor', default=10000, type=int, help='scaling factor for normalization; default=10000')
arg_obj.add_argument('--exp_min', metavar='exp_min', default=0.1, type=float, help='expression minimum; default 0.1')
arg_obj.add_argument('--numGenes', metavar='numGenes', default=500, type=int, help='number of variable genes; default 500')
arg_obj.add_argument('--max_dim', metavar='max_dim', default = 30, type=int, help='number of PCs to calculate; default 30')
arg_obj.add_argument('-nn', '--neighbor', action='append',type=int, 
                     help='Number of nearest neighbors to use for NN metric calculation')
arg_obj.add_argument('--color_file', metavar='color_file', help='file for colors')
arg_obj.add_argument('--clusterLogDir',metavar = 'clusterLogDir', help='clusterLogDir; if not given, outDir/clusterLogs/')
arg_obj.add_argument('outDir',metavar = 'outDir', help='output directory')
arg_obj.add_argument('prefix',metavar = 'prefix', help='prefix for output files')
arg_obj.add_argument('-seed', metavar = 'seed', type=int, help='seed; if not given, uses per-script defaults')
arg_obj.add_argument('-dryrun',action = 'store_true', help='if given, write commands file(s), but do not run them')
arg_obj.add_argument('-addTime',action = 'store_true', help='if given, add /usr/bin/time in front of commands')
args = arg_obj.parse_args()


#Verifying Input Arguments
if not isfile(args.gxc_file) or not isfile(args.meta_file):
    raise Exception("Input files must exist.")

if args.outDir[-1]!='/':
    args.outDir+='/'
mkdir(args.outDir)

if args.clusterLogDir is None:
    clusterLogs_dir = args.outDir+'clusterLogs/'
    mkdir(clusterLogs_dir)
else:
    clusterLogs_dir = args.clusterLogDir

if args.addTime:
    time_str = '/usr/bin/time '
else:
    time_str = ''

#Printing Arguments
print("Argument List:")
args_print = vars(args)
for k in args_print:
    print("%s = %s"%(k,args_print[k]))
print('\n')


submitted_jobs = []

print_time('Variable Genes')
cmd = 'Rscript %s %s %s %s %s %s --numGenes %s --norm_scale_factor %s --exp_min %s%s&> %s%s_varGenes.out'%(varGenes_script, args.gxc_file, args.meta_file, args.sample_col, args.outDir, args.prefix, args.numGenes, args.norm_scale_factor, args.exp_min, ' ' if args.seed is None else ' --seed %s '%(args.seed), clusterLogs_dir, args.prefix)
cmd = time_str + cmd
if args.dryrun:
    print(cmd)
else:
    ret = getOneShellOutput(cmd)
    print("\nOUTPUT:\n%s\n"%(ret[0]))
    print("ERROR:\n%s\n"%(ret[1]))

print_time('PCA')
cmd = 'Rscript %s %s%s_norm_scaled_transpose.rds %s %s --max_dim %s%s&> %s%s_PCA_embeddings.out'%(PCA_script, args.outDir, args.prefix, args.outDir, args.prefix, str(args.max_dim), ' ' if args.seed is None else ' --seed %s '%(args.seed), clusterLogs_dir, args.prefix)
cmd = time_str + cmd
if args.dryrun:
    print(cmd)
else:
    ret = getOneShellOutput(cmd)
    print("\nOUTPUT:\n%s\n"%(ret[0]))
    print("ERROR:\n%s\n"%(ret[1]))

print_time('Harmony')
cmd = 'Rscript %s %s%s_PCA_embeddings.rds %s %s %s %s_PCA%s&> %s%s_PCA_Harmony_embeddings.out'%(Harmony_script, args.outDir, args.prefix, args.meta_file, args.sample_col, args.outDir, args.prefix,  ' ' if args.seed is None else ' --seed %s '%(args.seed), clusterLogs_dir, args.prefix)
cmd = time_str + cmd
if args.dryrun:
    print(cmd)
else:
    ret = getOneShellOutput(cmd)
    print("\nOUTPUT:\n%s\n"%(ret[0]))
    print("ERROR:\n%s\n"%(ret[1]))

print_time('NN graph')
cmd = 'Rscript %s %s%s_PCA_Harmony_embeddings.rds %sNNgraph/ %s_PCA_Harmony %s%s&> %s%s_PCA_Harmony_NNgraph.out'%(NNgraph_script, args.outDir, args.prefix, args.outDir, args.prefix, '--knn '+' --knn '.join([str(x) for x in args.neighbor]),  ' ' if args.seed is None else ' --seed %s '%(args.seed), clusterLogs_dir, args.prefix)
cmd = time_str + cmd
if args.dryrun:
    print(cmd)
else:
    ret = getOneShellOutput(cmd)
    print("\nOUTPUT:\n%s\n"%(ret[0]))
    print("ERROR:\n%s\n"%(ret[1]))

print_time('UMAP')
cmd = 'Rscript %s %s%s_PCA_Harmony_embeddings.rds %s \'RNA PCA + Harmony\' %sUMAP/ %s_PCA_Harmony --meta_cols %s --color_file %s%s&> %s%s_PCA_Harmony_UMAP.out'%(UMAP_script, args.outDir, args.prefix, args.meta_file, args.outDir, args.prefix, args.sample_col,  args.color_file,  ' ' if args.seed is None else ' --seed %s '%(args.seed), clusterLogs_dir, args.prefix)
cmd = time_str + cmd
if args.dryrun:
    print(cmd)
else:
    ret = getOneShellOutput(cmd)
    print("\nOUTPUT:\n%s\n"%(ret[0]))
    print("ERROR:\n%s\n"%(ret[1]))

print_time('Done.')

