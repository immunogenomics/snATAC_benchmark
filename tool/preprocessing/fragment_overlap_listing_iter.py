import argparse
from argparse import ArgumentParser
from os.path import isfile, isdir, basename, dirname, realpath
from os import listdir
import re
import time


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

def makeListFile(toDo_file, toDo_arr):
    
    with open(toDo_file,'w') as f:
        for fl in toDo_arr:
            if not isfile(fl):
                raise Exception("ERROR: output file doesn't exist: %s"%(fl))
            
            f.write(fl+'\n')


def exists(f):
    if not isfile(f):
        raise Exception("File doesn't exist: %s"%(f))

def iter_exists(l):
    for x in l:
        exists(x)


"""Current script directory"""
homeDir = dirname(realpath(__file__))+"/"

files_required = []
"""Scripts"""
fragOver_script = homeDir+'fragment_overlap.R'
concat_matrices_script = homeDir+'concat_matrices.R'
iter_exists([fragOver_script, concat_matrices_script])


"""Input Arguments"""
arg_obj = ArgumentParser()
arg_obj.add_argument('listing',metavar = 'listing', help='listing file - 2 tabbed columns: donorID fragments')
arg_obj.add_argument('outDir', metavar = 'outDir', help='output directory')
arg_obj.add_argument('prefix', metavar='prefix', help='Prefix for output files')
arg_obj.add_argument('byType', metavar = 'byType', help='feature type for output name (e.g., peak, basicGAS)')
arg_obj.add_argument('byFile', metavar = 'byFile', help='feature file')
arg_obj.add_argument('-feat_colNum', metavar = 'feat_colNum', type=int, help='if given, feature file column number to use as feature name; if not given, default to col1:col2-col3')
arg_obj.add_argument('-concat', metavar = 'concat', help='output directory for concatenate matrices over donors; waits for jobs to finish.')
arg_obj.add_argument('-NOrm1', action='store_true', help='if given, do NOT remove the "-1" at the end of the cell barcode; note it will verify that all cell barcodes have it before removing any!')
arg_obj.add_argument('-addTime',action = 'store_true', help='if given, add /usr/bin/time in front of commands')
args = arg_obj.parse_args()



"""Verifying Input Arguments
A lot of this is done in the ArgumentParser"""
if not isfile(args.listing) or not isfile(args.byFile):
    raise Exception("Input files must exist.")
    
if args.addTime:
    time_str = '/usr/bin/time '
else:
    time_str = ''

if args.concat is not None:
    if args.concat[-1]!='/':
        args.concat+='/'
    mkdir(args.concat)
    
if args.outDir[-1]!='/':
    args.outDir+='/'
mkdir(args.outDir)

clusterLogs = args.outDir+'clusterLogs/'
mkdir(clusterLogs)


"""Printing Arguments"""
print("Argument List:")
args_print = vars(args)
for k in args_print:
    print("%s = %s"%(k,args_print[k]))
print('\n')



output_cells_files = []
output_donors_files = []
IN = open(args.listing,'r')
for i,line in enumerate(IN):
    tabs = line.strip().split('\t')
    
    donorID = tabs[0]
    fragments_file = tabs[1]
    
    if not isfile(fragments_file):
        print("SKIPPING: fragment file does not exist! %s"%(fragments_file))
        continue
    
    outPrefix = args.outDir+args.prefix+'_'+donorID+'_'+args.byType
    
    output_cell_file = outPrefix+'Xcells.rds'
    output_cells_files.append(output_cell_file)
    output_donors_files.append(outPrefix+'Xdonor.rds')
    
    cmd = time_str+"Rscript %s %s %s %s %s --donor_name %s"%(fragOver_script,args.byFile,fragments_file,args.outDir,
                                                             args.prefix+'_'+donorID+'_'+args.byType,donorID)
    if args.feat_colNum is not None:
        cmd += ' --feat_colNum %s'%(args.feat_colNum)
    if args.NOrm1 is True:
        cmd += ' --NOrm1'
    
    #print(cmd)
    ret = getOneShellOutput(cmd)
    print("OUTPUT:\n%s\n"%(ret[0]))
    print("ERROR:\n%s\n\n"%(ret[1]))
    
IN.close()


if args.concat is not None:
    
    concat_prefix = args.concat+args.prefix+"_" + args.byType
    
    cells_file_list = concat_prefix+'Xcells_list.txt'
    makeListFile(cells_file_list, output_cells_files)
    donors_file_list = concat_prefix+'Xdonor_list.txt'
    makeListFile(donors_file_list, output_donors_files)
    
    cmd = time_str+"Rscript %s %s row %s %s"%(concat_matrices_script,cells_file_list,args.concat,
                                              args.prefix+'_'+args.byType+'Xcells')
    #print(cmd)
    ret = getOneShellOutput(cmd)
    print("OUTPUT:\n%s\n"%(ret[0]))
    print("ERROR:\n%s\n\n"%(ret[1]))

    cmd = time_str+"Rscript %s %s row %s %s"%(concat_matrices_script,donors_file_list,args.concat,
                                              args.prefix+'_'+args.byType+'Xdonors')
    #print(cmd)
    ret = getOneShellOutput(cmd)
    print("OUTPUT:\n%s\n"%(ret[0]))
    print("ERROR:\n%s\n\n"%(ret[1]))
    
    if not isfile(concat_prefix+"Xcells.rds") or not isfile(concat_prefix+"Xdonors.rds"):
        raise Exception("ERROR: concatenation files don't exist.")
    

