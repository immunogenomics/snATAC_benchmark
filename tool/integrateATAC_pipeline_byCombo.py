#Imports
import argparse
from argparse import ArgumentParser
from os.path import isfile, isdir, basename, dirname, realpath
from os import listdir
import re
from submitERISslurm import SubmitJobCompact, WhichJobsAreRunning, WaitForJobs, CheckQueues
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
        if ret[0]!='' or ret[1]!='':
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

def reord_feat(ll, ft):
    if ft in ll:
        ll.remove(ft)
        ll.append(ft)
    return ll



#Current script directory
homeDir = dirname(realpath(__file__))+"/"

#Scripts
ArchR_cCRE_script = homeDir+'ArchR/ArchR_addcCREMat.R'
ArchR_peak_script = homeDir+'ArchR/ArchR_addPeakMat.R'
ArchR_basicGAS_script = homeDir+'ArchR/ArchR_addGASMat.R'
ArchR_itLSI_script = homeDir+'ArchR/ArchR_itLSI.R'
reform_emb_script = homeDir+'preprocessing/reformat_embeddings.R'
varGenes_script = homeDir+'PCA/variableGenes.R'
varPeaks_script = homeDir+'PCA/variablePeaks.R'
PCA_script = homeDir+'PCA/PCA_embedding.R'
Signac_LSI_script = homeDir+'Signac/Signac_LSI_embedding.R'
Signac_rLSI_script = homeDir+'Signac/Signac_rLSI_embedding.R'
Seurat_rCCA_script = homeDir+'Seurat/Seurat_rCCA_embedding.R'

conv_RDS_MM_script = homeDir+'preprocessing/convert_matrices.R'
conv_TXT_RDS_script = homeDir+'preprocessing/convertTXTtoRDS.R'
SnapATAC2_script = homeDir+'SnapATAC2/SnapATAC2_embedding.py'
PeakVI_ad_script = homeDir+'PeakVI/make_annData.py'
PeakVI_script = homeDir+'PeakVI/PeakVI_embedding.py'
CellSpace_getEmb_script = homeDir+'CellSpace/CellSpace_embedding_convert.R'

Harmony_script = homeDir+'Harmony/Harmony_correction.R'

UMAP_script = homeDir+'metrics/makeUMAP.R'
NNgraph_script = homeDir+'metrics/makeNNgraph.R'
NNmetrics_script = homeDir+'metrics/NNmetrics.R'
LISImetrics_script = homeDir+'metrics/LISImetrics.R'
iter_exists([ArchR_cCRE_script, ArchR_peak_script, ArchR_basicGAS_script, ArchR_itLSI_script, reform_emb_script, 
             varGenes_script, varPeaks_script, PCA_script, Signac_LSI_script, Signac_rLSI_script, Seurat_rCCA_script, 
             conv_RDS_MM_script, conv_TXT_RDS_script, SnapATAC2_script, PeakVI_ad_script, PeakVI_script, 
             CellSpace_getEmb_script, Harmony_script, UMAP_script, NNgraph_script, NNmetrics_script, LISImetrics_script])


#Input arguments
arg_obj = ArgumentParser()
overall = arg_obj.add_argument_group('Overall')
overall.add_argument('outDir',metavar = 'outDir', help='output directory')
overall.add_argument('dictionary_file',metavar = 'dictionary_file', 
                     help='dictionary file of dataset, feature, type, file')
overall.add_argument('-ds', '--dataset', action='append',type=str, 
                     help='dataset identifier(s) found in dictionary file')
overall.add_argument('-ft', '--feature', action='append',type=str, 
                     help='feature identifier(s) found in dictionary file')
overall.add_argument('-md', '--method', action='append',type=str, 
                     choices=['PCA','LSI','rLSI','itLSI','SA2','CS','pVI'], 
                     help='method identifier(s)')
overall.add_argument('-nn', '--neighbor', action='append',type=int, 
                     help='Number of nearest neighbors to use for NN metric calculation')
overall.add_argument('-LISI_skipVal', metavar='LISI_skipVal', help='values to skip for LISI calculations in style of: dataset,val;dataset,val')
overall.add_argument('-sample_col',metavar = 'sample_col', default='sample', help='Sample column in metadata; default: sample')
overall.add_argument('-phenotype_col',metavar = 'phenotype_col', default='phenotype', help='Phenotype column in metadata; default: phenotype')
overall.add_argument('-max_dim',metavar = 'max_dim', type=int, help='max dimensions; if not given, uses per-script defaults')
overall.add_argument('-binarize',action = 'store_true', help='if given, binarize input [non-GAS] matrices')
overall.add_argument('-color_file', metavar='color_file', help='file for colors')
overall.add_argument('-prefix',metavar = 'prefix', help='optional prefix for output files')
overall.add_argument('-addTime',action = 'store_true', help='if given, add /usr/bin/time in front of commands')
overall.add_argument('-dryrun',action = 'store_true', help='if given, write commands file(s), but do not run them')
overall.add_argument('-slurm', metavar='slurm', help='SLURM queue,wall time (min),memory (MB)')
overall.add_argument('-slurm_addFlags', metavar='slurm_addFlags', help='additional flags to give to SLURM, like -A account')
overall.add_argument('-seed',metavar = 'seed', type=int, help='seed; if not given, uses per-script defaults')

Signac = arg_obj.add_argument_group('Signac/Seurat')
Signac.add_argument('-Signac_NOrmDim1', action='store_true', help='if given, do not remove first dimension in Signac LSI and rLSI - often correlated with fragment count')
Signac.add_argument('-rLSI_kWeight', metavar='rLSI_kWeight', type=int, default=100, help='k.weight for Signac rLSI IntegrateEmbeddings; default 100')
Signac.add_argument('-rCCA_varFeat_method', metavar='rCCA_varFeat_method', default='vst', choices=['vst','mvp','disp'], help='FindVariableFeatures selection.method for Seurat rCCA in GAS features - options: vst, mvp, disp; default vst')

CellSpace = arg_obj.add_argument_group('CellSpace')
CellSpace.add_argument('-genome', metavar='genome', choices=['hg38','hg19'], default='hg38', help='genome version to make fasta file; default: hg38')
CellSpace.add_argument('-CS_path', metavar='CS_path', help='filepath for CellSpace command line')

args = arg_obj.parse_args()


#Verifying Input Arguments
if not isfile(args.dictionary_file):
    raise Exception("Input file must exist.")

if args.dataset is None or args.feature is None or args.method is None or args.neighbor is None:
    raise Exception('Dataset, feature, method, and NN must all have at least one input.')

#Remove duplicates - note that this does not maintain order
args.dataset = list(set(args.dataset))
args.feature = list(set(args.feature))
args.method = list(set(args.method))
args.neighbor = list(set(args.neighbor))

if args.seed is None:
    seed_str = ' '
else:
    seed_str = ' --seed %s '%(args.seed)

if args.max_dim is None:
    dim_str = ' '
else:
    dim_str = ' --max_dim %s '%(args.max_dim)

if args.binarize:
    binary_str = ' --binarize '
else:
    binary_str = ' '

if args.addTime:
    time_str = '/usr/bin/time '
else:
    time_str = ''

if 'itLSI' in args.method and len(args.feature)>1:
    args.feature = reord_feat(args.feature, 'basicGAS')
    print("\nNOTE: Because all features for the itLSI method reference the same ArchR project per dataset, only one feature can be run per project/dataset at a time. Therefore, a SLURM dependency flag will be included within dataset.\n")

useDefaults = args.slurm is None
if not useDefaults:
    queue, Wtime, mem = args.slurm.split(',')
    useDefaults = not CheckQueues(queue, Wtime, mem)
if useDefaults:
    queue='normal'
    Wtime='600'
    mem='20000'

core_str = ' '

if args.outDir[-1]!='/':
    args.outDir+='/'
mkdir(args.outDir)

clusterLogs_dir = args.outDir+'clusterLogs/'
mkdir(clusterLogs_dir)

cmdFiles_dir = args.outDir+'cmdFiles/'
mkdir(cmdFiles_dir)

#Printing Arguments
print("Argument List:")
args_print = vars(args)
for k in args_print:
    print("%s = %s"%(k,args_print[k]))
print('\n')


#Create dictionary 
dic = {}

IN = open(args.dictionary_file,'r')
for i,line in enumerate(IN):
    tabs = line.strip().split('\t')
    
    if len(tabs)!=4:
        IN.close()
        raise Exception("Need four tabs in dictionary file.")
    if i==0:
        if tabs!=['dataset','feature','type','file']:
            IN.close()
            raise Exception("Header row should be dataset, feature, type, and file")
    else:
        if not isfile(tabs[3]) and not isdir(tabs[3]):
            #archrProject is a directory
            print("SKIPPING %s: file or directory does not exist."%(i))
            continue
        
        dic[(tabs[0],tabs[1],tabs[2])]=tabs[3]
IN.close() 


#overall list
submitted_jobs = []

#file check str
bash_testFile_cmd = 'if [ ! -f %s ]; then\n\techo "%s file not found!"\nfi'

for ds in args.dataset:
    print('Dataset: '+ds)
    
    #verify ds in dictionary
    print('Verify dataset in dictionary')
    if ds not in {t[0] for t in dic.keys()}:
        print('SKIPPING %s: Dataset is not in dictionary.'%(ds))
        continue
    
    #verify meta, NNgraphHarmony, NNgraphSeurat in dictionary
    print('Verify feature-nonspecific inputs in dictionary')
    toCheck = [(ds,'NA','meta')]+[(ds,'gene','NNgraphHarmony_'+str(nn)) for nn in args.neighbor]+[(ds,'gene','NNgraphSeurat_'+str(nn)) for nn in args.neighbor]
    if not all(key in dic.keys() for key in toCheck):
        print('SKIPPING %s: Dataset does not have all necessary files (meta, NNgraphsHarmony_NN, NNgraphsSeurat_NN) in dictionary.'%(ds))
        continue
    
    meta_fl = dic[(ds,'NA','meta')]

    mkdir(args.outDir+'%s/'%(ds))
    mkdir(clusterLogs_dir+'%s/'%(ds))

    itLSI_dep_job_id = None
        
    for ft in args.feature:
        print('Feature: '+ft)
        
        #verify ft in dictionary
        print('Verify feature in dictionary')
        if ft not in {t[1] for t in dic.keys()}:
            print('SKIPPING %s: Feature is not in dictionary.'%(ft))
            continue
        
        for md in args.method:
            print('Method: '+md)
            
            #verify input files from dic
            print('Verify method inputs in dictionary')
            if md=='itLSI':
                if ft in ['archrGAS','tile']:
                    if (ds,'NA','archrProject') not in dic.keys():
                        print('SKIPPING %s %s %s: no input archrProject.'%(ds,ft,md))
                        continue
                else:
                    if not all(key in dic.keys() for key in [(ds,ft,'ext'),(ds,'NA','archrProject')]):
                        print('SKIPPING %s %s %s: no input file(s).'%(ds,ft,md))
                        continue
                    input_fl = dic[(ds,ft,'ext')]
            else:
                if (ds,ft,'matrix') not in dic.keys():
                    print('SKIPPING %s %s: no input matrix file.'%(ds,ft))
                    continue
                input_fl = dic[(ds,ft,'matrix')]

            if md=='rLSI' and ft in ['basicGAS','archrGAS']:
                print('NOTE: using rCCA instead of rLSI for GAS features')
                md = 'rCCA'

            mkdir(args.outDir+'%s/%s/'%(ds,md))
            this_outDir = args.outDir+'%s/%s/%s/'%(ds,md,ft)
            mkdir(this_outDir)
            this_prefix = '%s_%s'%(ds,ft)
            if args.prefix is not None:
                this_prefix = '%s_%s'%(args.prefix,this_prefix)

            mkdir(clusterLogs_dir+'%s/%s/'%(ds,md))
            this_clusterLogsDir = clusterLogs_dir+'%s/%s/%s/'%(ds,md,ft)
            mkdir(this_clusterLogsDir)
            this_logPrefix = '%s%s'%(this_clusterLogsDir,this_prefix)
            md_outFile = '%s_%s_embeddings.out'%(this_logPrefix,md)
            emb_fl = this_outDir+this_prefix+'_%s_embeddings.rds'%(md)

            #Make output cmd file
            print('Start command file')
            this_cmd_fl = '%s%s_%s_cmds.sh'%(cmdFiles_dir,this_prefix,md)
            OUT = open(this_cmd_fl,'w')
            
            match md:
                case "PCA":
                    print("PCA command")

                    if ft in ['basicGAS','archrGAS']:
                        varGenes_options = seed_str
                        if ft=='archrGAS':
                            varGenes_options = ' --is_norm'+varGenes_options
                        cmd = 'Rscript %s %s %s %s %s %s%s&> %s_%s_varGenes.out'%(varGenes_script, input_fl, meta_fl, 
                                                                                  args.sample_col, this_outDir, this_prefix, 
                                                                                  varGenes_options, this_logPrefix, md)
                        OUT.write(time_str+cmd+'\n\n')
                    else:
                        if args.binarize:
                            varPeaks_options = seed_str
                        else:
                            varPeaks_options = ' --NObinarize '+seed_str
                        cmd = 'Rscript %s %s %s %s %s %s%s&> %s_%s_varPeaks.out'%(varPeaks_script, input_fl, meta_fl, 
                                                                                  args.sample_col, this_outDir, this_prefix,
                                                                                  varPeaks_options, this_logPrefix, md)
                        OUT.write(time_str+cmd+'\n\n')
                    
                    cmd = 'Rscript %s %s%s_norm_scaled_transpose.rds %s %s%s&> %s'%(PCA_script, this_outDir, this_prefix, 
                                                                                    this_outDir, this_prefix, dim_str[:-1]+seed_str, 
                                                                                    md_outFile)
                    OUT.write(time_str+cmd+'\n\n')
                case "LSI":
                    if ft in ['basicGAS','archrGAS']:
                        print('SKIPPING: LSI for GAS features defaults to PCA')
                        OUT.close()
                        ret = getOneShellOutput("rm %s"%(this_cmd_fl))
                        ret = getOneShellOutput("rmdir %s"%(this_outDir))
                        ret = getOneShellOutput("rmdir %s"%(this_clusterLogsDir))
                        continue
                    else:
                        print("LSI command")
                    
                    LSI_options = binary_str[:-1]+dim_str[:-1]+seed_str
                    if args.Signac_NOrmDim1:
                        LSI_options = ' --NOrmDim1'+LSI_options
                    cmd = 'Rscript %s %s %s %s %s --sample_col %s%s&> %s'%(Signac_LSI_script, input_fl, meta_fl, 
                                                                           this_outDir, this_prefix, args.sample_col, 
                                                                           LSI_options, md_outFile)
                    OUT.write(time_str+cmd+'\n\n')
                case "rLSI":
                    if ft in ['basicGAS','archrGAS']:
                        print("SKIPPING: Should not have gotten here!")
                        OUT.close()
                        ret = getOneShellOutput("rm %s"%(this_cmd_fl))
                        ret = getOneShellOutput("rmdir %s"%(this_outDir))
                        ret = getOneShellOutput("rmdir %s"%(this_clusterLogsDir))
                        continue
                    else:
                        print('rLSI command')
                        
                        LSI_options = binary_str[:-1]+dim_str[:-1]+seed_str
                        if args.Signac_NOrmDim1:
                            LSI_options = ' --NOrmDim1'+LSI_options
                        cmd = 'Rscript %s %s %s %s %s --sample_col %s --kWeight %s%s&> %s'%(Signac_rLSI_script, input_fl, 
                                                                                            meta_fl, this_outDir, 
                                                                                            this_prefix, args.sample_col,
                                                                                            args.rLSI_kWeight, LSI_options, 
                                                                                            md_outFile)
                        OUT.write(time_str+cmd+'\n\n')
                case "rCCA":
                    if ft in ['basicGAS','archrGAS']:
                        print("rCCA command")

                        rLSI_options = dim_str[:-1]+seed_str
                        if ft=='archrGAS':
                            rLSI_options = ' --is_norm'+rLSI_options
                        cmd = 'Rscript %s %s %s %s %s --sample_col %s --varFeat_method %s%s&> %s'%(Seurat_rCCA_script, input_fl,
                                                                                                   meta_fl, this_outDir, this_prefix,
                                                                                                   args.sample_col, 
                                                                                                   args.rCCA_varFeat_method, 
                                                                                                   rLSI_options, md_outFile)
                        OUT.write(time_str+cmd+'\n\n')
                    else:
                        print("SKIPPING: Should not have gotten here!")
                        OUT.close()
                        ret = getOneShellOutput("rm %s"%(this_cmd_fl))
                        ret = getOneShellOutput("rmdir %s"%(this_outDir))
                        ret = getOneShellOutput("rmdir %s"%(this_clusterLogsDir))
                        continue
                case "itLSI":
                    print("itLSI command")
                    
                    archrProj_dir = dic[(ds,'NA','archrProject')]
                    if archrProj_dir[-1]!='/':
                        archrProj_dir+='/'
                    archrExt_dir = dirname(archrProj_dir[:-1])+'/extraction/'
                    mkdir(archrExt_dir)
                    
                    featMat_outFile = '%s_%s_featMat.out'%(this_logPrefix,md)

                    #Not binarizing the featMat; bringing in binarization at itLSI step
                    match ft:
                        case "peak":
                            cmd = 'Rscript %s %s %s %s%s&> %s'%(ArchR_peak_script, archrProj_dir, input_fl, 
                                                                archrExt_dir, seed_str, featMat_outFile)
                            OUT.write(time_str+cmd+'\n\n')
                            
                            matrix_name = 'PeakMatrix'
                            LSI_name = 'PeakLSI'
                            GAS_opt=False
                        case "cCRE":
                            cmd = 'Rscript %s %s %s cCREmatrix %s%s&> %s'%(ArchR_cCRE_script, archrProj_dir,
                                                                           input_fl, archrExt_dir, seed_str,
                                                                           featMat_outFile)
                            OUT.write(time_str+cmd+'\n\n')
                            
                            matrix_name = 'cCREmatrix'
                            LSI_name = 'cCRElsi'
                            GAS_opt=False
                        case "basicGAS":
                            cmd = 'Rscript %s %s %s %s%s&> %s'%(ArchR_basicGAS_script, archrProj_dir,
                                                                input_fl, archrExt_dir, seed_str, featMat_outFile)
                            OUT.write(time_str+cmd+'\n\n')
                            
                            matrix_name = 'GeneExpressionMatrix'
                            LSI_name = 'GeneExpLSI'
                            GAS_opt=True
                        case "tile":
                            #featMat should already exist in archrProject
                            matrix_name = 'TileMatrix'
                            LSI_name = 'TileLSI'
                            GAS_opt=False
                        case "archrGAS":
                            #featMat should already exist in archrProject
                            matrix_name = 'GeneScoreMatrix'
                            LSI_name = 'GeneScoreLSI'
                            GAS_opt=True
                    
                    if GAS_opt:
                        itLSI_options = ' --GAS_opt'+dim_str[:-1]+seed_str
                    else:
                        itLSI_options = binary_str[:-1]+dim_str[:-1]+seed_str
                    cmd = 'Rscript %s %s %s --matrix_name %s --LSI_name %s%s&> %s' %(ArchR_itLSI_script, archrProj_dir, 
                                                                                     archrExt_dir, matrix_name, LSI_name,
                                                                                     itLSI_options, md_outFile)
                    OUT.write(time_str+cmd+'\n\n')
                    core_str = ' -n 8 '
                    
                    itLSI_fl = '%s%s_%s%sIterativeLSI.rds'%(archrExt_dir,basename(archrProj_dir[:-1]),matrix_name,
                                                            '_GASopt_' if GAS_opt else '_')
                    cmd = 'Rscript %s %s %s %s %s_%s_embeddings &> %s_%s_reformatArchR.out'%(reform_emb_script, itLSI_fl, meta_fl, 
                                                                                             this_outDir, this_prefix, md, 
                                                                                             this_logPrefix, md)
                    OUT.write(time_str+cmd+'\n\n')
                case "SA2":
                    print("SA2 command")
                    
                    SA2_options = binary_str
                    if ft in ['basicGAS','archrGAS'] and args.binarize:
                        print("We are not binarizing GAS matrices.")
                        SA2_options = ' '

                    SA2_inputDir = '%sinputMM/'%(this_outDir)
                    cmd = 'Rscript %s %s %s %s --meta_cols %s,%s%s&> %s_%s_convMM.out'%(conv_RDS_MM_script, input_fl,
                                                                                        meta_fl, SA2_inputDir,
                                                                                        args.sample_col,
                                                                                        args.phenotype_col, 
                                                                                        SA2_options, this_logPrefix,md)
                    OUT.write(time_str+cmd+'\n\n')
                    
                    SA2_colFile = '%smeta_colnames.tsv'%(SA2_inputDir)
                    SA2_options = ' -colFile '+SA2_colFile+dim_str[:-1].replace('--','-')+seed_str.replace('--','-')
                    if ft=='basicGAS':
                        SA2_options = ' -GAS counts -n_feat 2000'+SA2_options
                    elif ft=='archrGAS':
                        SA2_options = ' -GAS norm -n_feat 2000'+SA2_options

                    cmd = 'python %s %s %s %s%s&> %s'%(SnapATAC2_script, SA2_inputDir, this_outDir, this_prefix, 
                                                       SA2_options, md_outFile)
                    OUT.write(time_str+cmd+'\n\n')

                    SA2_emb_txt_fl = this_outDir+this_prefix+'_SA2_embeddings.txt'
                    cmd = 'Rscript %s %s %s %s %s &> %s_SA2_convRDS.out'%(conv_TXT_RDS_script, SA2_emb_txt_fl, meta_fl,
                                                                          this_outDir, basename(SA2_emb_txt_fl).rsplit('.',1)[0], 
                                                                          this_logPrefix)
                    OUT.write(time_str+cmd+'\n\n')
                case "CS":
                    if ft in ['basicGAS','archrGAS']:
                        print('SKIPPING: CellSpace not valid for GAS features')
                        OUT.close()
                        ret = getOneShellOutput("rm %s"%(this_cmd_fl))
                        ret = getOneShellOutput("rmdir %s"%(this_outDir))
                        ret = getOneShellOutput("rmdir %s"%(this_clusterLogsDir))
                        continue
                    else:
                        if args.CS_path is None or not isdir(args.CS_path):
                            print('SKIPPING: CellSpace path not a valid file')
                            OUT.close()
                            ret = getOneShellOutput("rm %s"%(this_cmd_fl))
                            ret = getOneShellOutput("rmdir %s"%(this_outDir))
                            ret = getOneShellOutput("rmdir %s"%(this_clusterLogsDir))
                            continue
                        else:
                            print("CS command")
                    
                    if args.binarize:
                        varPeaks_options = seed_str
                    else:
                        varPeaks_options = ' --NObinarize '+seed_str
                    cmd = 'Rscript %s %s %s %s %s %s%s&> %s_%s_varPeaks.out'%(varPeaks_script, input_fl, meta_fl, 
                                                                              args.sample_col, this_outDir, this_prefix,
                                                                              varPeaks_options, this_logPrefix, md)
                    OUT.write(time_str+cmd+'\n\n')
                    
                    CS_inputDir = '%sinputMM/'%(this_outDir)
                    CS_varPeak_fl = '%s%s_varFeat.rds'%(this_outDir,this_prefix)
                    CS_options = ' --varFeat_file %s --mat_transpose --writeCells --feat_form %s%s'%(CS_varPeak_fl, 
                                                                                                    args.genome, binary_str)
                    cmd = 'Rscript %s %s %s %s%s&> %s_%s_convMM.out'%(conv_RDS_MM_script, input_fl, meta_fl, CS_inputDir, 
                                                                      CS_options, this_logPrefix, md)
                    OUT.write(time_str+cmd+'\n\n')

                    cmd = 'export PATH=%s:$PATH;%sCellSpace -output %s%s_CellSpace_model -cpMat %smatrix_transpose.mtx -peaks %sfeatures.fa -thread 5%s&> %s'%(args.CS_path, time_str, this_outDir, this_prefix, CS_inputDir, CS_inputDir, ' ' if args.max_dim is None else ' -dim %s '%(args.max_dim), md_outFile)
                    OUT.write(cmd+'\n\n')
                    core_str = ' -n 5 '

                    CS_emb_tsv_fl = '%s%s_CellSpace_model.tsv'%(this_outDir, this_prefix)
                    cmd = 'Rscript %s %s %scells.txt %s %s %s --cellspace_name %s &> %s_%s_getEmb.out'%(CellSpace_getEmb_script, 
                                                                                                        CS_emb_tsv_fl, CS_inputDir, 
                                                                                                        meta_fl, this_outDir, 
                                                                                                        this_prefix, this_prefix, 
                                                                                                        this_logPrefix, md)
                    OUT.write(time_str+cmd+'\n\n')
                case "pVI":
                    if ft in ['basicGAS','archrGAS']:
                        print('SKIPPING: PeakVI not valid for GAS features')
                        OUT.close()
                        ret = getOneShellOutput("rm %s"%(this_cmd_fl))
                        ret = getOneShellOutput("rmdir %s"%(this_outDir))
                        ret = getOneShellOutput("rmdir %s"%(this_clusterLogsDir))
                        continue
                    else:
                        print("pVI command")
                    
                    pVI_inputDir = '%sinputMM/'%(this_outDir)
                    pVI_options = ' --rmHyp --conv_meta_txt --feat_form bed'+binary_str
                    cmd = 'Rscript %s %s %s %s%s&> %s_%s_convMM.out'%(conv_RDS_MM_script, input_fl, meta_fl, pVI_inputDir,
                                                                      pVI_options, this_logPrefix,md)
                    OUT.write(time_str+cmd+'\n\n')
                    
                    conv_meta_txt = pVI_inputDir+'meta.txt' 
                    cmd = 'python %s %s %s %s %s &> %s_pVI_annData.out'%(PeakVI_ad_script, pVI_inputDir, conv_meta_txt, 
                                                                         this_outDir, this_prefix, this_logPrefix)
                    OUT.write(time_str+cmd+'\n\n')

                    pVI_ad_fl = "%s%s_AnnData.h5ad"%(this_outDir, this_prefix)
                    pVI_options = dim_str[:-1].replace('--','-')+seed_str.replace('--','-')
                    cmd = 'python %s %s %s %s -batch_key %s%s&> %s'%(PeakVI_script, pVI_ad_fl, this_outDir, this_prefix, 
                                                                     args.sample_col, pVI_options, md_outFile)
                    OUT.write(time_str+cmd+'\n\n')
                    
                    pVI_emb_txt_fl = this_outDir+this_prefix+'_pVI_embeddings.txt'
                    cmd = 'Rscript %s %s %s %s %s --addHyp &> %s_pVI_convRDS.out'%(conv_TXT_RDS_script, pVI_emb_txt_fl, 
                                                                                   meta_fl, this_outDir, 
                                                                                   basename(pVI_emb_txt_fl).rsplit('.',1)[0], 
                                                                                   this_logPrefix)
                    OUT.write(time_str+cmd+'\n\n')
                case _:
                    print("SKIPPING: method choice invalid") #shouldn't happen with input checks
                    OUT.close()
                    ret = getOneShellOutput("rm %s"%(this_cmd_fl))
                    ret = getOneShellOutput("rmdir %s"%(this_outDir))
                    ret = getOneShellOutput("rmdir %s"%(this_clusterLogsDir))
                    continue

            this_prefix = this_prefix+'_'+md
            this_logPrefix = '%s%s'%(this_clusterLogsDir,this_prefix)
            
            for hidx in range(0,2,1):
                #Harmony or not
                if hidx==1:
                    print("Harmony command")
                    cmd = 'Rscript %s %s %s %s %s %s%s&> %s_Harmony_embeddings.out'%(Harmony_script, emb_fl, meta_fl, 
                                                                                     args.sample_col, this_outDir, this_prefix, 
                                                                                     seed_str, this_logPrefix)
                    OUT.write(time_str+cmd+'\n\n')
                    this_prefix = this_prefix+'_Harmony'
                    this_logPrefix = '%s%s'%(this_clusterLogsDir,this_prefix)
                    emb_fl = this_outDir+this_prefix+'_embeddings.rds'

                
                #Post-processing/metrics
                
                print("LISI metric command")
                skip_val_lisi = []
                if args.LISI_skipVal is not None:
                    for xx in args.LISI_skipVal.split(';'):
                        yy = re.match(r'(\w+),(\w+)',xx)
                        if yy is None:
                            print("SKIPPING: LISI value %s is invalid"%(xx))
                        else:
                            if yy.groups()[0]==ds:
                                skip_val_lisi.append(yy.groups()[1])

                LISI_options = ' '
                if len(skip_val_lisi)!=0:
                    LISI_options = ' --skip_val_lisi %s '%(','.join(skip_val_lisi))
                
                cmd = 'Rscript %s %s %s %s %sLISImetrics/ %s%s&> %s_LISImetrics.out'%(LISImetrics_script, emb_fl, meta_fl, 
                                                                                      ','.join([args.sample_col,args.phenotype_col]), 
                                                                                      this_outDir, this_prefix, LISI_options,
                                                                                      this_logPrefix)
                OUT.write(time_str+cmd+'\n\n')
                LISI_outFile = '%sLISImetrics/%s_LISImetrics.rds'%(this_outDir, this_prefix)
                OUT.write(bash_testFile_cmd%(LISI_outFile, 'LISI metric')+'\n\n')

                
                print("NN graph command")
                knn_str = '--knn '+' --knn '.join([str(x) for x in args.neighbor])
                cmd = 'Rscript %s %s %sNNgraph/ %s %s &> %s_NNgraph.out'%(NNgraph_script, emb_fl, this_outDir, this_prefix,
                                                                          knn_str, this_logPrefix)
                OUT.write(time_str+cmd+'\n\n')


                for nn in args.neighbor:
                    ATAC_NNgraph_fl = this_outDir+'NNgraph/'+this_prefix+'_NN%sgraph.rds'%(str(nn))
                    for cidx in range(0,2,1):
                        if cidx==0:
                            print(str(nn)+" NN metric command w/ RNA Harmony")
                            RNA_NNgraph_fl = dic[(ds,'gene','NNgraphHarmony_'+str(nn))]
                            NNmetric_prefix = this_prefix+'_RNAPCAHarmony'+str(nn)
                            
                        elif cidx==1:
                            print(str(nn)+" NN metric command w/ RNA Seurat")
                            RNA_NNgraph_fl = dic[(ds,'gene','NNgraphSeurat_'+str(nn))]
                            NNmetric_prefix = this_prefix+'_RNASeurat'+str(nn)

                        cmd = 'Rscript %s %s %s --knn %s %sNNmetrics/ %s &> %s_NNmetrics.out'%(NNmetrics_script, RNA_NNgraph_fl, 
                                                                                               ATAC_NNgraph_fl, str(nn), this_outDir, 
                                                                                               NNmetric_prefix, 
                                                                                               this_clusterLogsDir+NNmetric_prefix)
                        OUT.write(time_str+cmd+'\n\n')
                        NN_outFile = '%sNNmetrics/%s_NNmetrics.rds'%(this_outDir, NNmetric_prefix)
                        OUT.write(bash_testFile_cmd%(NN_outFile, 'NN metric')+'\n\n')


                print("UMAP command")
                UMAP_title = ' '.join([ds,ft,md])
                if hidx==1:
                    UMAP_title += ' Harmony'
                cols_toPlot = [args.sample_col, args.phenotype_col, 'lisi_'+args.sample_col, 'lisi_'+args.phenotype_col]
                UMAP_options = ' --LISI_file %s --meta_cols \'%s\''%(LISI_outFile, ','.join(cols_toPlot))
                if args.color_file is not None:
                    UMAP_options = ' --color_file %s'%(args.color_file)+UMAP_options
                cmd = 'Rscript %s %s %s \'%s\' %sUMAP/ %s%s&> %s_UMAP.out'%(UMAP_script, emb_fl, meta_fl, UMAP_title, this_outDir,
                                                                            this_prefix, UMAP_options+seed_str, this_logPrefix)
                OUT.write(time_str+cmd+'\n\n')
                UMAP_outFile = '%sUMAP/%s_UMAP.rds'%(this_outDir,this_prefix)
                OUT.write(bash_testFile_cmd%(UMAP_outFile, 'UMAP')+'\n\n')

                            
            #close cmds file!
            print('End command file')
            OUT.close()

            
            #handling itLSI dependencies
            if md=='itLSI' and itLSI_dep_job_id is not None:
                dep_str = ' -d afterok:%s '%(itLSI_dep_job_id)
            else:
                dep_str = ' '
            
            #job string
            job_name = basename(this_cmd_fl).rsplit('.',1)[0]
            cmd = time_str+'sh %s'%(this_cmd_fl)
            queue_call = 'sbatch -p %s -t %s --mem=%sM%s--job-name=%s -o %s-%%j.out -e %s-%%j.out%s-D %s%s--wrap="%s"'%(queue, Wtime, mem, core_str, job_name, job_name, job_name, dep_str, this_clusterLogsDir, ' ' if args.slurm_addFlags is None else ' %s '%(args.slurm_addFlags), cmd)
            
            if not args.dryrun:
                print("Submit to cluster & print job ID & add to submitted list")
                job_id = SubmitJobCompact(queue_call) #prints job_string within
                submitted_jobs.append(job_id)
                print("\n"+"\t".join([ds,ft,md,job_id])+'\n')

                if md=='itLSI':
                    itLSI_dep_job_id = job_id
            else:
                print("Job that would be submitted")
                print(queue_call)

                if md=='itLSI':
                    print("\nNOTE: there is a possible itLSI dependency missing here!\n")

if not args.dryrun:
    print("Submitted job list:\n%s\n"%('\n'.join(submitted_jobs)))


