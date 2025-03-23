import logging
from subprocess import call
from os.path import isfile, basename, dirname


def SubmitJobCompact(job_string):
    print(job_string)
    
    from subprocess import Popen, PIPE
    import re
    job_process = Popen(job_string, shell = True, stdout = PIPE, stderr = PIPE, cwd = None)
    out_stream, err_stream = job_process.communicate()
    job_process.wait()
    out_stream=out_stream.decode('UTF-8')
    err_stream=err_stream.decode('UTF-8')
    
    logging.info(out_stream)
    print(str(out_stream).strip())
    
    if len(err_stream) > 0:
        raise Exception('Error submitting job: %s \nException: %s' %(job_string, err_stream))
    
    job_id = re.search(r'Job <([0-9]+)> is submitted to queue', out_stream)
    
    if not job_id is None:
        return str(job_id.group(1))
    else:
        raise Exception('Job ID could not be retrieved from LSF system. This is likely an error.')


def SubmitJob(job_string, log_file = 'job_log.txt', memory = 6, cwd = None, num_cores = 1, depends = None, wd = None, extra_flags = ''):
    """Submits a job to the LSF queueing system. The string should be
    exactly what would be typed at the Unix shell in the current
    directory.
    
    Returns the job ID given by the LSF system.
    """
    
    from subprocess import Popen, PIPE
    import re
    from string import replace

    
    job_string = replace(job_string, '$','\$')
    job_string = replace(job_string, '"','\"')

    
    if depends is not None:
        if type(depends) is str or type(depends) is int:
            depends = [str(depends)]

        job_depend_str=' && '.join(['done(%s)' %(_) for _ in depends])
        
        extra_flags += ' -w "%s" ' %job_depend_str
    else:
        extra_flags= ' '
    if wd is not None:
        working_dir=' -cwd %s '%wd
    else:
        working_dir=' '
    
    if name=='':
        name = job_string.split(" ")[1].split("/")[-1]
    out_file = name+"-%J.out"
    err_file = name+"-%J.err"
    
    queue_call = 'bsub -q normal -M %s -J %s -o %s -e %s%s%s"%s"'%(memory*1000, name, out_file, err_file, working_dir, extra_flags, 
                                                                   job_string)
        
    print(queue_call)
    

    job_process = Popen(queue_call, shell = True, stdout = PIPE, stderr = PIPE, cwd = None)
    out_stream, err_stream = job_process.communicate()
    job_process.wait()
    out_stream=out_stream.decode('UTF-8')
    err_stream=err_stream.decode('UTF-8')

    logging.info(out_stream)
    print(str(out_stream).strip())
    
    if len(err_stream) > 0:
        raise Exception('Error submitting job: %s \nException: %s' %(job_string, err_stream))
    
    job_id = re.search(r'Job <([0-9]+)> is submitted to queue', out_stream)
    
    if not job_id is None:
        return str(job_id.group(1))
    else:
        raise Exception('Job ID could not be retrieved from LSF system. This is likely an error.')
    

def SubmitJobsAndWait(list_of_jobs, num_cores = 1, print_jobs = True, wd = None, memory = 6, **kwargs):
    """Submits a list of string specifiying jobs. These should be written as they'd be typed
    into the shell.
    kwargs go straight to the 'SubmitJob" function."""    
    
    from time import sleep

    if len(list_of_jobs) == 0:
        return 0

    submitted_jobs = set()
    for job_string in list_of_jobs:
        print('Submitting to ERISone: \n %s' %job_string)
        submitted_jobs.add(SubmitJob(job_string, num_cores = num_cores, memory = memory, wd = wd, **kwargs))
    
    sleep(10)
    while True:
        
        running_jobs = WhichJobsAreRunning()
        
        #This intersects the set of running and submitted jobs
        #If the set is empty, this means that all the jobs have finished
        if len(running_jobs & submitted_jobs) == 0:
            break
        else:
            sleep(10)


def WhichJobsAreRunning():
    
    from subprocess import Popen, PIPE
    import re
    
    
    job_process = Popen('bjobs', shell = True, stdout = PIPE, stderr = PIPE)
    
    out_stream, err_stream = job_process.communicate()
    out_stream=out_stream.decode('UTF-8')
    err_stream=err_stream.decode('UTF-8')

    
    job_process.wait()
    
    job_ids = []
    
    
    for l in out_stream.split('\n'):
        parse = re.search(r"([0-9]+)\s+",l)
        if parse is not None:
            job_ids.append(parse.group(1))
    
    return set(job_ids)


def WaitForJobs(job_id_list):
    """Similar to SubmitJobsAndWait, but gets a list of job numbers for jobs that are assumed to have
    already been submitted."""
    
    from time import sleep
    
    if job_id_list is None:
        return 0

    submitted_jobs = set(job_id_list)
    print("Waiting for jobs:\n%s\n"%("\n".join(submitted_jobs)))
    
    sleep(10)
    while True:

        running_jobs = WhichJobsAreRunning()
        
        #This intersects the set of running and submitted jobs
        #If the set is empty, this means that all the jobs have finished
        if len(running_jobs & submitted_jobs) == 0:
            break
        else:
            sleep(10)


def getOneShellOutput(cmd):
    from subprocess import call, PIPE, Popen
    job_process = Popen(cmd, shell=True, stdout = PIPE, stderr = PIPE)
    out_stream, err_stream = job_process.communicate()
    return [out_stream.strip().decode("UTF-8"),err_stream.decode("UTF-8")]


def CheckQueues(queue, Wtime, mem, avail_queues = {'normal':[21600,32000],'bigmem':[7200,400000]}):
    #NOTE: they do not list a max mem for bigmem

    allgood = True
    
    if queue not in avail_queues.keys():
        print("Queue doesn't exist on ERIStwo")
        allgood = False
    if int(Wtime)>avail_queues.get(queue)[0]:
        print("Time limit must be less than or equal to %s for queue %s"%(Wtime,avail_queues.get(queue)[0]))
        allgood = False
    if int(mem)>avail_queues.get(queue)[1]:
        print("Memory limit must be less than or equal to %s for queue %s"%(mem,avail_queues.get(queue)[1]))
        allgood = False

    return allgood

    
