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
    print(str(err_stream).strip())
    
    job_id = re.search(r'Submitted batch job ([0-9]+)', out_stream)
    
    if not job_id is None:
        return str(job_id.group(1))
    else:
        raise Exception('Job ID could not be retrieved from SLURM system. This is likely an error.')


def WhichJobsAreRunning():
    
    from subprocess import Popen, PIPE
    import re
    
    
    job_process = Popen('sacct -s R', shell = True, stdout = PIPE, stderr = PIPE)
    
    out_stream, err_stream = job_process.communicate()
    out_stream=out_stream.decode('UTF-8')
    err_stream=err_stream.decode('UTF-8')

    
    job_process.wait()
    
    job_ids = []
    
    
    for l in out_stream.split('\n'):
        parse = re.search(r"^([0-9]+)\s+",l)
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


def CheckQueues(queue, Wtime, mem, 
                avail_queues = {'short':[180,360000],'normal':[7200,360000],'bigmem':[11520,1000000],'long':[23040,360000]}):

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

    
