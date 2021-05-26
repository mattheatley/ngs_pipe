import os, sys, subprocess, argparse, time
from core import CAPTURE

parser = argparse.ArgumentParser(description='Sub Buddy', prog=os.path.basename(__file__), usage='%(prog)s [options]', epilog='to assist re-submitting SLURM tasks')
parser.add_argument('-u', metavar='<name>', type=str, help='specify user name (default: find current user)')
parser.add_argument('-i', metavar='</path/to/id_file>', type=str, help='specify path to slurm id file')
user, id_file = vars(parser.parse_args()).values() # define user inputs

subprocess.call('echo BUDDY STARTED `date`', shell=True) # log script start

while True: # continuously check tasks

    HEADERS, *running_tasks = [ line.strip().split('|') for line in CAPTURE(f'squeue -u {user} --format "%i|%j|%R" | grep -v "buddy"').split('\n') ] # count tasks excluding subbuddy & pals

    if len(running_tasks) == 0: break # exit if no running tasks
    else: # examine any running tasks

        *_, NAME, NODELIST = HEADERS # specify headers
        running_tasks = { job_id: {NAME: job_name, NODELIST: nodelist_reason} for job_id, job_name, nodelist_reason in running_tasks } # organise task info

        with open(f'{id_file}', 'a+') as id_update: # append id file

            id_update.seek(0) # reset file position
            submitted_tasks = { job_id: script for job_id, script, *_ in [ line.strip().split('\t') for line in id_update.readlines() ] } # extract submitted ids

            stuck_tasks = [ job_id for job_id, info in running_tasks.items() if 'failed' in info[NODELIST] ] # find any stuck tasks
            relevant_tasks, tasks_to_resubmit = [ set(ids).intersection(submitted_tasks) for ids in [running_tasks, stuck_tasks] ] # find any relevant tasks

            if not relevant_tasks: break # exit if no relevant running tasks
            elif tasks_to_resubmit: # re-submit any relevant stuck tasks
                
                for job_id in sorted(tasks_to_resubmit):
                    
                    sh_file = submitted_tasks[job_id] # specify sh script to resubmit
                    task_name = running_tasks[job_id][NAME] # specify name of resubmitted task

                    subprocess.call(f'scancel {job_id}', shell=True) # cancel task
                    resub_id = CAPTURE(f'sbatch -d singleton {sh_file}') # resubmit script

                    subprocess.call(f'echo "{task_name} cancelled ({job_id}) & resubmitted ({resub_id})"', shell=True) # log cancelation & resubmission
                    print(resub_id, sh_file, '(RESUBMITTED)', sep='\t', file=id_update) # record resubmitted id

        time.sleep(600) # check every 10 minutes

subprocess.call('echo BUDDY COMPLETED `date`', shell=True) # log finish
