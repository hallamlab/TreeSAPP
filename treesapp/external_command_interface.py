import os
import sys
import logging
import subprocess
import multiprocessing

from tqdm import tqdm

from treesapp import logger

LOGGER = logging.getLogger(logger.logger_name())


def launch_write_command(cmd_list, collect_all=True):
    """
    Wrapper function for opening subprocesses through subprocess.Popen()

    :param cmd_list: A list of strings forming a complete command call
    :param collect_all: A flag determining whether stdout and stderr are returned
    via stdout or just stderr is returned leaving stdout to be written to the screen
    :return: A string with stdout and/or stderr text and the returncode of the executable
    """
    stdout = ""
    if collect_all:
        proc = subprocess.Popen(' '.join(cmd_list),
                                shell=True,
                                preexec_fn=os.setsid,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT)
        stdout = proc.communicate()[0].decode("utf-8")
    else:
        proc = subprocess.Popen(' '.join(cmd_list),
                                shell=True,
                                preexec_fn=os.setsid)
        proc.wait()

    # Ensure the command completed successfully
    if proc.returncode != 0:
        LOGGER.error(cmd_list[0] + " did not complete successfully! Command used:\n" +
                     ' '.join(cmd_list) + "\nOutput:\n" + stdout)
        sys.exit(19)

    return stdout, proc.returncode


def run_apply_async_multiprocessing(func, arguments_list: list, num_processes: int, pbar_desc: str,
                                    disable=False) -> list:
    if len(arguments_list) == 0:
        return []
    pool = multiprocessing.Pool(processes=num_processes)

    jobs = []
    result_list_tqdm = []
    pbar = tqdm(jobs, total=len(arguments_list), desc=pbar_desc, ncols=120, disable=disable)

    def update(*a):
        pbar.update()

    for args in arguments_list:
        jobs.append(pool.apply_async(func=func, args=(*args,), callback=update))
    pool.close()

    for job in pbar:
        result_list_tqdm.append(job.get())

    pbar.close()

    return result_list_tqdm
