__author__ = 'Connor Morgan-Lang'

import os
import sys
import logging
import subprocess


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
        logging.error(cmd_list[0] + " did not complete successfully! Command used:\n" +
                      ' '.join(cmd_list) + "\nOutput:\n" + stdout)
        sys.exit(19)

    return stdout, proc.returncode


def setup_progress_bar(num_items):
    if num_items > 50:
        progress_bar_width = 50
        step_proportion = float(num_items) / progress_bar_width
    else:
        progress_bar_width = num_items
        step_proportion = 1

    sys.stdout.write("[%s]" % (" " * progress_bar_width))
    sys.stdout.write("\b" * (progress_bar_width+1))
    sys.stdout.flush()

    return step_proportion
