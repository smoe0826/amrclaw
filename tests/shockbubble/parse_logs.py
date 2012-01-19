#!/usr/bin/env python

r"""Parse a log file generated by run_thread_tests.py"""

import sys
import os
import re
import glob
import math
import getopt
    
import matplotlib.pyplot as plt

import run_thread_tests

help_message = r"""
Parse log files generated by run_thread_tests.py

Command line options:
    -v, --verbose - Verbose output (default == False)
    -p, --path=(string) - Path to directory containing log files 
                          (default == LOG_FILE_BASE in run_thread_tests.py)
    -f, --format=(string) - Format of plot output, can be anything that
                            matplotlib can understand (default == 'pdf')
    -h, --help - Display help message
"""

log_regex = re.compile(r"\*{3}\sOMP_NUM_THREADS\s=\s.*|\*{7}\stick\stiming\s=\s+.*\ss{1}")
time_regex = re.compile(r"\*{3}\sOMP_NUM_THREADS\s=\s.*|\s+User\sTime\s\(seconds\)\:\s.*")

def parse_log_file(path,verbose=True):
    # Open contents of log file
    log_contents = open(path,'r').read()
    
    # Loop through regular expression finds
    threads = []
    times = []
    for match in log_regex.finditer(log_contents):
        if "OMP_NUM_THREADS" in match.group():
            threads.append(int(match.group().strip("*** OMP_NUM_THREADS =")))
            if verbose:
                print "Threads = %s" % threads[-1]
        elif "tick timing" in match.group():
            times.append(float(match.group().strip("******* tick timing =")[:-1]))
            if verbose:
                print "Time = %s" % times[-1]
    
    if not len(threads) == len(times):
        print "Parsing may not have been successful, len(threads) != len(times)."
        print "  path = %s" % path
        print "  threads = %s" % threads
        print "  times = %s" % times
        sys.exit(2)
    
    return threads,times

def parse_time_file(path,verbose=True):
    
    threads = []
    times = []
    
    # Open contents of log file
    log_contents = open(path,'r').read()
    
    # Loop through regular expression finds
    for match in time_regex.finditer(log_contents):
        if "OMP_NUM_THREADS" in match.group():
            threads.append(int(match.group().strip("*** OMP_NUM_THREADS =")))
            if verbose:
                print "Threads = %s" % threads[-1]
        elif "User time" in match.group():
            times.append(float(match.group().strip().strip("User time (seconds): ")))
            if verbose:
                print "Time = %s" % times[-1]
    
    if not len(threads) == len(times):
        print "Parsing may not have been successful, len(threads) != len(times)."
        print "  path = %s" % path
        print "  threads = %s" % threads
        print "  times = %s" % times
        sys.exit(2)
    
    return threads,times

def create_timing_plots(log_paths,plot_file=None,verbose=False):
    # Setup plots
    rows = int(math.ceil(len(log_paths) / 2.0))
    fig = plt.figure(figsize=(10,12))
    
    for (i,path) in enumerate(log_paths):
        # Parse the log file
        if verbose:
            print os.path.basename(path)[:-4]
        if os.path.basename(path)[0:3] == "log":
            threads,times = parse_log_file(path,verbose=verbose)
        elif os.path.basename(path)[0:4] == "time":
            threads,times = parse_time_file(path,verbose=verbose)
        if verbose:
            print threads,times

        # Plot this run
        axes = fig.add_subplot(rows,2,i+1)
        axes.plot(threads,times,'or-')
        
        # Labeling
        axes.set_xbound(threads[0]-0.5,threads[-1]+0.5)
        if os.path.basename(path)[0:3] == "log":
            axes.set_title(os.path.basename(path).strip('log_')[:-4])
        elif os.path.basename(path)[0:4] == "time":
            axes.set_title(os.path.basename(path).strip('time_')[:-4])
        axes.set_xlabel('Number of Threads')
        axes.set_xticks(threads)
        axes.set_xticklabels(threads)
        axes.set_ylabel('Time (s)')

    plt.tight_layout()    
    if plot_file is not None:
        plt.savefig(plot_file)
    else:
        plt.show()

class Usage(Exception):
    def __init__(self,msg):
        self.msg = msg   

if __name__ == "__main__":
    # Parse command line arguments
    try:
        try:
            opts,args = getopt.getopt(sys.argv[1:],"hvp:f:",['help','verbose','path=','format='])
        except getopt.error, msg:
            raise Usage(msg)
        
        # Default values
        verbose = False
        log_dir = run_thread_tests.LOG_PATH_BASE
        format = 'pdf'
    
        # Option parsing
        for option,value in opts:
            if option in ("-v","--verbose"):
                verbose = True
            if option in ("-p",'--path'):
                log_dir = value
            if option in ('-f','--format'):
                format = value
            if option in ("-h","--help"):
                raise Usage(help_message)
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        # print >> sys.stderr, "\t for help use --help"
        sys.exit(2)
        
    # Find log files
    log_files = glob.glob(os.path.join(log_dir,"log*.txt"))
    if len(log_files) == 0:
        print >> sys.stderr, "Did not find any log files at:"
        print >> sys.stderr, "\t%s" % log_dir
        sys.exit(2)
        
    # Find and parse all log files found in log_dir
    create_timing_plots(log_files,plot_file='./tick_plots.%s' % format,verbose=verbose)
    
    # Only use the timing files if we are not on Darwin (time does not work as awesome there)
    if not os.uname()[0] == 'Darwin':
        time_files = glob.glob(os.path.join(log_dir,"time*.txt"))
        if len(log_files) == 0:
            print >> sys.stderr, "Did not find any time files at:"
            print >> sys.stderr, "\t%s" % log_dir
            sys.exit(2)
        create_timing_plots(time_files,plot_file='./time_plots.%s' % format,verbose=verbose)

        