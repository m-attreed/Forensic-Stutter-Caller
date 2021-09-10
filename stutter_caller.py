# -*- coding: utf-8 -*-
"""
This program adds stutter call information to a GeneMarker LIMS report export
file. The program takes as input the Genemarker LIMS report file and a file
containing profile data. It expects both files to be in tab-separated values
(tsv) format.

The program appends an additional value to the end of each line of the input
file. These values can be as follows: 1, 2, del1, del2, X, db1, db2, b1, b2,
hb1, hb2, f1, and f2. 1 and 2 refer to the parent peaks (profile alleles)
of the profile at a particular locus i.e. 12, 16. The del1 and del2 values
indicate the parent peak was not called due to a reason outlined in
the stutter study protocol. An X indicates the allele call is not
a parent peak nor a called stutter peak. The remaining values refer to stutter
with the number at the end indicating which parent peak the stutter originated
from. A db is a double back stutter; b is back stutter;
hb is halfback stutter; and f is forward stutter.

The program writes the output to a new file in the tab-separated values (tsv)
format. The output excludes the Allelic Ladder and Amp Neg data.

The program is designed to run from the folder containing both the LIMS report
file and the profiles file. The names of these files are set in the software
there is no command line input.
"""


import os
import sys
import csv
import time
import copy

import strlibrary


# Start a timer to measure speed
start = time.perf_counter()


def main():
    """
    The main loop first sets the input data to be used by the analyzeAll
    function.
    """
    os.chdir(os.path.dirname(sys.argv[0]))

    input_directory = "Helen_GM_AT"
    profiles_file_name = 'profiles_3500.tsv'

    files = os.listdir(input_directory)

    for file in files:
        if not file.startswith("."):
            if os.path.isdir(input_directory + "/" + file):
                continue
            else:
                profile_db = strlibrary.ProfileDB(profiles_file_name)
                profile_db.addMixes("Mixtures.tsv")

                report_db = strlibrary.ReportDB(input_directory + "/" + file)
                report_db.mark_parent_peaks(profile_db)
                report_db.mark_stutter(profile_db)
                report_db.mark_pullup()
                report_db.write_output()


    # stop the timer to report the total time it took to complete
    # and report this time to the command line
    finish = time.perf_counter()

    print(f'Finished in {round(finish - start, 2)} second(s)')


if __name__ == '__main__':
    main()
