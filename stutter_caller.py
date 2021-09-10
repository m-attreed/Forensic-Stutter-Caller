
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
