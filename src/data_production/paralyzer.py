import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
from src.data_production.chunk import chunk
from pathlib import Path
from datetime import datetime, timedelta
import concurrent.futures
from src.utilities.misc import T_print


def paralyze(runs, total_production_years):
    for run in runs:
        with concurrent.futures.ProcessPoolExecutor() as executor:
            results = []
            for i, y in enumerate(total_production_years[:]):
                if i%3 == 0:
                    production_years = total_production_years[i:i+3]
                    test_name = "{}_{}".format(production_years[0], production_years[-1])
                    f_path = Path("data/results/parallel_processes/{}/{}".format(run, test_name))
                    results.append(executor.submit(chunk, f_path, production_years, test_name))
            
            for f in concurrent.futures.as_completed(results):
                print(f.result())

        total_p_days = 0
        total_impactors = 0
        total_p = 0 
        chunk_times = []
        for f in results:
            p_days_run, chunk_time, num_ps, num_impactors = f.result()
            total_p_days += p_days_run
            chunk_times.append(chunk_time.total_seconds())
            total_p += num_ps
            total_impactors += num_impactors
        end = datetime.now()

        with open("data/results/parallel_processes/{}/report.txt".format(run), "w")as txt:
            T_print("Full Sim Ran in {} seconds".format((end-start).total_seconds()), txt)
            T_print("A total of {} ParticleDays were run".format(total_p_days), txt)
            T_print("Average Chunk Time: {}".format(sum(chunk_times)/len(chunk_times)), txt)

            T_print("{} Impactors from {} Total Particles".format(total_impactors, total_p), txt)


if __name__ == '__main__':
    total_production_years = list(range(1900, 2030))
    start = datetime.now()
    #runs = ["sawg3_{}".format(i) for i in range(100)]
    runs = ["test"]
    paralyze(runs, total_production_years)