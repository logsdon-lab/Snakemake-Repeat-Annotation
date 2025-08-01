from os.path import join, dirname


SAMPLES = config["samples"].keys()
OUTPUT_DIR = config.get("output_dir", "results")
LOG_DIR = config.get("log_dir", "logs")
BMK_DIR = config.get("benchmark_dir", "benchmarks")


wildcard_constraints:
    sm="|".join(SAMPLES),
    fname=r"[^/]*",
