"""
Lightweight script to verify the Hail environment is functional.

Initialises Hail, runs two in-memory checks (range Table and
Balding-Nichols MatrixTable), prints a success message, then stops Hail.
Exits with code 1 on any failure so it can be used as a pre-flight check
before expensive pipelines.

Usage examples:

    # Default settings (4 cores, 4 g driver memory, GRCh38)
    python script-utils/check_hail_env.py

    # Custom settings
    python script-utils/check_hail_env.py --n-cores 8 --driver-memory 8g --reference-genome GRCh37
"""

import argparse
import logging
import sys

import hail as hl

from utils.generic import hail_version

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("check_hail_env")
logger.setLevel(logging.INFO)


def main(args):
    logger.info(f"Hail version: {hail_version()}")

    logger.info(
        f"Initialising Hail (n_cores={args.n_cores}, "
        f"driver_memory={args.driver_memory}, "
        f"reference_genome={args.reference_genome}) ..."
    )
    hl.init(
        local=f"local[{args.n_cores}]",
        spark_conf={"spark.driver.memory": args.driver_memory},
        default_reference=args.reference_genome,
        quiet=True,
    )

    try:
        # Check 1: range Table
        logger.info("Running Table check (range_table(100)) ...")
        ht = hl.utils.range_table(100)
        n_rows = ht.count()
        assert n_rows == 100, f"Table check failed: expected 100 rows, got {n_rows}"
        logger.info(f"  Table check passed ({n_rows} rows)")

        # Check 2: Balding-Nichols MatrixTable
        logger.info(
            "Running MatrixTable check (balding_nichols_model: 3 pops, 10 samples, 100 variants) ..."
        )
        mt = hl.balding_nichols_model(n_populations=3, n_samples=10, n_variants=100)
        n_variants, n_samples = mt.count()
        assert n_variants == 100, (
            f"MatrixTable check failed: expected 100 variants, got {n_variants}"
        )
        assert n_samples == 10, (
            f"MatrixTable check failed: expected 10 samples, got {n_samples}"
        )
        logger.info(f"  MatrixTable check passed ({n_variants} variants x {n_samples} samples)")

        print("\u2713 Hail environment OK")
    finally:
        hl.stop()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Verify that the Hail environment is functional."
    )
    parser.add_argument(
        "--n-cores",
        help="Number of cores for local Hail session (default: 4)",
        type=int,
        default=4,
    )
    parser.add_argument(
        "--driver-memory",
        help="Driver memory for Hail session, e.g. 4g (default: 4g)",
        type=str,
        default="4g",
    )
    parser.add_argument(
        "--reference-genome",
        help="Default reference genome (default: GRCh38)",
        type=str,
        default="GRCh38",
    )

    args = parser.parse_args()

    try:
        main(args)
    except Exception as e:
        logger.error(f"Hail environment check failed: {e}")
        sys.exit(1)
