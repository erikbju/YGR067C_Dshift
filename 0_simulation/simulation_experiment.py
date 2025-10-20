from datetime import datetime
import logging
import socket
import os
from tqdm import tqdm
from random import sample, choice, seed
import argparse
import yaml
import json
import sys
from itertools import repeat
import pandas as pd

import multiprocessing as mp
from functools import partial

from cobra.io import read_sbml_model
from cobra.exceptions import OptimizationError

sys.path.append("/usr/src/lgemplus/")
from lgemcore.pathways import calculate_pathway

from utils.results_comparison import results_comparison
from utils.pathway_calculations import calculate_lgem_result, calculate_fba_result

# Seed random for reproducibility
seed(67)

def create_logger(metadata, level=logging.DEBUG):
    logger = logging.getLogger(__name__)
    logging.basicConfig(
        format="[%(asctime)s] %(name)s:%(levelname)s\t%(message)s",
        filename=os.path.join(metadata["output-directory"], "run.log"),
        filemode="w",
        encoding="utf-8",
        level=level,
    )
    return logger


def log_script_info(args):
    # Get the current time
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Get the machine name
    machine_name = socket.gethostname()

    # Create a dictionary with the logging information
    log_info = {
        "timestamp": current_time,
        "machine_name": machine_name,
        "arguments": vars(args),
    }

    return log_info


def validate_settings(settings):
    required_keys = [
        "name",
        "description",
        "inputs",
        "parameters",
        "output-root-directory",
    ]
    required_inputs = [
        "glucose-theory-directory",
        "ethanol-theory-directory",
        "fba-model-xml",
        "gene-kos",
        "reaction-kos",
        "hypothesis-strain",
        "hypothesis-reaction-set",
        "metabolomics-glucose-csv",
        "metabolomics-ethanol-csv",
        "transcriptomics-glucose-csv",
        "transcriptomics-ethanol-csv",
    ]
    required_parameters = [
        "min-reaction-subset-size",
        "max-reaction-subset-size",
        "extract-transcriptome",
        "extract-metabolome",
        "flux-threshold",
        "number-of-simulations",
    ]

    settings = metadata["settings"]

    # Top-level keys
    for key in required_keys:
        if key not in settings:
            raise KeyError(f"Key '{key}' missing at top level of settings YAML file")

    # Inputs keys
    for key in required_inputs:
        if key not in settings["inputs"]:
            raise KeyError(f"Key '{key}' missing in 'inputs' of settings YAML file")

    # Parameters keys
    for key in required_parameters:
        if key not in settings["parameters"]:
            raise KeyError(f"Key '{key}' missing in 'parameters' of settings YAML file")


def strip_whitespace(obj):
    if isinstance(obj, dict):
        return {k: strip_whitespace(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [strip_whitespace(v) for v in obj]
    elif isinstance(obj, str):
        return obj.strip()
    else:
        return obj


def create_directory(base_path, datetime_str):
    # Convert the datetime string to a datetime object
    dt = datetime.strptime(datetime_str, "%Y-%m-%d %H:%M:%S")

    # Create the directory path
    dir_path = os.path.join(base_path, f"full-{dt.strftime('%Y%m%d_%H%M%S')}")

    # Create the directory
    os.makedirs(dir_path, exist_ok=True)

    return dir_path


def create_output_directory(metadata):
    return create_directory(
        metadata["settings"]["output-root-directory"], metadata["runtime"]["timestamp"]
    )



## Credit for this function to leimao@github
def run_imap_multiprocessing(
    func, argument_list, num_processes, maxtasksperchild=None, **tqdm_kwargs
):
    with mp.Pool(processes=num_processes, maxtasksperchild=maxtasksperchild) as pool:
        result_list_tqdm = []
        for result in tqdm(
            pool.imap(func=func, iterable=argument_list),
            total=len(argument_list),
            **tqdm_kwargs,
        ):
            result_list_tqdm.append(result)

        return result_list_tqdm


def simulate_lgem_results(
    medium, metadata, res, cpus_to_use, max_tasks_per_child, logger=None
):
    # Get all reactions from file
    with open(
        os.path.join(
            metadata["settings"]["inputs"][f"{medium}-theory-directory"],
            "all_reactions.txt",
        ),
        "r",
    ) as fi:
        all_reactions = fi.read().splitlines()
    # print(len(all_reactions))

    # Get hypothesis reactions from file
    with open(metadata["settings"]["inputs"]["hypothesis-reaction-set"], "r") as fi:
        hypothesis_reactions = fi.read().splitlines()

    # Define partial functions
    hyp_map = partial(
        calculate_lgem_result,
        reaction_list=hypothesis_reactions,
        metadata=metadata,
        medium=medium,
        logger=logger,
    )

    all_map = partial(
        calculate_lgem_result,
        reaction_list=all_reactions,
        metadata=metadata,
        medium=medium,
        logger=logger,
    )

    hyp_map_i = lambda i: (i, hypothesis_map())
    all_map_i = lambda i: (i, all_map())

    # Conduct base strain perturbations
    res["base-all"] = run_imap_multiprocessing(
        all_map,
        range(metadata["settings"]["parameters"]["number-of-simulations"]),
        num_processes=cpus_to_use,
        maxtasksperchild=max_tasks_per_child,
        desc=f"lgem-{medium}-base-all",
    )

    # Conduct deletant perturbations
    strain = "-".join(metadata["settings"]["inputs"]["hypothesis-strain"])
    res[strain + "-hyp"] = run_imap_multiprocessing(
        hyp_map,
        range(metadata["settings"]["parameters"]["number-of-simulations"]),
        num_processes=cpus_to_use,
        maxtasksperchild=max_tasks_per_child,
        desc=f"lgem-{medium}-{strain}-hyp",
    )


def simulate_fba_results(
    fba_model, medium, metadata, results, cpus_to_use, max_tasks_per_child, logger=None
):
    # Set all reactions based on FBA model
    all_reactions = [r.id for r in fba_model.reactions]

    # Get hypothesis reactions from file
    with open(metadata["settings"]["inputs"]["hypothesis-reaction-set"], "r") as fi:
        hypothesis_reactions = fi.read().splitlines()

    # Define results function
    pathway_format = []
    pathway_format += (
        ["genes"] if metadata["settings"]["parameters"]["extract-transcriptome"] else []
    )
    pathway_format += (
        ["metabolites"]
        if metadata["settings"]["parameters"]["extract-metabolome"]
        else []
    )

    # Define partial functions
    hyp_map = partial(
        calculate_fba_result,
        fba_model=fba_model,
        reaction_list=hypothesis_reactions,
        metadata=metadata,
        medium=medium,
        logger=logger,
    )

    all_map = partial(
        calculate_fba_result,
        fba_model=fba_model,
        reaction_list=all_reactions,
        metadata=metadata,
        medium=medium,
        logger=logger,
    )

    hyp_map_i = lambda i: (i, hypothesis_map())
    all_map_i = lambda i: (i, all_map())

    # Conduct base strain perturbations
    # results["base-all"] = []
    # calculate_n_results(all_reactions, metadata, results["base-all"])
    results["base-all"] = run_imap_multiprocessing(
        all_map,
        range(metadata["settings"]["parameters"]["number-of-simulations"]),
        num_processes=cpus_to_use,
        maxtasksperchild=max_tasks_per_child,
        desc=f"fba-{medium}-base-all",
    )

    # Conduct deletant perturbations
    strain = "-".join(metadata["settings"]["inputs"]["hypothesis-strain"])
    # results[strain + "-hyp"] = []
    # calculate_n_results(hypothesis_reactions, metadata, results[strain + "-hyp"])
    results[strain + "-hyp"] = run_imap_multiprocessing(
        hyp_map,
        range(metadata["settings"]["parameters"]["number-of-simulations"]),
        num_processes=cpus_to_use,
        maxtasksperchild=max_tasks_per_child,
        desc=f"fba-{medium}-{strain}-hyp",
    )


if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(prog="SimulationExperimentRunner")
    parser.add_argument("settings_yaml")  # positional argument for settings yaml
    parser.add_argument(
        "-f", "--fba-only", action="store_true", help="Run only FBA experiments"
    )
    parser.add_argument(
        "-l", "--lgem-only", action="store_true", help="Run only LGEM experiments"
    )
    args = parser.parse_args()

    if args.fba_only and args.lgem_only:
        print("Error: Only one of `-f` or `-l` can be selected", file=sys.stderr)
        sys.exit(1)

    # Add runtime metadata
    metadata = {"runtime": log_script_info(args)}

    # Load settings from .yaml
    with open(args.settings_yaml, "r") as fi:
        metadata["settings"] = yaml.safe_load(fi)

    # Validate settings
    validate_settings(metadata["settings"])

    # Strip leading and trailing whitespac
    metadata["settings"] = strip_whitespace(metadata["settings"])

    # Create output directory
    metadata["output-directory"] = create_output_directory(metadata)
    os.makedirs(
        os.path.join(
            metadata["output-directory"],
            "raw_sim",
        )
    )
    os.makedirs(
        os.path.join(
            metadata["output-directory"],
            "sheets_as_csv",
        )
    )

    # Initiate logger
    logger = create_logger(metadata)

    # Determine number of CPUs to use
    CPUS_TO_USE = max(1, mp.cpu_count() - 1)
    MAX_TASKS_PER_CHILD = 10
    logger.debug(f"Using {CPUS_TO_USE} CPUs")

    # Copy settings yaml
    with open(os.path.join(metadata["output-directory"], "settings.yaml"), "w") as fo:
        with open(args.settings_yaml, "r") as fi:
            fo.write(fi.read())

    # Write runtime metadata to file
    with open(
        os.path.join(metadata["output-directory"], "runtime-metadata.json"), "w"
    ) as fo:
        json.dump(metadata["runtime"], fo)

    # Load FBA model
    fba_model = read_sbml_model(metadata["settings"]["inputs"]["fba-model-xml"])
    fba_model.solver.configuration.timeout = 5

    # Setup results dictionary
    results = {
        "lgem": {"glucose": {}, "ethanol": {}},
        "fba": {"glucose": {}, "ethanol": {}},
    }

    # LGEM+ experiments
    if args.fba_only:
        logger.debug("Skipping LGEM+ experiments")
    else:
        for medium in ("glucose", "ethanol"):
            # Simulate
            simulate_lgem_results(
                medium,
                metadata,
                results["lgem"][medium],
                CPUS_TO_USE,
                MAX_TASKS_PER_CHILD,
                logger=logger,
            )
            # Compare
            T, M = results_comparison(
                results["lgem"][medium],
                "LGEM",
                metadata["settings"]["inputs"][f"transcriptomics-{medium}-csv"],
                metadata["settings"]["inputs"][f"metabolomics-{medium}-csv"],
                theory_directory=metadata["settings"]["inputs"][
                    f"{medium}-theory-directory"
                ],
                fba_model_file=metadata["settings"]["inputs"]["fba-model-xml"],
            )
            T.to_csv(
                os.path.join(
                    metadata["output-directory"],
                    "sheets_as_csv",
                    f"lgem_{medium}_tra.csv",
                )
            )
            M.to_csv(
                os.path.join(
                    metadata["output-directory"],
                    "sheets_as_csv",
                    f"lgem_{medium}_met.csv",
                )
            )

    # FBA experiments
    if args.lgem_only:
        logger.debug("Skipping FBA experiments")
    else:
        for medium in ("glucose", "ethanol"):
            # Simulate
            simulate_fba_results(
                fba_model,
                medium,
                metadata,
                results["fba"][medium],
                CPUS_TO_USE,
                MAX_TASKS_PER_CHILD,
                logger=logger,
            )

            # Compare
            T, M = results_comparison(
                results["fba"][medium],
                "FBA",
                metadata["settings"]["inputs"][f"transcriptomics-{medium}-csv"],
                metadata["settings"]["inputs"][f"metabolomics-{medium}-csv"],
                fba_model_file=metadata["settings"]["inputs"]["fba-model-xml"],
            )
            T.to_csv(
                os.path.join(
                    metadata["output-directory"],
                    "sheets_as_csv",
                    f"fba_{medium}_tra.csv",
                )
            )
            M.to_csv(
                os.path.join(
                    metadata["output-directory"],
                    "sheets_as_csv",
                    f"fba_{medium}_met.csv",
                )
            )

    # Write out results to file
    for sim in ("lgem", "fba"):
        for medium in ("glucose", "ethanol"):
            with open(
                os.path.join(
                    metadata["output-directory"],
                    "raw_sim",
                    f"{sim}-{medium}-results.json",
                ),
                "w",
            ) as fo:
                json.dump(results[sim][medium], fo)

    # Combine results
    with pd.ExcelWriter(
        os.path.join(metadata["output-directory"], "all_results.xlsx")
    ) as writer:
        for csv in filter(
            lambda f: f.endswith(".csv"),
            os.listdir(os.path.join(metadata["output-directory"], "sheets_as_csv")),
        ):
            df = pd.read_csv(
                os.path.join(metadata["output-directory"], "sheets_as_csv", csv)
            )
            df.to_excel(writer, sheet_name=csv.replace(".csv", ""))
