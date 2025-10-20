from tqdm import tqdm
from random import sample, choice, seed
import sys
from itertools import repeat

from cobra.io import read_sbml_model
from cobra.exceptions import OptimizationError

sys.path.append("/usr/src/lgemplus/")
from lgemcore.pathways import calculate_pathway

seed(67)

def calculate_lgem_result(
    i,
    reaction_list,
    metadata,
    medium="glucose",
    pathway_format=["metabolites", "genes"],
    logger=None,
):
    sample_size = min(
        choice(
            range(
                metadata["settings"]["parameters"]["min-reaction-subset-size"],
                metadata["settings"]["parameters"]["max-reaction-subset-size"] + 1,
            )
        ),
        len(reaction_list),
    )
    ko_rxns = sample(reaction_list, sample_size)
    # print(sample_size, ": ", ko_rxns, sep="", file=sys.stdout)
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
    # # Write to the log file the reactions being knocked out
    # if logger is not None:
    #     logger.debug(f"Simulation {i} (LGEM+, {medium}): KO reactions: {', '.join(ko_rxns)}")
    return tuple(
        [i]
        + list(
            calculate_pathway(
                [],
                ko_rxns,
                metadata["settings"]["inputs"][f"{medium}-theory-directory"],
                pathway_format=pathway_format,
            )
        )
    )


def calculate_n_lgem_results(reaction_list, metadata, results_list, medium="glucose"):
    for i in tqdm(range(metadata["settings"]["parameters"]["number-of-simulations"])):
        sample_size = min(
            choice(
                range(
                    metadata["settings"]["parameters"]["min-reaction-subset-size"],
                    metadata["settings"]["parameters"]["max-reaction-subset-size"] + 1,
                )
            ),
            len(reaction_list),
        )
        ko_rxns = sample(reaction_list, sample_size)
        # print(sample_size, ": ", ko_rxns, sep="", file=sys.stdout)
        results_list.append(
            tuple(
                [i]
                + list(
                    calculate_pathway(
                        [],
                        ko_rxns,
                        metadata["settings"]["inputs"][f"{medium}-theory-directory"],
                        pathway_format=pathway_format,
                    )
                )
            )
        )


def calculate_pathway_fba(
    fba_base_model,
    ko_rxns,
    medium="glucose",
    pathway_format=["metabolites", "genes"],
    flux_threshold=1.0e-9,
):
    results = []
    try:
        S = calculate_fluxes_fba(fba_base_model, ko_rxns, medium)

        if "metabolites" in pathway_format:
            results.extend(
                zip(
                    repeat("metabolite"),
                    extract_metabolites_fba(
                        fba_base_model, ko_rxns, S, flux_threshold=flux_threshold
                    ),
                )
            )

        if "genes" in pathway_format:
            results.extend(
                zip(
                    repeat("gene"),
                    extract_genes_fba(
                        fba_base_model, ko_rxns, S, flux_threshold=1.0e-9
                    ),
                )
            )
    except OptimizationError:
        None

    return ([], ko_rxns, results)


def calculate_fluxes_fba(fba_base_model, ko_rxns, medium="glucose"):
    with fba_base_model:
        for r in ko_rxns:
            fba_base_model.reactions.get_by_id(r).knock_out()
        if medium == "ethanol":
            medium = fba_base_model.medium
            medium["r_1714"] = 0.0  # Switch off glucose (it is on by default)
            medium["r_1761"] = 1.0  # Switch on ethanol
            fba_base_model.medium = medium
        return fba_base_model.optimize(raise_error=True)


def extract_metabolites_fba(fba_base_model, ko_rxns, solution, flux_threshold=1.0e-9):
    S = solution
    if S.status != "optimal":
        return []
    rxns = S.fluxes[S.fluxes.abs() > flux_threshold].index
    return list(
        set(
            [
                (m.name, m.annotation.get("kegg.compound"))
                for r in rxns
                for m in fba_base_model.reactions.get_by_id(r).metabolites
                if r not in ko_rxns
            ]
        )
    )


def extract_genes_fba(fba_base_model, ko_rxns, solution, flux_threshold=1.0e-9):
    S = solution
    if S.status != "optimal":
        return []
    rxns = S.fluxes[S.fluxes.abs() > flux_threshold].index
    return list(
        set(
            [
                (g.id, g.name)
                for r in rxns
                for g in fba_base_model.reactions.get_by_id(r).genes
                if r not in ko_rxns
            ]
        )
    )


def calculate_fba_result(
    i, reaction_list, fba_model, metadata, medium="glucose", logger=None
):
    sample_size = min(
        choice(
            range(
                metadata["settings"]["parameters"]["min-reaction-subset-size"],
                metadata["settings"]["parameters"]["max-reaction-subset-size"] + 1,
            )
        ),
        len(reaction_list),
    )
    ko_rxns = sample(reaction_list, sample_size)
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
    # # print(sample_size, ": ", ko_rxns, sep="", file=sys.stdout)
    # if logger is not None:
    #     logger.debug(f"Simulation {i} (FBA, {medium}): KO reactions: {', '.join(ko_rxns)}")
    return tuple(
        [i]
        + list(
            calculate_pathway_fba(
                fba_model,
                ko_rxns,
                medium=medium,
                pathway_format=pathway_format,
                flux_threshold=metadata["settings"]["parameters"]["flux-threshold"],
            )
        )
    )


def calculate_n_fba_results(
    fba_model, reaction_list, metadata, results_list, medium="glucose"
):
    for i in tqdm(range(metadata["settings"]["parameters"]["number-of-simulations"])):
        results_list.append(
            calculate_fba_result(i, reaction_list, fba_model, metadata, medium=medium)
        )
