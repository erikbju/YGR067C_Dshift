import os
import sys
import json
import pandas as pd
import numpy as np
from argparse import ArgumentParser
from cobra.io import read_sbml_model

sys.path.append("/usr/src/lgemplus/")
from lgemcore.io import logical_model_from_sbml
from lgemcore.pathways import calculate_pathway
from utils.pathway_calculations import calculate_pathway_fba


def results_list_to_df(l):
    return pd.DataFrame(
        columns=["i", "omic", "value1", "value2"],
        data=[(r[0], o[0], o[1][0], o[1][1]) for r in l for o in r[3]],
    )


def df_to_simp(df):
    return df[["i", "omic", "value1"]].drop_duplicates().reset_index(drop=True)


def simpdf_to_met_and_tran(df):
    return (
        df[df["omic"] == omic].groupby("value1")["i"].count().sort_values()
        for omic in ["metabolite", "gene"]
    )


def results_comparison(
    results,
    simulator,
    transcriptomics_file,
    metabolomics_file,
    theory_directory=None,
    fba_model_file="/usr/src/lgemplus/model-files/yeast9.xml",
):
    ref, comp = results.keys()
    refdf, compdf = map(results_list_to_df, results.values())

    met, tra = simpdf_to_met_and_tran(df_to_simp(compdf))
    n_successful = tra.max()

    if simulator.upper() == "LGEM":
        if theory_directory is None:
            raise ValueError("theory_directory must be provided for LGEM simulations")
        y9lm = logical_model_from_sbml(fba_model_file)
        gene_lookup = {g.word: g.sbml_id for g in y9lm.genes}
        kegg_lookup = {m.word: m.kegg_id for m in y9lm.metabolites}
        no_perturb = calculate_pathway(
            [], [], theory_directory, pathway_format=["genes", "metabolites"]
        )
    elif simulator.upper() == "FBA":
        fba_model = read_sbml_model(fba_model_file)
        gene_lookup = {g.id: g.id for g in fba_model.genes}
        kegg_lookup = {
            m.name: m.annotation.get("kegg.compound") for m in fba_model.metabolites
        }
        no_perturb = calculate_pathway_fba(
            fba_model, [], pathway_format=["genes", "metabolites"]
        )
    else:
        raise ValueError("Simulator must be 'LGEM' or 'FBA'")

    wttra = pd.DataFrame(index=[t[1][0] for t in no_perturb[2] if t[0] == "gene"])
    wttra["wt"] = 1

    diff = wttra.join(tra, how="outer").fillna(0)
    diff.index = [i.replace("__45__", "_").replace("_ARG5", "") for i in diff.index]
    diff["i"] = diff["i"] / n_successful
    diff["diff"] = diff.i - diff.wt

    wtmet = pd.DataFrame(index=[t[1][0] for t in no_perturb[2] if t[0] == "metabolite"])
    wtmet["wt"] = 1

    diffmet = wtmet.join(met, how="outer").fillna(0)
    diffmet["i"] = diffmet["i"] / n_successful
    diffmet["diff"] = diffmet.i - diffmet.wt

    all_genes = pd.DataFrame(
        data=[{"sbml_id": gene_lookup[g], "word": g} for g in gene_lookup]
    ).set_index("sbml_id")
    preds_tra = all_genes.join(diff, how="outer").fillna(0)

    all_met = pd.DataFrame(
        data=[{"keggid": kegg_lookup[k], "word": k} for k in kegg_lookup]
    ).set_index("keggid")
    preds_met = all_met.join(diffmet, on="word", how="outer").fillna(0)
    preds_met = preds_met.loc[list(filter(lambda i: i is not None, preds_met.index)), :]

    tra_true = pd.read_csv(transcriptomics_file, delimiter="\t").set_index(
        "Gene_stable_ID"
    )
    met_true = pd.read_csv(metabolomics_file, delimiter="\t").set_index("KEGG")

    joint_tra = preds_tra.join(tra_true, how="left").drop_duplicates()
    joint_met = preds_met.join(met_true, how="left").drop_duplicates()

    joint_tra["pred_dir"] = joint_tra["diff"].map(np.sign)
    joint_tra["act_dir"] = joint_tra["log2FoldChange"].map(np.sign)

    joint_met["pred_dir"] = joint_met["diff"].map(np.sign)
    joint_met["act_dir"] = joint_met["l2fc"].map(np.sign)

    return (joint_tra, joint_met)


if __name__ == "__main__":
    parser = ArgumentParser(prog="ExperimentResultsComparison")
    parser.add_argument("results_json")  # positional argument for results json
    parser.add_argument("-t", "--theory_directory")  # argument theory directory
    parser.add_argument("-x", "--model_xml")  # argument for model .xml file
    parser.add_argument("-T", "--transcriptomics")
    parser.add_argument("-M", "--metabolomics")
    args = parser.parse_args()

    with open(args.results_json, "r") as fi:
        results = json.load(fi)

    simulator = os.path.basename(args.results_json).split("-")[0].upper()

    T, M = results_comparison(
        results,
        simulator,
        args.transcriptomics,
        args.metabolomics,
        fba_model_file="/usr/src/lgemplus/model-files/yeast9.xml",
    )

    T.to_csv("lgem_ethanol_tra.csv")
    M.to_csv("lgem_ethanol_met.csv")
