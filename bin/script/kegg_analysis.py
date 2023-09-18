import argparse
import json
import scipy.stats as stats
import statsmodels.api as sm
import numpy as np
from tqdm import tqdm
import pandas as pd
import os
import shutil

if __name__ == "__main__":
    # usage: python kegg_analysis.py --input ../已完成数据/2.MetaboliteIdentification/identification.raw.intensity.csv --compound ../database/compounds.json
    # --outdir ../已完成数据/2.MetaboliteIdentification/KEGG
    parser = argparse.ArgumentParser("KEGG analysis")
    parser.add_argument("--input", help="indensity file", required=True)
    parser.add_argument("--compound", help="compound file", required=True)
    parser.add_argument("--outdir", help="output dir", required=True)
    parser.add_argument("--mode", help="mode, default is all", required=False, default="all")

    args = parser.parse_args()

    if args.mode == "all":
        exclude = []
    elif args.mode == "plant":
        exclude = ["Organismal Systems", "Human Diseases", "Drug Development"]
    else:
        raise Exception(f"unknown mode {args.mode}")

    indensity_data = pd.read_csv(args.input)
    indensity_data = indensity_data[indensity_data.Level3.notna()]
    pathway_statistics = {}
    for _, row in tqdm(indensity_data.iterrows()):
        level1 = row.Level1.split("|")
        level2 = row.Level2.split("|")
        level3 = row.Level3.split("|")
        for idx, pathway in enumerate(row.Pathway.split("|")):
            if pathway in pathway_statistics:
                pathway_statistics[pathway]["Compound"].append(row.KEGG)
                pathway_statistics[pathway]["Feature"].append(row.ID)
            else:
                pathway_statistics[pathway] = {
                    "Level1": level1[idx],
                    "Level2": level2[idx],
                    "Pathway": level3[idx],
                    "KEGG": pathway,
                    "Compound": [row.KEGG],
                    "Feature": [row.ID]
                }
    for pathway, item in pathway_statistics.items():
        item['NumberCompound'] = len(set(item['Compound']))
        item['NumberFeature'] = len(set(item['Feature']))
        item['Compound'] = ";".join(set(item['Compound']))
        item['Feature'] = ";".join(set(item['Feature']))

    pathway_statistics = pd.DataFrame(pathway_statistics.values())
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(os.path.join(args.outdir, "kegg_map"), exist_ok=True)
    if pathway_statistics.shape[0] > 0:
        pathway_statistics = pathway_statistics[
            ["Level2", "Level1", "Pathway", "KEGG", "NumberCompound", "NumberFeature", "Compound", "Feature"]]

        level2_data = []
        for level, group in pathway_statistics.groupby("Level2"):
            level2 = group.Level2.iloc[0]
            level1 = group.Level1.iloc[0]
            pathway = ";".join(group.Pathway)
            kegg = ";".join(group.KEGG)
            compound = ";".join(set(";".join(group.Compound).split(";")))
            feature = ";".join(set(";".join(group.Feature).split(";")))
            num_compound = len(compound.split(";"))
            num_feature = len(feature.split(";"))
            level2_data.append({"Level2": level2, "Level1": level1, "Pathway": pathway, "KEGG": kegg, "NumberCompound": num_compound,
                                "NumberFeature": num_feature, "Compound": compound, "Feature": feature})
        level2_data = pd.DataFrame(level2_data)
        level2_data = level2_data[~level2_data.Level1.isin(exclude)]
        level2_data.to_csv(os.path.join(args.outdir, "kegg.level2.csv"), index=False)

        ## enrichment analysis

        pathway_with_compound = {}
        for line in open(args.compound):
            item = json.loads(line.strip())
            for pathway in item['pathways'].keys():
                if pathway in pathway_with_compound:
                    pathway_with_compound[pathway].add(item['entry_id'])
                else:
                    pathway_with_compound[pathway] = set([item['entry_id']])

        detect_compound = indensity_data.KEGG.unique()
        total_compounds = [json.loads(line.strip())['entry_id'] for line in open(args.compound)]

        pathway_statistics['Background'] = None
        pathway_statistics['Pvalue'] = None
        for idx, row in tqdm(pathway_statistics.iterrows()):
            detect_compound_in_pathway = row.NumberCompound
            detect_compound_in_pathway_not_detect = len(
                set(pathway_with_compound[row.KEGG]) - set(row.Compound.split(";")))
            detect_compound_not_in_pathway = len(detect_compound) - detect_compound_in_pathway
            total_compound_not_in_pathway_not_detect = len(
                set(total_compounds) - set(detect_compound) - set(pathway_with_compound[row.KEGG]))
            contingency_table = np.array([[detect_compound_in_pathway, detect_compound_in_pathway_not_detect],
                                          [detect_compound_not_in_pathway, total_compound_not_in_pathway_not_detect]])
            odds_ratio, pvalue = stats.fisher_exact(contingency_table, alternative="greater")
            pathway_statistics.loc[idx, 'Background'] = len(pathway_with_compound[row.KEGG])
            pathway_statistics.loc[idx, 'Pvalue'] = pvalue

        pathway_statistics['FDR'] = sm.stats.multipletests(pathway_statistics.Pvalue, method="fdr_bh")[1]

        pathway_statistics = pathway_statistics[
            ["Level2", "Level1", "Pathway", "KEGG", "NumberCompound", "NumberFeature", "Background", "Pvalue", "FDR", "Compound", "Feature"]]

        ## kegg map dir
        pathway_statistics=pathway_statistics[~pathway_statistics.Level1.isin(exclude)]
        kegg_map_dir = os.path.join(os.path.dirname(args.compound), "KEGG")
        for idx, row in pathway_statistics.iterrows():
            if os.path.exists(os.path.join(kegg_map_dir, row.KEGG + ".png")):
                shutil.copyfile(os.path.join(kegg_map_dir, row.KEGG + ".png"), os.path.join(args.outdir, "kegg_map", row.KEGG + ".png"))
        pathway_statistics.to_csv(os.path.join(args.outdir, "kegg.csv"), index=False)
    else:
        open(os.path.join(args.outdir, "kegg.level2.csv"), "w").write(",".join(
            ["Level2", "Level1", "Pathway", "KEGG", "NumberCompound", "NumberFeature", "Compound", "Feature"]) + "\n")
        open(os.path.join(args.outdir, "kegg.csv"), "w").write(",".join(
            ["Level2", "Level1", "Pathway", "KEGG", "NumberCompound", "NumberFeature", "Background", "Pvalue", "FDR", "Compound", "Feature"]) + "\n")

