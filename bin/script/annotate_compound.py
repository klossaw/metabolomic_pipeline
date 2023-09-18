# -*- coding: UTF-8 -*-
import pandas as pd
from tqdm import tqdm
from sqlite3 import connect
import json
import warnings
import argparse

warnings.filterwarnings("ignore")


def name_search(candidates, name, source):
    for _, candidate in candidates.iterrows():
        for candidate_name in candidate['name'].split("|"):
            if candidate_name == name:
                if source == "KEGG":
                    return candidate.entry, candidate.pathway, ""
                elif source == "HMDB":
                    return candidate.accession, candidate['super_class'], candidate['class']
    return "", "", ""


if __name__ == "__main__":
    ## usage: python annotate_compound.py --db ../database/metabolics.db --kegg ../database/pathway_tree2.json --peak-neg ../CD_neg.xlsx --peak-pos ../CD_pos.xlsx
    # --intensity-pos ../已完成数据/metaX/pos/data/pos-raw-metaboAnalystInput.csv --intensity-neg ../已完成数据/metaX/neg/data/neg-raw-metaboAnalystInput.csv
    # --output ../已完成数据/2.MetaboliteIdentification/identification.raw.intensity.csv --stats ../已完成数据/2.MetaboliteIdentification/identification.stat.csv
    parser = argparse.ArgumentParser("annotate compounds")
    parser.add_argument("--db", help="metabolics database", required=True)
    parser.add_argument("--kegg", help="kegg pathway", required=True)
    parser.add_argument("--peak-neg", help="input neg peak file, sheet name='Compounds'", required=True)
    parser.add_argument("--peak-pos", help="input pos peak file, sheet name='Compounds'", required=True)
    parser.add_argument("--intensity-pos", help="metaX pos intensity file", required=True)
    parser.add_argument("--intensity-neg", help="metaX neg intensity file", required=True)
    parser.add_argument("--sample-file", help="sample file", required=True)
    parser.add_argument("--output", help="output file, csv format", required=True)
    parser.add_argument("--stats", help="output statistics file, csv format", required=True)
    parser.add_argument("--mode", help="annotation mode, defualt is all", required=False, default="all")

    args = parser.parse_args()
    if args.mode == "all":
        exclude = []
    elif args.mode == "plant":
        exclude = ["Organismal Systems", "Human Diseases", "Drug Development"]
    else:
        raise Exception(f"unknown mode {args.mode}")
    ## use sqlite to accelerate search
    conn = connect(args.db)
    ## load sample file
    sample_file = pd.read_csv(args.sample_file, sep="\t", dtype=str)
    sample_file['class'].fillna("QC", inplace=True)
    sample_dict = dict(zip(sample_file['class'] + "_" + sample_file['order'].astype(str), sample_file['sample']))
    for file in [args.intensity_pos, args.intensity_neg]:
        with open(file) as f:
            if f.readline().startswith('"'):
                _tmp = pd.read_csv(file)
                _tmp.columns = [sample_dict.get(col, col) for col in _tmp.columns]
                _tmp.to_csv(file, index=False, header=True, sep=",")
    ## load CD_neg.xls
    total = []
    stats = []
    for peak_file, intensity_file, mode in [(args.peak_pos, args.intensity_pos, "pos"), (args.peak_neg, args.intensity_neg, "neg")]:

        print(f"load {peak_file}")
        cd_neg = pd.read_excel(peak_file, sheet_name="Compounds", engine='openpyxl', dtype=str)
        cd_neg.Formula = cd_neg.Formula.str.replace(" ", "")
        cd_neg.Name = cd_neg.Name.str.lower()
        cd_neg.fillna("", inplace=True)
        ## kegg annotation
        cd_neg['KEGG'] = ""
        cd_neg['Pathway'] = ""
        for idx, row in tqdm(cd_neg.iterrows(), desc=f"{peak_file} kegg annotation"):
            name, formula = row.Name, row.Formula
            candidates = pd.read_sql(f"SELECT * from kegg_compound WHERE formula = '{formula}'", conn)
            if candidates.shape[0] == 0:
                continue
            else:
                cd_neg.loc[idx, 'KEGG'], cd_neg.loc[idx, 'Pathway'], *_ = name_search(candidates, name, "KEGG")
        ## hmdb annotation
        cd_neg['HMDB'] = ""
        cd_neg['Superclass'] = ""
        cd_neg['Class'] = ""
        for idx, row in tqdm(cd_neg.iterrows(), desc=f"{peak_file} hmdb annotation"):
            name, formula = row.Name, row.Formula
            candidates = pd.read_sql(f"SELECT * from hmdb WHERE formula = '{formula}'", conn)
            if candidates.shape[0] == 0:
                continue
            else:
                cd_neg.loc[idx, 'HMDB'], cd_neg.loc[idx, 'Superclass'], cd_neg.loc[idx, 'Class'] = name_search(
                    candidates, name, "HMDB")

        ## Complementary annotation within KEGG and HMDB
        for idx, row in tqdm(cd_neg.iterrows(), desc=f"{peak_file} complementary annotation"):
            if row.KEGG == "" and row.HMDB == "":
                continue
            elif row.KEGG != "" and row.HMDB != "":
                continue
            elif row.KEGG == "":
                candidates = pd.read_sql(f"SELECT * from hmdb WHERE accession = '{row.HMDB}'", conn)
                if candidates.shape[0] > 0:
                    cd_neg.loc[idx, "KEGG"] = candidates.iloc[0,].kegg_id
                    kegg_candidates = pd.read_sql(
                        f"SELECT * from kegg_compound WHERE entry = '{candidates.iloc[0,].kegg_id}'", conn)
                    if kegg_candidates.shape[0] > 0:
                        cd_neg.loc[idx, "Pathway"] = kegg_candidates.iloc[0].pathway
            else:
                candidates = pd.read_sql(f"SELECT * from hmdb WHERE kegg_id = '{row.KEGG}'", conn)
                if candidates.shape[0] > 0:
                    cd_neg.loc[idx, "HMDB"] = candidates.iloc[0,].accession
                    cd_neg.loc[idx, 'Superclass'] = candidates.iloc[0,]['super_class']
                    cd_neg.loc[idx, 'Class'] = candidates.iloc[0,]['class']
        ## annotation with kegg pathway level
        pathway_tree = json.load(open(args.kegg))
        cd_neg['Level1'] = ""
        cd_neg['Level2'] = ""
        cd_neg['Level3'] = ""
        for idx, row in tqdm(cd_neg.iterrows(), desc=f"{peak_file} kegg pathway level"):
            if row.Pathway:
                _pathways = [_pathway for _pathway in row.Pathway.split("|") if pathway_tree.get(_pathway, {}).get('level1', "") not in exclude]
                cd_neg.loc[idx, 'Level1'] = "|".join(
                    [pathway_tree.get(_pathway, {}).get('level1', "") for _pathway in _pathways])
                cd_neg.loc[idx, 'Level2'] = "|".join(
                    [pathway_tree.get(_pathway, {}).get('level2', "") for _pathway in _pathways])
                cd_neg.loc[idx, 'Level3'] = "|".join(
                    [pathway_tree.get(_pathway, {}).get('level3', "") for _pathway in _pathways])
                cd_neg.loc[idx, "Pathway"] = "|".join(_pathways)
        output = cd_neg[
            ['Molecular Weight', 'RT [min]', 'Name', 'Formula', 'HMDB', 'Superclass', 'Class', 'KEGG', 'Level1',
             'Level2', 'Level3', 'Pathway']]
        output.columns = ["MZ", "RT", "Metabolite", 'Formula', 'HMDB', 'Superclass', 'Class', 'KEGG', 'Level1',
                          'Level2', 'Level3', 'Pathway']
        output['ID'] = output['RT'].astype(str) + "_" + output['MZ'].astype(str)
        output = output[output.columns.tolist()[-1:] + output.columns.tolist()[:-1]]

        # output.drop("MS2", axis=1, inplace=True)
        ## intensity data
        intensity_data = pd.read_csv(intensity_file, dtype=str)
        output = output.merge(intensity_data, left_on="ID", right_on="Sample")
        output['ID'] = f"{mode}-" + output['ID']
        output.drop("Sample", axis=1, inplace=True)
        output.fillna("", inplace=True)
        ## statistics for raw
        stats.append({"mode": mode, "All": output.shape[0], "Metabolite": output[output.Metabolite != ""].shape[0],
                      "HMDB": output[output.HMDB != ""].shape[0], "KEGG": output[output.KEGG != ""].shape[0]})
        
        total.append(output)

    total = pd.concat(total)
    total.to_csv(args.output, sep=",", index=False)
    stats = pd.DataFrame(stats)
    stats.to_csv(args.stats, sep=",", index=False)

