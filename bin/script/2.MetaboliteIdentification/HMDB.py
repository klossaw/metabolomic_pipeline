# -*- coding: utf-8 -*-
import pandas as pd
import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import argparse
import matplotlib

matplotlib.use('AGG')
plt.figure(constrained_layout=True)


def HMDB_all_bar_plot(data, outdir):
    super_class = list(set(data["Superclass"]))
    count = [len(meta_id_data[meta_id_data["Superclass"] == x]) for x in super_class]
    plt.figure()
    p, ax = plt.subplots(figsize=(9, 5))
    ax.bar(x=super_class, height=count, color=plt.cm.get_cmap("tab20")(range(len(super_class))))
    ax.set_ylabel("Number of identified feature")
    ax.set_xlabel("HMDB super class")
    plt.xticks(super_class, fontsize=8, ha='right', rotation=30)
    plt.savefig(f"{os.path.join(outdir, 'countClass.png')}", bbox_inches='tight')
    plt.close()


def HMDB_countMZ_RT_plot(data, outdir):
    superclass = list(set(data["Superclass"]))
    color_dict = {superclass[x]: x for x in range(len(superclass))}
    re_color_list = [list(color_dict.keys())[x] for x in color_dict.values()]

    RT = data["RT"].values
    MZ = data["MZ"].values
    super_class = data["Superclass"].values
    color = [color_dict[super_class[x]] for x in range(len(super_class))]
    plt.figure()
    p, ax = plt.subplots(figsize=(11, 6))
    scatter = ax.scatter(x=RT, y=MZ, c=color, cmap=plt.cm.Paired)
    ax.set_ylabel("MZ")
    ax.set_xlabel("RT")
    handles, labels = scatter.legend_elements(prop="colors", num=len(superclass))
    labels_0 = re_color_list
    ax.legend(handles, labels_0, loc=[1.05, 0], title="Super class")
    plt.savefig(f"{os.path.join(outdir, 'countMZ-RT.png')}", bbox_inches='tight')
    plt.close()


def HMDB_countMZ_RT_facet_plot(data: pd.DataFrame, outdir):
    # 数据分割
    data_pos = data[data["ID"].str.contains("pos", regex=False)].sort_values("Superclass")
    data_neg = data[data["ID"].str.contains("neg", regex=False)].sort_values("Superclass")

    # 阳性数据整理
    RT_pos = data_pos["RT"].values
    MZ_pos = data_pos["MZ"].values
    super_class_pos = data_pos["Superclass"].values
    color_pos = [color_dict[super_class_pos[x]] for x in range(len(super_class_pos))]

    # 阴性数据整理
    RT_neg = data_neg["RT"].values
    MZ_neg = data_neg["MZ"].values
    super_class_neg = data_neg["Superclass"].values
    color_neg = [color_dict[super_class_neg[x]] for x in range(len(super_class_neg))]

    fig = plt.figure(figsize=(9, 6))
    gs = GridSpec(2, 6, figure=fig)
    ax = fig.add_subplot(gs[0, :-2])
    ax1 = fig.add_subplot(gs[1, :-2])
    scatter1 = ax.scatter(x=RT_neg, y=MZ_neg, c=color_neg, cmap=plt.cm.tab20c, label="negative")
    scatter2 = ax1.scatter(x=RT_pos, y=MZ_pos, c=color_pos, cmap=plt.cm.tab20c, label="positive")
    ax.set_title("negative")
    ax.set_ylabel("MZ")
    ax1.set_title("positive")
    handles1, labels1 = scatter1.legend_elements(prop="colors", num=len(superclass))
    handles2, labels2 = scatter2.legend_elements(prop="colors", num=len(superclass))
    labels_0 = re_color_list
    ax1.set_xlabel("Retention Time(min)")
    ax1.set_ylabel("MZ")

    plt.figlegend(handles1, labels_0, bbox_to_anchor=(1.02, 0.9), title="Super class", frameon=False, fontsize=9,
                  labelspacing=1.2)
    plt.savefig(f"{os.path.join(outdir, 'countMZ-RT-facet.png')}", bbox_inches='tight')
    plt.close()


def HMDB_countRT_plot(data: pd.DataFrame, outdir):
    # 获取各个时间段的各样本的count
    RT_count_df = pd.DataFrame()
    for i in superclass:
        superclass_data = data[data["Superclass"] == i]
        superclass_count_list = [
            len(superclass_data[(superclass_data["RT"] < x + 1) & (superclass_data["RT"] > x)]) if x < 8 else len(
                superclass_data[superclass_data["RT"] > x]) for x in range(9)]
        RT_count_df[i] = superclass_count_list

    x_ticks = [f"{i}-{i + 1}" if i < 8 else f">{i}" for i in range(9)]
    # 开始画图
    plt.figure()
    p, ax = plt.subplots(figsize=(10, 6))
    for i in range(len(superclass)):
        exec(f"var_{i} = RT_count_df[superclass[{i}]].values")
        if i == 0:
            ax.bar(x_ticks, eval(f"var_0"), label=superclass[i])
            bottom = eval(f"var_0")
        else:
            ax.bar(x_ticks, eval(f"var_{i}"), bottom=bottom, label=superclass[i])
            bottom = bottom + eval(f"var_{i}")
    ax.set_xlabel("Retention times(min)")
    ax.set_ylabel("Number of identified feature")
    plt.xticks(rotation=45, ha="right")
    ax.legend(loc=(1.05, 0), title="Super class")
    plt.savefig(f"{os.path.join(outdir, 'countRT.png')}", bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    ## usage: python 2.MetaboliteIdentification/HMDB.py --input ../已完成数据/2.MetaboliteIdentification/identification.raw.intensity.csv
    # --outdir ../已完成数据/2.MetaboliteIdentification/HMDB
    parser = argparse.ArgumentParser("plot HMDB class distribution")
    parser.add_argument("--input", help="annotated peak file", required=True)
    parser.add_argument("--outdir", help="", required=True)

    args = parser.parse_args()
    meta_id_data = pd.read_csv(args.input)
    # filter unknown data
    meta_id_data = meta_id_data[meta_id_data["Superclass"].notna()]

    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    data = meta_id_data
    superclass = list(set(data["Superclass"]))
    global color_dict, re_color_list

    color_dict = {superclass[x]: x for x in range(len(superclass))}
    re_color_list = [list(color_dict.keys())[x] for x in color_dict.values()]

    HMDB_all_bar_plot(meta_id_data, outdir)
    HMDB_countMZ_RT_plot(meta_id_data, outdir)
    HMDB_countRT_plot(meta_id_data, outdir)
    HMDB_countMZ_RT_facet_plot(meta_id_data, outdir)

