# 2. MetaboliteIdentification
## compound annotation
python annotate_compound.py --db ../database/metabolics.db --kegg ../database/pathway_tree2.json --peak-neg ../CD_neg.xlsx --peak-pos ../CD_pos.xlsx \
--indensity-pos ../已完成数据/metaX/pos/data/pos-raw-metaboAnalystInput.csv --indensity-neg ../已完成数据/metaX/neg/data/neg-raw-metaboAnalystInput.csv \
--output ../已完成数据/2.MetaboliteIdentification/identification.raw.intensity.csv --stats ../已完成数据/2.MetaboliteIdentification/identification.stat.csv

## plot HMDB
python 2.MetaboliteIdentification/HMDB.py --input ../已完成数据/2.MetaboliteIdentification/identification.raw.intensity.csv --outdir ../已完成数据/2.MetaboliteIdentification/HMDB

## KEGG statistics

python kegg_analysis.py --input ../已完成数据/2.MetaboliteIdentification/identification.raw.intensity.csv --outdir ../已完成数据/2.MetaboliteIdentification/KEGG --compound ../database/compounds.json

Rscript 2.MetaboliteIdentification/KEGG.R -l ../已完成数据/2.MetaboliteIdentification/KEGG/kegg.level2.csv -k ../已完成数据/2.MetaboliteIdentification/KEGG/kegg.csv -o ../已完成数据/2.MetaboliteIdentification/KEGG

# 3.MetaboliteQuantification
## Feature2quant
Rscript 3.MetaboliteQuantification/Feature2quant.R --intensity ../已完成数据/metaX/neg/data/neg-norm-metaboAnalystInput.csv --pca_rds ../已完成数据/metaX/neg/data/neg-PCA.rds --sample_file ../已完成数据/metaX/neg/sampleList.txt -o ../已完成数据/3.MetaboliteQuantification/Feautre2quant/neg --annotation ../已完成数据/2.MetaboliteIdentification/identification.raw.intensity.csv --type neg
Rscript 3.MetaboliteQuantification/Feature2quant.R --intensity ../已完成数据/metaX/pos/data/pos-norm-metaboAnalystInput.csv --pca_rds ../已完成数据/metaX/pos/data/pos-PCA.rds --sample_file ../已完成数据/metaX/pos/sampleList.txt -o ../已完成数据/3.MetaboliteQuantification/Feautre2quant/pos --annotation ../已完成数据/2.MetaboliteIdentification/identification.raw.intensity.csv --type pos

