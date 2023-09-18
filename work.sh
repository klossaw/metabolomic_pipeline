THREAD=8

## get samples
sampleInfoFile="sample_info.xls"
sampleListTMP=`cut -f 2 ${sampleInfoFile} | sed 1d | sed ':a;N;s/\n/", "/;ta;'`
sampleList="\"${sampleListTMP}\""

## get project running directory
projectDir=`pwd`

## config snakemake file
programDir="/home/project/pipelines/Metabolomics/Metabolomics_pipe_V1.3"
pipeFileName="Metabolomics_Snakefile.py"
cp ${programDir}/${pipeFileName} ${projectDir}/${pipeFileName}
sed -i "1iProgramDir = \"${programDir}\"" ${pipeFileName}
sed -i "1iTHREAD = ${THREAD}" ${pipeFileName}
sed -i "1iSAMPLES = [${sampleList}]" ${pipeFileName}
sed -i "1iPROJECT = \"${projectDir}\"" ${pipeFileName}

## run the pipeline
rawDir="raw"
perl ${programDir}/bin/prepare_data.pl ${sampleInfoFile} ${rawDir} 
conda activate /opt/anaconda3/envs/Metabolomics/
/home/project/.local/bin/snakemake -s ${pipeFileName} --dag | dot -Tpdf > dag.pdf
nohup /home/project/.local/bin/snakemake -s ${pipeFileName} -p -j ${THREAD} >> run.log 2>&1 &
