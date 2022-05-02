#extensions:
CSV_EXT=csv
LST_EXT=lst
TXT_EXT=txt

#directories:
OUT_DIR=outputs
INP_DIR=inputs
BIN_DIR=bin
ASS_DIR=${INP_DIR}/Assembly
PRO_DIR=${INP_DIR}/BioProject
SAM_DIR=${INP_DIR}/BioSample

#scripts:
get_ids=get_uids
get_header=header_csv
get_summary=fill_csv

UIDS_FILE=${INP_DIR}/ids.csv
#ASSEMBLY=${ASS_DIR}/*.${CSV_EXT}
ASSEMBLY_IDS = $(shell cat ${UIDS_FILE} | tail -n +4 | awk -F ',' '{print $$1}')
#BIOSAMPLE=${SAM_DIR}/*.${CSV_EXT}
BIOSAMPLE_IDS = $(shell cat ${UIDS_FILE} | tail -n +4 | awk -F ',' '{print $$2}')
#BIOPROJECT= ${PRO_DIR}/*.${CSV_EXT}
BIOPROJECT_IDS = $(shell cat ${UIDS_FILE} | tail -n +4 | awk -F ',' '{print $$3}')

PROTEIN_FILES = $(wildcard ${INP_DIR}/*protein.faa.gz)
BLASTP_FILES= ${PROTEIN_FILES:${INP_DIR}/%protein.faa.gz=${OUT_DIR}/%protein_blastp.tab}

NUC_FILES = $(wildcard ${INP_DIR}/*genomic.fna.gz)
BLASTN_FILES = ${NUC_FILES:${INP_DIR}/%genomic.fna.gz=${OUT_DIR}/%genomic_blastn.tab}
BED_FILES = ${NUC_FILES:${INP_DIR}/%genomic.fna.gz=${OUT_DIR}/%genomic_intersect.bed}
BED_FASTA = ${NUC_FILES:${INP_DIR}/%genomic.fna.gz=${OUT_DIR}/%genomic_intersect.fasta}



.PHONY: all display

VAR=PATH

display:
	echo ${VAR}=${${VAR}}


matrixn.${CSV_EXT}: ${BLASTN_FILES}
	ls outputs/*blastn.tab | head -n 1 | xargs -I {} cat {} | grep -oP 'Query: \K\S+' | tr '\n' ',' | head -c -1 > $@; \
	echo '' >> $@; \
	for BLASTN in outputs/*blastn.tab; \
	do \
		cat $${BLASTN} | grep -oP '^# \K\d+ (?=hits)' | tr '\n' ',' | tr -d '[:space:]' | head -c -1 >> $@; \
		echo '' >> $@; \
	done; \

matrixp.${CSV_EXT}: ${BLASTP_FILES}
	ls outputs/*blastp.tab | head -n 1 | xargs -I {} cat {} | grep -oP 'Query: \K\S+' | tr '\n' ',' | head -c -1 > $@; \
	echo '' >> $@; \
	for BLASTP in outputs/*blastp.tab; \
	do \
		cat $${BLASTP} | grep -oP '^# \K\d+ (?=hits)' | tr '\n' ',' | tr -d '[:space:]' | head -c -1 >> $@; \
		echo '' >> $@; \
	done; \

getfasta: ${BED_FASTA}

${OUT_DIR}/%genomic_intersect.fasta:
	gunzip -c ${INP_DIR}/$*genomic.fna.gz > ${INP_DIR}/$*genomic.tmp.fna; \
	bedtools getfasta -name -fi ${INP_DIR}/$*genomic.tmp.fna -bed ${OUT_DIR}/$*genomic_intersect.bed -fo $@; \
	rm ${INP_DIR}/$*genomic.tmp.fna; \
	rm ${INP_DIR}/$*genomic.fna.fai;

intersect: ${BED_FILES}

${OUT_DIR}/%genomic_intersect.bed:
	- grep -v '^#' outputs/$*genomic_blastn.tab | \
	awk '{OFS="\t"; if ($$9<=$$10) {print $$2,$$9,$$10,$$1; }}' | \
	bedtools intersect -a stdin -b inputs/$*genomic.gff.gz >> $@;

blastp: ${BLASTP_FILES}
	
${OUT_DIR}/%protein_blastp.tab: ${INP_DIR}/VFDB_setB_pro.fas 
	gunzip -c ${INP_DIR}/$*protein.faa.gz | \
	blastp -num_threads 16 -qcov_hsp_perc 0.8 -evalue 0.001 -query $< -subject - -outfmt 7 -out $@;

blastn: ${BLASTN_FILES}
	
${OUT_DIR}/%genomic_blastn.tab: ${INP_DIR}/VFDB_setB_nt.fas 
	gunzip -c ${INP_DIR}/$*genomic.fna.gz | \
	blastn -qcov_hsp_perc 0.8 -evalue 0.001 -query $< -subject - -outfmt 7 -out $@;


get_data: ${INP_DIR}/Assembly.${CSV_EXT}
	- for extension in  _genomic.gff.gz _genomic.gtf.gz _genomic.fna.gz _protein.faa.gz _translated_cds.faa.gz; \
	do \
		cat ${INP_DIR}/Assembly.${CSV_EXT} | \
		tail -n +3 | \
		awk -F, "match(\$$33, /\/GCA_/) {print \$$33 substr(\$$33, RSTART) \"$${extension}\"}" \
		| xargs -I {} wget -P ${INP_DIR} -nc {};\
		sleep 1; \
	done;
	
#VFDB_setB_nt.fas.gz VFDB_setB_pro.fas.gz
${INP_DIR}/VFDB_setB_%.fas:
	wget http://mgc.ac.cn/VFs/Down/VFDB_setB_$*.fas.gz -O VFDB_setB_$*.fas.gz; \
	gunzip -c VFDB_setB_$*.fas.gz > $@;



#Usage: make inputs/ids.csv SQUERY='40324 AND latest[filter]'
${INP_DIR}/ids.${CSV_EXT}: ${BIN_DIR}/${get_ids}
	echo "# `date`"> $@
	echo "# Assembly_search_query=${SQUERY}" >> $@
	./$< ${SQUERY} >> $@
	
	
${INP_DIR}/Summary.${CSV_EXT}: ${INP_DIR}/Assembly.${CSV_EXT} ${INP_DIR}/BioSample.${CSV_EXT} ${INP_DIR}/BioProject.${CSV_EXT}
	echo "# `date`"> $@
	paste -d ',' $< $(word 2,$^) $(word 3,$^) >> $@


${INP_DIR}/Assembly.${CSV_EXT}: ${BIN_DIR}/${get_header} ${BIN_DIR}/${get_summary} ${INP_DIR}/ids.csv
	echo "# `date`"> $@
	cat ${UIDS_FILE} | tail -n +4 | cut -d , -f 1 | head -n 1 \
	| xargs -I {} efetch -db Assembly -id {} -format docsum \
	| ./$< >> $@
	cat ${UIDS_FILE} | tail -n +4 | cut -d , -f 1 | while read ID; do \
	efetch -db Assembly -id $${ID} -format docsum  \
	| ./$(word 2,$^) >> $@; \
	done

${INP_DIR}/BioSample.${CSV_EXT}: ${BIN_DIR}/${get_header} ${BIN_DIR}/${get_summary} ${INP_DIR}/ids.csv
	echo "# `date`"> $@
	cat ${UIDS_FILE} | tail -n +4 | cut -d , -f 2 | head -n 1 \
	| xargs -I {} efetch -db BioSample -id {} -format docsum \
	| ./$< >> $@
	cat ${UIDS_FILE} | tail -n +4 | cut -d , -f 2 | while read ID; do \
	efetch -db BioSample -id $${ID} -format docsum  \
	| ./$(word 2,$^) >> $@; \
	done;
	
	
${INP_DIR}/BioProject.${CSV_EXT}: ${BIN_DIR}/${get_header} ${BIN_DIR}/${get_summary} ${INP_DIR}/ids.csv
	echo "# `date`"> $@
	cat ${UIDS_FILE} | tail -n +4 | cut -d , -f 3 | head -n 1 \
	| xargs -I {} efetch -db BioProject -id {} -format xml \
	| ./$< >> $@
	cat ${UIDS_FILE} | tail -n +4 | cut -d , -f 3 | while read ID; do \
	efetch -db BioProject -id $${ID} -format xml  \
	| ./$(word 2,$^) >> $@; \
	done;
	
	
	
	


