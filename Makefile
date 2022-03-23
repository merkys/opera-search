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

UIDS_FILE=${INP_DIR}/*.${TXT_EXT}
#ASSEMBLY=${ASS_DIR}/*.${CSV_EXT}
ASSEMBLY_IDS = $(shell cat ${UIDS_FILE} | tail -n +4 | awk -F ',' '{print $$1}')
#BIOSAMPLE=${SAM_DIR}/*.${CSV_EXT}
BIOSAMPLE_IDS = $(shell cat ${UIDS_FILE} | tail -n +4 | awk -F ',' '{print $$2}')
#BIOPROJECT= ${PRO_DIR}/*.${CSV_EXT}
BIOPROJECT_IDS = $(shell cat ${UIDS_FILE} | tail -n +4 | awk -F ',' '{print $$3}')




.PHONY: all display

VAR=PATH

display:
	echo ${VAR}=${${VAR}}

#Usage: make get_ids SQUERY=pubmed_search_query
get_ids: ${BIN_DIR}/${get_ids}
	FNAME=`echo ${SQUERY} | sed -e 's/[]\$\!\*\?\;\&\>\<)\(}\{\[]\t /_/g'`; \
	echo "# `date`"> ${INP_DIR}/$${FNAME}`date +%Y-%m-%d`.${TXT_EXT}; \
	echo "# Assembly_search_query=${SQUERY}" >> ${INP_DIR}/$${FNAME}`date +%Y-%m-%d`.${TXT_EXT}; \
	./$< ${SQUERY} >> ${INP_DIR}/$${FNAME}`date +%Y-%m-%d`.${TXT_EXT};
	


${INP_DIR}/Assembly.${CSV_EXT}: ${BIN_DIR}/${get_header} ${BIN_DIR}/${get_summary} 
	echo "# `date`"> $@; \
	echo ${ASSEMBLY_IDS} \
	| grep -o "^\w*\b" \
	| xargs -I {} efetch -db Assembly -id {} -format docsum \
	| ./$< >> $@; \
	for ID in ${ASSEMBLY_IDS}; do \
	efetch -db Assembly -id $${ID} -format docsum  \
	| ./$(word 2,$^) >> $@; \
	done;

${INP_DIR}/BioSample.${CSV_EXT}: ${BIN_DIR}/${get_header} ${BIN_DIR}/${get_summary}
	echo "# `date`"> $@; \
	echo ${BIOSAMPLE_IDS} \
	| grep -o "^\w*\b" \
	| xargs -I {} efetch -db BioSample -id {} -format docsum \
	| ./$< >> $@; \
	for ID in ${BIOSAMPLE_IDS}; do \
	efetch -db BioSample -id $${ID} -format docsum  \
	| ./$(word 2,$^) >> $@; \
	done;
	
	
${INP_DIR}/BioProject.${CSV_EXT}: ${BIN_DIR}/${get_header} ${BIN_DIR}/${get_summary} 
	echo "# `date`"> $@; \
	echo ${BIOPROJECT_IDS} \
	| grep -o "^\w*\b" \
	| xargs -I {} efetch -db BioProject -id {} -format xml \
	| ./$< >> $@; \
	for ID in ${BIOPROJECT_IDS}; do \
	efetch -db BioProject -id $${ID} -format xml  \
	| ./$(word 2,$^) >> $@; \
	done;
	
	
	
	


