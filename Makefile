#extensions:
CSV_EXT=csv
LST_EXT=lst

#directories:
OUT_DIR=outputs
INP_DIR=inputs
BIN_DIR=bin

#scripts:
get_ids=get_uids
get_header=header_csv
get_summary=fill_csv

ASSEMBLY=${INP_DIR}/Assembly/*.${CSV_EXT}
ASSEMBLY_IDS = $(shell cat ${ASSEMBLY} | grep -o "^[^\#]*" | sort -u)
BIOPROJECT= ${INP_DIR}/*.${CSV_EXT}
BIOPROJECT_IDS = $(shell cat ${BIOPROJECT} | grep -o "^[^\#]*" | sort -u)
BIOSAMPLE=${INP_DIR}/Assembly/*.${CSV_EXT}
BIOSAMPLE_IDS = $(shell cat ${BIOSAMPLE} | grep -o "^[^\#]*" | sort -u)




.PHONY: all display

VAR=PATH

display:
	echo ${VAR}=${${VAR}}

#Usage:
# make get_ids SQUERY=pubmed_search_query
get_ids: ${BIN_DIR}/${get_ids}
	FNAME=`echo ${SQUERY} | sed -e 's/[]\$\!\*\?\;\&\>\<)\(}\{\[]\t /_/g'`; \
	echo "# `date`"> ${INP_DIR}/$${FNAME}`date +%Y-%m-%d`.${CSV_EXT}; \
	echo "# Assembly_search_query=${SQUERY}" >> ${INP_DIR}/$${FNAME}`date +%Y-%m-%d`.${CSV_EXT}; \
	./$< ${SQUERY} >> ${INP_DIR}/$${FNAME}`date +%Y-%m-%d`.${CSV_EXT};

get_ass get_assembly get_assembly_summary: ${BIN_DIR}/${get_header}
	tail -n +4 ${INP_DIR}/Stenotrophomonas\ maltophilia\ latest2022-03-21.csv \
	| awk -F ',' 'print $$2}' \
	| `while read CMD; do \
    	echo $${CMD} \
	done;`
	


