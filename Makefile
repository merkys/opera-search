
#extensions:
CSV_EXT=csv
LST_EXT=lst
TXT_EXT=txt

#directories:
OUT_DIR=outputs
INP_DIR=inputs
BIN_DIR=bin

#scripts:
get_ids=get_uids
#get_header=header_csv
#get_summary=fill_csv
#get_max=get_max

UIDS_FILE=${INP_DIR}/ids.csv

ASSEMBLY_IDS = $(shell cat ${UIDS_FILE} | tail -n +4 | awk -F ',' '{print $$1}')
ASSEMBLY_FILES = ${ASSEMBLY_IDS:%=${INP_DIR}/xmls/Assembly_%.xml}

BIOSAMPLE_IDS = $(shell cat ${UIDS_FILE} | tail -n +4 | awk -F ',' '{print $$2}')
BIOSAMPLE_FILES = ${BIOSAMPLE_IDS:%=${INP_DIR}/xmls/BioSample_%.xml}

BIOPROJECT_IDS = $(shell cat ${UIDS_FILE} | tail -n +4 | awk -F ',' '{print $$3}')
BIOPROJECT_FILES = ${BIOPROJECT_IDS:%=${INP_DIR}/xmls/BioProject_%.xml}

PROTEIN_FILES = $(wildcard ${INP_DIR}/*protein.faa.gz)
BLASTP_FILES= ${PROTEIN_FILES:${INP_DIR}/%protein.faa.gz=${OUT_DIR}/%protein_blastp.tab}

NUC_FILES = $(wildcard ${INP_DIR}/*genomic.fna.gz)
FASTA_FILES = ${NUC_FILES:${INP_DIR}/%genomic.fna.gz=${INP_DIR}/genomes/unzipped_genomic_fasta/%genomic.fna}
BLASTN_FILES = ${NUC_FILES:${INP_DIR}/%genomic.fna.gz=${OUT_DIR}/%genomic_blastn.tab}
BED_FILES = ${NUC_FILES:${INP_DIR}/%genomic.fna.gz=${OUT_DIR}/%genomic_intersect.bed}
BED_FASTA = ${NUC_FILES:${INP_DIR}/%genomic.fna.gz=${OUT_DIR}/%genomic_intersect.fasta}



.PHONY: all display

VAR=PATH

display:
	echo ${VAR}=${${VAR}}

get_ready:
	mkdir 

busco: ${FASTA_FILES}}
	busco -c 30 -i inputs/genomes/unzipped_genomic_fasta  -l xanthomonadales_odb10 -o xanthomonadales2 -m genome

Ranalysis/inputs/busco_values.csv:
	echo "Genome,C,S,D,F,M,n" > $@; 
	for BUSCO in xanthomonadales/*/*.txt; \
	do \
		echo -n "$(basename $${BUSCO})	" >> $@; \
		cat $${BUSCO} | perl -nle 'print "$$1,$$2,$$3,$$4,$$5,$$6" if /C:(.*)%\[S:(.*)%,D:(.*)%],F:(.*)%,M:(.*)%,n:(.*)/' >> $@; \
	done;

#extract_fasta: ${FASTA_FILES}
${INP_DIR}/genomes/unzipped_genomic_fasta/%genomic.fna:
	gunzip -c inputs/$*genomic.fna.gz > $@;

Ranalysis/inputs/matrixn.${CSV_EXT}: ${BLASTN_FILES}
	echo  -n ',' > $@;
	ls outputs/*blastn.tab | head -n 1 | xargs -I {} cat {} | grep -oP 'Query: \K\S+' | tr '\n' ',' | head -c -1 >> $@; \
	#echo '' >> $@; \
	for BLASTN in outputs/*blastn.tab; \
	do \
		echo -n "$${BLASTN}," >> $@; \
		cat $${BLASTN} | grep -oP '^# \K\d+ (?=hits)' | tr '\n' ',' | tr -d '[:space:]' | head -c -1 >> $@; \
		echo '' >> $@; \
	done; \

Ranalysis/inputs/matrixp.${CSV_EXT}: ${BLASTP_FILES}
	echo -n ',' > $@
	ls outputs/*blastp.tab | head -n 1 | xargs -I {} cat {} | grep -oP 'Query: \K\S+' | tr '\n' ',' | head -c -1 >> $@; \
	for BLASTP in outputs/*blastp.tab; \
	do \
		echo -n "$${BLASTP}," >> $@; \
		cat $${BLASTP} | grep -oP '^# \K\d+ (?=hits)' | tr '\n' ',' | tr -d '[:space:]' | head -c -1 >> $@; \
		echo '' >> $@; \
	done; \

getfasta: ${BED_FASTA}

${OUT_DIR}/%genomic_intersect.fasta: ${FASTA_FILES}
	#gunzip -c ${INP_DIR}/$*genomic.fna.gz > ${INP_DIR}/$*genomic.tmp.fna; \
	bedtools getfasta -name -fi ${INP_DIR}/genomes/unzipped_genomic_fasta/$*genomic.fna -bed ${OUT_DIR}/$*genomic_intersect.bed -fo $@; \
	rm ${INP_DIR}/genomes/unzipped_genomic_fasta/$*genomic.fna.fai;

intersect: ${BED_FILES}

${OUT_DIR}/%genomic_intersect.bed:
	- grep -v '^#' outputs/$*genomic_blastn.tab | \
	awk '{OFS="\t"; if ($$9<=$$10) {print $$2,$$9,$$10,$$1; }}' | \
	bedtools intersect -a stdin -b inputs/genomes/$*genomic.gff.gz >> $@;

blastp: ${BLASTP_FILES}
	
${OUT_DIR}/%protein_blastp.tab: ${INP_DIR}/objectives/proteins/*
	gunzip -c ${INP_DIR}/$*protein.faa.gz | \
	blastp -num_threads 16 -qcov_hsp_perc 0.8 -evalue 0.001 -query $< -subject - -outfmt 7 -out $@;

blastn: ${BLASTN_FILES}
	
${OUT_DIR}/%genomic_blastn.tab: ${INP_DIR}/objectives/nucleotides/*
	gunzip -c ${INP_DIR}/$*genomic.fna.gz | \
	blastn -qcov_hsp_perc 0.8 -evalue 0.001 -query $< -subject - -outfmt 7 -out $@;


#get_data: ${INP_DIR}/Assembly.${CSV_EXT}
#	- for extension in  _genomic.gff.gz _genomic.gtf.gz _genomic.fna.gz _protein.faa.gz _translated_cds.faa.gz; \
#	do \
#		cat ${INP_DIR}/Assembly.${CSV_EXT} | \
#		tail -n +3 | \
#		awk -F, "match(\$$33, /\/GCA_/) {print \$$33 substr(\$$33, RSTART) \"$${extension}\"}" \
#		| xargs -I {} wget -P ${INP_DIR} -nc {};\
#		sleep 1; \
#	done;
	
#VFDB_setB_nt.fas.gz VFDB_setB_pro.fas.gz
download_vf:
	wget http://mgc.ac.cn/VFs/Down/VFDB_setB_nt.fas.gz -O VFDB_setB_nt.fas.gz; \
	gunzip -c VFDB_setB_nt.fas.gz > inputs/objectives/nucleotides/VFDB_setB_nt.fas; \

	wget http://mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz -O VFDB_setB_pro.fas.gz; \
	gunzip -c VFDB_setB_pro.fas.gz > inputs/objectives/proteins/VFDB_setB_pro.fas;




#Usage: make inputs/ids.csv SQUERY='40324[Taxonomy ID] AND latest[filter]'
${INP_DIR}/ids.${CSV_EXT}: ${BIN_DIR}/${get_ids}
	echo "# `date`"> $@
	echo "# Assembly_search_query=${SQUERY}" >> $@
	./$< ${SQUERY} >> $@
	
	
Ranalysis/inputs/Metadata.${CSV_EXT}: ${UIDS_FILE}
	echo "# `date`">$@; \
	echo "Genbank	Id	RsUid	GbUid	AssemblyAccession	LastMajorReleaseAccession	ChainId	AssemblyName	Taxid	Organism	SpeciesTaxid	SpeciesName	AssemblyType	AssemblyStatus	Isolate	Sub_type	Sub_value	Coverage	ContigN50	ScaffoldN50	Title	Host	Strain	Isolation_source	Collection_date	Geo_loc_name" >> $@; \
	cat ${UIDS_FILE} | tail -n +4 | while read ROW; do \
		cat inputs/xmls/Assembly_`echo $${ROW} | cut -d "," -f 1`.xml | xtract -pattern DocumentSummary -def "null" -tab "\t" -element \
			Genbank Id RsUid GbUid AssemblyAccession LastMajorReleaseAccession ChainId  AssemblyName Taxid Organism SpeciesTaxid SpeciesName AssemblyType AssemblyStatus Isolate Sub_type Sub_value \
			Coverage ContigN50 ScaffoldN50 | tr "\n" "\t" >> $@; \
		cat inputs/xmls/BioSample_`echo $${ROW} | cut -d "," -f 2`.xml | xtract -pattern DocumentSummary -def "null" -tab "\t" -block SampleData -def "null" -tab "\t" -element Title | tr "\n" "\t" >> $@; \
		cat inputs/xmls/BioSample_`echo $${ROW} | cut -d "," -f 2`.xml | xtract -pattern DocumentSummary \
			-block Attribute -if @attribute_name -equals host  -element Attribute | tr -d "\n" >> $@; \
		echo -n "\t">> $@; \
		cat inputs/xmls/BioSample_`echo $${ROW} | cut -d "," -f 2`.xml | xtract -pattern DocumentSummary \
			-block Attribute -if @attribute_name -equals strain -element Attribute | tr -d "\n" >> $@; \
		echo -n "\t">> $@; \
		cat inputs/xmls/BioSample_`echo $${ROW} | cut -d "," -f 2`.xml | xtract -pattern DocumentSummary \
			-block Attribute -if @attribute_name -equals isolation_source -element Attribute | tr -d "\n" >> $@; \
		echo -n "\t">> $@; \
		cat inputs/xmls/BioSample_`echo $${ROW} | cut -d "," -f 2`.xml | xtract -pattern DocumentSummary \
			-block Attribute -if @attribute_name -equals collection_date -element Attribute | tr -d "\n" >> $@; \
		echo -n "\t">> $@; \
		cat inputs/xmls/BioSample_`echo $${ROW} | cut -d "," -f 2`.xml | xtract -pattern DocumentSummary \
			-block Attribute -if @attribute_name -equals geo_loc_name -element Attribute | tr -d "\n" >> $@;  \
		echo -n "\n" >> $@; \
	done;


#${INP_DIR}/Summary.${CSV_EXT}: ${INP_DIR}/Assembly.${CSV_EXT} ${INP_DIR}/BioSample.${CSV_EXT} ${INP_DIR}/BioProject.${CSV_EXT}
#	echo "# `date`"> $@
#	paste -d ',' $< $(word 2,$^) $(word 3,$^) >> $@


#${INP_DIR}/Assembly.${CSV_EXT}: ${BIN_DIR}/${get_max} ${BIN_DIR}/${get_header} ${INP_DIR}/ids.csv
#	#ROWNUM=`cat ${UIDS_FILE} | tail -n +4 | cut -d , -f 1 | while read ID; do \
	#efetch -db Assembly -id $${ID} -format docsum \
	#| ./$<; \
	#done \
	#| awk '(NR==1){Max=$$1;MaxNr=1};(NR>=2){if(Max<$$1) {Max=$$1; MaxNr=NR}} END {print MaxNr}'`; \
	#echo "hi" > $${ROWNUM};
#	ROWNUM=987; \
#	echo "# `date`"> $@
#	cat ${UIDS_FILE} | head -n 987 | tail -n 1 | cut -d , -f 1 \
#	| xargs -I {} efetch -db Assembly -id {} -format docsum \
#	| ./$(word 2,$^) >> $@
	#cat ${UIDS_FILE} | tail -n +4 | cut -d , -f 1 | while read ID; do \
	#efetch -db Assembly -id $${ID} -format docsum  \
	#| ./$(word 3,$^) >> $@; \
	#done

#${INP_DIR}/BioSample.${CSV_EXT}: ${BIN_DIR}/${get_header} ${BIN_DIR}/${get_summary} ${INP_DIR}/ids.csv
#	echo "# `date`"> $@
#	cat ${UIDS_FILE} | tail -n +4 | cut -d , -f 2 | head -n 1 \
#	| xargs -I {} efetch -db BioSample -id {} -format docsum \
#	| ./$< >> $@
#	cat ${UIDS_FILE} | tail -n +4 | cut -d , -f 2 | while read ID; do \
#	efetch -db BioSample -id $${ID} -format docsum  \
#	| ./$(word 2,$^) >> $@; \
#	done;
	
#${INP_DIR}/BioProject.${CSV_EXT}: ${BIN_DIR}/${get_header} ${BIN_DIR}/${get_summary} ${INP_DIR}/ids.csv
#	echo "# `date`"> $@
#	cat ${UIDS_FILE} | tail -n +4 | cut -d , -f 3 | head -n 1 \
#	| xargs -I {} efetch -db BioProject -id {} -format xml \
#	| ./$< >> $@
#	cat ${UIDS_FILE} | tail -n +4 | cut -d , -f 3 | while read ID; do \
#	efetch -db BioProject -id $${ID} -format xml  \
#	| ./$(word 2,$^) >> $@; \
#	done;
	
get_xmls: ${ASSEMBLY_FILES} ${BIOSAMPLE_FILES} ${BIOPROJECT_FILES}

${INP_DIR}/xmls/Assembly_%.xml:
	echo "# `date`"> $@
	efetch -db Assembly -id $* -format docsum >> $@;

${INP_DIR}/xmls/BioSample_%.xml:
	echo "# `date`"> $@
	efetch -db BioSample -id $* -format docsum >> $@;

${INP_DIR}/xmls/BioProject_%.xml:
	echo "# `date`"> $@
	efetch -db BioProject -id $* -format xml >> $@;

Ranalysis/inputs/VFmap.csv:
	cat inputs/objectives/proteins/VFDB_setB_pro.fas | grep -oP "VFG\d+" > VFG.tmp; \
	cat inputs/objectives/proteins/VFDB_setB_pro.fas | grep -oP "VF\d+" > VF.tmp; \
	cat inputs/objectives/proteins/VFDB_setB_pro.fas | grep -oP "VFC\d+" > VFC.tmp; \
 	echo "VFGID,VFID,VFCID" > $@; \
	paste -d ',' VFG.tmp VF.tmp VFC.tmp | sort | uniq >> $@; \
	rm V*.tmp;
