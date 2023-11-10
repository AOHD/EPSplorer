#!/usr/bin/env bash 
start_time=$(date +%s)
#WD="/user_data/ahd/EPS_PIPELINE/EPSplorer"
#threads=20

while getopts ":w:t:h" opt; do
  case ${opt} in
    w)
      WD="$OPTARG"
      ;;
    t)
      threads="$OPTARG"
      ;;
    h)
      echo "Usage: $(basename $0) -w <working_directory> [-t <threads>] [-h]"
      exit 0
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

shift $((OPTIND-1))

# Check if required options are provided
if [ -z "$WD" ]; then
  echo "Error: Missing required option -w <working directory>"
  exit 1
fi

# Default output directory if not provided
if [ -z "$threads" ]; then
  output_dir=20
fi

ERROR_LOG=$WD/"error.log"
mkdir $WD/"temp"
TMPDIR=$WD/"temp"
exec 2>"$ERROR_LOG"
cd $WD
mkdir $WD/"data"
mkdir $WD/"data/prokka"
mkdir $WD/"data/databases"
mkdir $WD/"data/interproscan_results"
mkdir $WD/"data/output_proximity_filtration"
mkdir $WD/"data/psiblast_results"
mkdir $WD/"figures"
mkdir $WD/"figures/operons"

##Command which runs prokka --kingdom Bacteria on all .fasta files of $WD/prokka/
echo Annotating genomes with Prokka
conda activate prokka_env
for file in $WD/genomes/*.fasta; do prokka --kingdom Bacteria --cpus $threads --outdir $WD/data/prokka/$(basename $file .fasta)/ --prefix $(basename $file .fasta) $file; done
conda deactivate
##Command which takes all the .faa files in the prokka output folder and changes the > + 8 letter string before each Prokka ID to the name of the .fasta file
for file in $WD/genomes/*.fasta; do sed -i "s/>.*_/>$(basename $file .fasta)_/g" $WD/data/prokka/$(basename $file .fasta)/*.faa; done

##Command which takes all the .gff files in the prokka output folder and changes the ID= + 8 letter string before each Prokka ID to the name of the .fasta file
for file in $WD/genomes/*.fasta; do sed -i "s/ID=.*_/ID=$(basename $file .fasta)_/g" $WD/data/prokka/$(basename $file .fasta)/*.gff; done

##Makeblastdb command which makes a blast database for each .fasta file in $WD/data/prokka/, 
##based on the .faa file in the prokka output folder. Each database should be in its own folder in $WD/databases/
echo Making blast databases
conda activate blast_env
for file in $WD/genomes/*.fasta; do makeblastdb -in $WD/data/prokka/$(basename $file .fasta)/*.faa -dbtype prot -out $WD/data/databases/$(basename $file .fasta)/$(basename $file .fasta); done

#Execute psiblast.sh
echo Running PSI-BLAST
data="$WD/data/databases/*"
query="$WD/queries_psiblast"

##Query path definitions, non-MSA queries
operon_fasta=(
acetan.faa
alginate.faa
amylovoran.faa
cellulose.faa
ColA.faa
curdlan.faa
diutan.faa
gellan1.faa
gellan2.faa
HA_Pasteurella.faa
HA_streptococcus.faa
NulO_merged.faa
pel_merged.faa
B_subtilis_EPS.faa
pnag_ica.faa
pnag_pga.faa
psl.faa
rhizobium_eps.faa
s88.faa
salecan.faa
stewartan.faa
succinoglycan.faa
vps.faa
xanthan.faa
burkholderia_eps.faa
levan.faa
synechan.faa
methanolan.faa
galactoglucan.faa
B_fragilis_PS_A.faa
B_fragilis_PS_B.faa
B_pseudomallei_EPS.faa
cepacian.faa
E_faecalis_PS.faa
emulsan.faa
EPS273.faa
GG.faa
glucorhamnan.faa
L_johnsonii_ATCC_33200_EPS_A.faa
L_johnsonii_ATCC_11506_EPS_B.faa
L_johnsonii_ATCC_2767_EPS_C.faa
L_lactis_EPS.faa
L_plantarum_HePS.faa
phosphonoglycan.faa
)
## PSI-BLAST of all operons against all databases in $WD/queries_psiblast/databases/
for database in ${data[@]}; do
mkdir $WD/data/psiblast_results/$(basename $database)
for operon in ${operon_fasta[@]}; do
psiblast -query $query/$operon -db $database/$(basename $database) -out $WD/data/psiblast_results/$(basename $database)/$operon -evalue 0.0001 -qcov_hsp_perc 20 -max_hsps 10 -max_target_seqs 100000 -outfmt 6 -num_iterations 20 -comp_based_stats 1 -num_threads $threads
echo $operon BLASTed
done
echo $database completed
done
conda deactivate

##Move all .faa files in $WD/data/prokka/*/* to $WD/data/prokka/
for file in $WD/data/prokka/*/*.faa; do
cp $file $WD/data/prokka
echo $file
done


##Concatenate all .faa files with the same name in $WD/data/psiblast_results/*/*.faa
echo Creating concatenated GFF file
mkdir $WD/data/psiblast_results/concatenated
for file in $WD/data/psiblast_results/*/*.faa; do cat $file >> $WD/data/psiblast_results/concatenated/$(basename $file .faa).faa; done

##Run $WD/scripts/generate_gff.R
conda activate R_env
Rscript $WD/scripts/generate_gff.R "$WD"


##Run $WD/scripts/proximity_main.R
echo Running proximity filtration
Rscript $WD/scripts/proximity_main.R "$WD"

##InterProScan analysis
echo Running InterProScan
# Folder with all subsetted fasta files
FASTA=$WD/data/output_proximity_filtration/fasta_output
# Results folder
RESULTS=$WD/data/interproscan_results

# Loading interproscan

for j in $FASTA/*; do
# For each polysaccharide
j=`basename $j`
echo $j
mkdir -p $RESULTS/$j
mkdir -p $RESULTS/$j\_fasta_removed
for i in $FASTA/$j/*; do
	# For each MAG id
	# Removing .faa from the name
	i=`basename $i | rev | cut -c 5- | rev`
	echo $i
	if [[ ! -s $RESULTS/$j\_fasta_removed/$i.gff3 ]]; then
	# Running interproscan
	## Only Pfam database
	## Only GFF3 file is written out
	interproscan.sh -cpu $threads -appl Pfam -f GFF3 -o $RESULTS/$j/$i.gff3 -i $FASTA/$j/$i.faa
	
	# Remove trailing fasta sequence in gff3 file and leading three lines
	F=$RESULTS/$j/$i.gff3
	sed '/^##FASTA$/,$d' ${F} | sed '1,/^##interproscan*/d' > $RESULTS/$j\_fasta_removed/$i.gff3
	fi
done
done


##Run $WD/scripts/ips_main.R
echo Generating gene arrow plots
Rscript $WD/scripts/ips_main.R "$WD"

echo Generating overview excel file
Rscript $WD/scripts/overview.R "$WD"


exec 2>&1

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

module purge
conda deactivate
echo Finished in $elapsed_time seconds! Have a nice day!


