#!/usr/bin/bash -l
#SBATCH --job-name=ips_install
#SBATCH --output=/home/bio.aau.dk/zs85xk/sbatch/job_%j_%x.out
#SBATCH --error=/home/bio.aau.dk/zs85xk/sbatch/job_%j_%x.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=general
#SBATCH --cpus-per-task=20
#SBATCH --mem=10G
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aohd@bio.aau.dk

# Exit on first error and if any variables are unset
set -eu
cd /home/bio.aau.dk/zs85xk/databases
#set up your environment and variables setting which versions of things you want
#mkdir interproscan
cd interproscan
#mamba create -y -n interproscan
#mamba activate interproscan
#mamba install -y -c bioconda interproscan
#mamba list -n interproscan
major_version=5.59
minor_version=91.0
pftools_version=3.2.12

#echo downloading everything you need, which will take some time
#wget http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${major_version}-${minor_version}/interproscan-${major_version}-${minor_version}-64-bit.tar.gz.md5
#wget http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${major_version}-${minor_version}/interproscan-${major_version}-${minor_version}-64-bit.tar.gz
#wget https://github.com/sib-swiss/pftools3/archive/refs/tags/v${pftools_version}.tar.gz
#md5sum -c interproscan-${major_version}-${minor_version}-64-bit.tar.gz.md5
#tar xvzf interproscan-${major_version}-${minor_version}-64-bit.tar.gz
#tar xzxf v${pftools_version}.tar.gz

#echo replacing the data folders in the conda environment with the new ones you downloaded
#rm -rf ${HOME}/.conda/envs/interproscan/share/InterProScan/data/
#cp -r interproscan-${major_version}-${minor_version}/data ${HOME}/.conda/envs/interproscan/share/InterProScan/

echo compile pftools and move the binaries you need to the conda bin
mkdir pftools3-3.2.12/build
cd pftools3-3.2.12/build
cmake .. -DCMAKE_INSTALL_PREFIX:PATH=/opt/interproscan/prosite/ -DCMAKE_BUILD_TYPE=Release -DUSE_GRAPHICS=OFF -DUSE_PDF=OFF -DUSE_AFFINITY=OFF
make
make test
rm ${HOME}/.conda/envs/interproscan/share/InterProScan/bin/prosite/pfscanV3 ${HOME}/.conda/envs/interproscan/share/InterProScan/bin/prosite/pfsearchV3
cp src/C/pfscanV3 ${HOME}/.conda/envs/interproscan/share/InterProScan/bin/prosite/
cp src/C/pfsearchV3 ${HOME}/.conda/envs/interproscan/share/InterProScan/bin/prosite/
cd ../..

echo Fix the sfld and superfamily hmm libraries as per the instructions of @pmjklemm, with minor tweak to include path
data_dir=${HOME}/.conda/envs/interproscan/share/InterProScan
for f in ${data_dir}/data/sfld/*/sfld.hmm ${data_dir}/data/superfamily/*/hmmlib_*[0-9]; do
cp $f $f.bak;
perl -lne '$_=~s/[\s\n]+$//g;if(/^(NAME|ACC) +(.*)$/){if(exists $d{$2}){$d{$2}+=1;$_.="-$d{$2}"}else{$d{$2}=0;}}print "$_"' $f > $f.tmp; mv $f.tmp $f;
done

echo go into the interproscan.properties file located in ${HOME}/.conda/envs/interproscan/share/InterProScan and make the following changes
#binary.prosite.psscan.pl.path=ps_scan.pl --> binary.prosite.psscan.pl.path=${bin.directory}/prosite/ps_scan.pl
#binary.prosite.pfscanv3.path=pfscanV3 --> binary.prosite.pfscanv3.path=${bin.directory}/prosite/pfscanV3
#binary.prosite.pfsearchv3.path=pfsearchV3 --> binary.prosite.pfsearchv3.path=${bin.directory}/prosite/pfsearchV3
#binary.pfscanv3.path=pfscanV3 --> binary.pfscanv3.path=${bin.directory}/prosite/pfscanV3
#binary.pfsearchv3.path=pfsearchV3 --> binary.pfsearchv3.path=${bin.directory}/prosite/pfsearchV3

#and HOPEFULLY it works! I was able to run interproscan on the test data without any errors.
mamba activate interproscan
test_data=${HOME}/.conda/envs/interproscan/share/InterProScan/test_all_appl.fasta
interproscan.sh -i ${test_data} -f tsv -dp
mamba deactivate

#I was able to get database matches for every analysis except Antifam and PIRSR. I'm unsure if this is the expected output since the readthedocs site isn't clear on what we should get from the test run https://interproscan-docs.readthedocs.io/en/latest/HowToRun.html#interproscan-test-run
awk '{print $4}' test_all_appl.fasta.tsv | sort | uniq -c