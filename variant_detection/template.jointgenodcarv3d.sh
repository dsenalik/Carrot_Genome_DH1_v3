#!/bin/bash
#SBATCH -p msn                #name of the queue you are submitting job to
#SBATCH -q msn                # qos use msn for msn notes, memlimit for public short etc.
#SBATCH -N 1                  #number of nodes in this job
#SBATCH -n TTTT                  #number of cores/tasks in this job; there are 20 physical cores on a node * 2 = 40 logical cores
#SBATCH -t 2-00:00:00         #time allocated for this job days-hours:mins:seconds 168h=7d
#SBATCH --mem-per-cpu 6GB     #default is 3.2GB
#SBATCH --mail-user=dsenalik@gmail.com   #enter your email address to receive emails
#SBATCH --mail-type=FAIL      #BEGIN,END,FAIL will receive an email when job starts, ends or fails
#SBATCH -o "%j.%N.stdout"     #standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "%j.%N.stderr"     #optional but it prints our standard error
#SBATCH --job-name=jobKKKK

echo "+++ $0 Starting  $(date)"
touch "flag.started"
set -o pipefail
set -u # error if reference uninitialized variable
shopt -s nullglob
set -e



# values here are replaced here by the makejob script
job="KKKK"
prefix="PPPP"
region="RRRR"
timeout=3600
regionnocolon=$(echo "$region" | tr ":" "_")
echo "+++ Job=\"$job\" prefix=\"$prefix\" region=\"$region\""



# load modules
#module load xxx



# configuration
# $TMPDIR points to a private per-job directory on the direct-attached SSD storage on the compute
#   node, rather than in the small /tmp on the login node
# use mounted directory for project directory, real path is not available from within the container
carrot_genome="/scinet01/gov/usda/ars/scinet/project/carrot_genome"
container_carrot_genome="/carrot_genome"
container_reference="$container_carrot_genome/ref/130.DHv2.fna"
outdir="$carrot_genome/job$job"
container_outdir="$container_carrot_genome/job$job"
gatkjointgenotypememory=24
db="${prefix}_${region}.genomicsdb"
container_db="$container_outdir/$db"
zip="${prefix}_${regionnocolon}.genomicsdb.zip"
vcf="${prefix}_${regionnocolon}.vcf.gz"
container_vcf="$container_outdir/${prefix}_${regionnocolon}.vcf.gz"
log="${prefix}_${regionnocolon}.jointgenotype.log"
md5="${prefix}_${regionnocolon}.md5"



# initialization
if [ ! -d "$outdir" ]; then
  echo "Error, directory \"$outdir\" does not exist"
  exit 1
fi
cd "$outdir"
echo "+++ Current directory is \"$(pwd)\""
echo "HOSTNAME = $HOSTNAME"
echo "TMPDIR = $TMPDIR"



# download genomicsdb zip file from our NAS server using FTP
if [ ! -s "$zip" ]; then
  echo "+++ Downloading genomicsdb  $(date)"
  echo '#!/usr/bin/expect
set timeout '$timeout'
spawn sftp ceres@redacted
expect "ceres@redacted'\''s password: "
send "redacted\n"
expect "sftp>"
send "get vcf/0079.13/PPPP/'$zip'\n"
expect "sftp>"
send "exit\n"
expect eof
' > job$job.$$.1.expect
  expect job$job.$$.1.expect
  rm -f job$job.$$.1.expect
fi

if [ ! -s "$zip" ]; then
  echo "Did not download zip file \"$zip\" to $(pwd)"
  exit 1
else
  echo "+++ Uncompressing genomicsdb  $(date)"
  unzip "$zip"
  r=$?; if [ "$r" != "0" ]; then echo "Uncompress failure r=$r"; exit $r; fi
fi
# remove compressed version to recover some disk space
rm -f "$zip"

echo "+++ Launching the GATK  $(date)"
touch "flag.running"
singularity exec \
  -B /scinet01/gov/usda/ars/scinet/project/carrot_genome:/carrot_genome \
  /scinet01/gov/usda/ars/scinet/project/carrot_genome/singularity/gatk-4.0.7.0.simg \
  gatk GenotypeGVCFs \
  --java-options "-Xmx${gatkjointgenotypememory}G" \
  --reference "$container_reference" \
  --variant "gendb://$container_db" \
  --use-new-qual-calculator \
  --only-output-calls-starting-in-intervals \
  --intervals "$region" \
  --output "$container_vcf" \
  &> "$log"
r=$?; if [ "$r" != "0" ]; then echo "GATK failure r=$r"; exit $r; fi

echo "+++ The GATK run has completed  $(date)"
gzip -9 "$log"
md5sum "$vcf"* > $md5
#  bgzip "$outdir/$uid.gvcf" && tabix -p vcf "$outdir/$uid.gvcf.gz"


echo "+++ Directory listing prior to upload:"
ls -ltr
echo

# upload vcf file and log to our NAS server using FTP
echo '#!/usr/bin/expect
set timeout '$timeout'
spawn sftp ceres@redacted
expect "ceres@redacted'\''s password: "
send "redacted\n"
expect "sftp>"
send "put '$log'.gz vcf/0079.13/PPPP/\n"
expect "sftp>"
send "put '$vcf' vcf/0079.13/PPPP/\n"
expect "sftp>"
send "put '$vcf'.tbi vcf/0079.13/PPPP/\n"
expect "sftp>"
send "put '$md5' vcf/0079.13/PPPP/\n"
expect "sftp>"
send "exit\n"
expect eof
' > job$job.$$.2.expect
expect job$job.$$.2.expect
rm -f job$job.$$.2.expect



# cleanup, delete the genomicsdb and uploaded results
echo "+++ Removing genomics db  $(date)"
rm -rf "$db"



touch "flag.completed"
echo "+++ $0 Done  $(date)"
exit 0

