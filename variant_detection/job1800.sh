#!/bin/bash
#SBATCH -p msn                #name of the queue you are submitting job to
#SBATCH -q msn                # qos use msn for msn notes, memlimit for public medium etc.
#SBATCH -N 1                  #number of nodes in this job
#SBATCH -n 4                  #number of cores/tasks in this job; there are 20 physical cores on a node * 2 = 40 logical cores
#SBATCH -t 5-00:00:00         #time allocated for this job days-hours:mins:seconds 168h=7d
#SBATCH --mem-per-cpu 6GB     #default is 3.2GB
#SBATCH --mail-user=example@zombo.com   #enter your email address to receive emails
#SBATCH --mail-type=FAIL      #BEGIN,END,FAIL will receive an email when job starts, ends or fails
#SBATCH -o "%j.%N.stdout"     #standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "%j.%N.stderr"     #optional but it prints our standard error
#SBATCH --job-name=80001

echo "+++ $0 Starting  $(date)"
touch "flag.started"
set -o pipefail
set -u # error if reference uninitialized variable
shopt -s nullglob
set -e



# job number e.g. 0007 is replaced here by the makejob script
job="1800"
echo "+++ Job = \"$job\""

uids="80001"
echo "+++ uids = \"$uids\""

ftpsource="bam/0082.04"
echo "+++ bam ftp download source = \"$ftpsource\""

ftpdest="vcf/0082.05"
echo "+++ gvcf ftp upload destination = \"$ftpdest\""



# load modules
module load tabix
module load samtools
module load picard/64



# configuration
# $TMPDIR points to a private per-job directory on the direct-attached SSD storage on the compute
#   node, rather than in the small /tmp on the login node
# use mounted directory for project directory, real path is not available from within the container
carrot_genome="/scinet01/gov/usda/ars/scinet/project/carrot_genome"
container_carrot_genome="/carrot_genome"
container_reference="$container_carrot_genome/ref/DCARv3d.fna"
outdir="$carrot_genome/job$job"
container_outdir="$container_carrot_genome/job$job"
gatkhaplotypecallermemory=24
ploidy=2



# initialization
if [ ! -d "$outdir" ]; then
  echo "Error, directory \"$outdir\" does not exist"
  exit 1
fi
cd "$outdir"
echo "+++ Current directory is \"$(pwd)\""
echo "HOSTNAME = \"$HOSTNAME\""
echo "TMPDIR = \"$TMPDIR\""
df -h "$TMPDIR"



for uid in $uids ; do

  # download bam file and index from our NAS server using FTP
  if [ ! -s "$uid.bam" ]; then
    echo "+++ Downloading bam file for \"$uid\""
    echo '#!/usr/bin/expect
set timeout 1800
spawn sftp ceres@ourserver.edu
expect "ceres@ourserver.edu'\''s password: "
send "supersecretpassword\n"
expect "sftp>"
send "get '$ftpsource'/'$uid'.ba?\n"
expect "sftp>"
send "get '$ftpsource'/'$uid'.md5\n"
expect "sftp>"
send "exit\n"
expect eof
' > job$job.$$.1.expect
    expect job$job.$$.1.expect
    rm -f job$job.$$.1.expect
  fi

  if [ ! -s "$uid.bam" ]; then
    echo "Did not download bam file \"$uid.bam\" to $(pwd)"
  else
    echo "+++ Checking bam file integrity  $(date)"
    md5sum -c "$uid.md5"
    rm -f "$uid.md5"

    echo "+++ Filtering bam file for sample $uid  $(date)"
    samtools view -b -q 30 -F 256 "$uid.bam" > "$uid.sF30.bam"
    rm -f "$uid.bam" "$uid.bai"

    echo "+++ Marking duplicates for sample $uid  $(date)"
    run_picard MarkDuplicates \
      INPUT="$uid.sF30.bam" \
      OUTPUT="$uid.sF30u.bam" \
      METRICS_FILE="$uid.metrics.txt" \
      CREATE_INDEX=true \
      TMP_DIR="$TMPDIR" \
      MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=900 \
      &> "$outdir/$uid.picard.log"
    rm -f "$uid.sF30.bam"
    gzip -9 "$outdir/$uid.picard.log"

    echo "+++ Launching the GATK for sample $uid  $(date)"
    singularity exec \
      -B /project/carrot_genome:/carrot_genome \
      /project/carrot_genome/singularity/gatk_4.1.4.0.sif \
      gatk HaplotypeCaller \
      --java-options "-Xmx${gatkhaplotypecallermemory}G" \
      -R "$container_reference" \
      -I "$container_outdir/$uid.sF30u.bam" \
      --emit-ref-confidence GVCF \
      --sample-ploidy $ploidy \
      --output "$container_outdir/$uid.gvcf" \
      &> "$outdir/$uid.haplotypecaller.log" &
  fi

done  # for uid
wait
echo "+++ All the GATK runs have completed  $(date)"
for uid in $uids ; do
  bgzip "$outdir/$uid.gvcf"
  tabix -p vcf "$outdir/$uid.gvcf.gz"
  gzip -9 "$outdir/$uid.haplotypecaller.log"
  md5sum "$outdir/$uid.gvcf.gz" > "$outdir/$uid".md5
  md5sum "$outdir/$uid.gvcf.gz.tbi" >> "$outdir/$uid".md5
done


echo "+++ Directory listing prior to upload:"
ls -ltr
echo

# upload vcf files to our NAS server using FTP
echo '#!/usr/bin/expect
set timeout 1800
spawn sftp ceres@ourserver.edu
expect "ceres@ourserver.edu'\''s password: "
send "supersecretpassword\n"
expect "sftp>"
send "put *.log.gz '$ftpdest'/logs/\n"
expect "sftp>"
send "put *.metrics.txt '$ftpdest'/logs/\n"
expect "sftp>"
send "put *.gvcf.gz* '$ftpdest'/\n"
expect "sftp>"
send "put *.tbi '$ftpdest'/\n"
expect "sftp>"
send "exit\n"
expect eof
' > job$job.$$.2.expect
expect job$job.$$.2.expect
rm -f job$job.$$.2.expect



# cleanup, delete the downloaded bam file and index
# and uploaded results
for uid in $uids ; do
  echo "+++ Removing $uid bam file"
  rm -f "$uid.sF30u.ba"? "$uid.gvcf.idx"
done



touch "flag.completed"
echo "+++ $0 Done  $(date)"
exit 0
