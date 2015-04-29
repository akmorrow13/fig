#!/bin/bash
#
# Copyright 2015 Frank Austin Nothaft
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

set -x

export SPARK_HOME=/root/spark
# useful other vars
export public_hostname=$(curl http://169.254.169.254/latest/meta-data/public-hostname)
export hdfs_root="hdfs://$public_hostname:9000"

# sets addl vars, incl AWS credentials
source /root/spark-ec2/ec2-variables.sh
export SPARK_MASTER="$MASTERS"

# get and install maven
wget ftp://apache.cs.utah.edu/apache.org/maven/maven-3/3.3.3/binaries/apache-maven-3.3.3-bin.tar.gz
tar xzvf apache-maven-3.3.3-bin.tar.gz

# build fig
./apache-maven-3.3.3/bin/mvn package -DskipTests

# start MR nodes
../ephemeral-hdfs/bin/stop-all.sh
sleep 10
../ephemeral-hdfs/bin/start-all.sh

# make a directory in hdfs
../ephemeral-hdfs/bin/hadoop fs -mkdir .

# pull NA12878 from 1000g
../ephemeral-hdfs/bin/hadoop distcp \
    s3n://bdg-eggo/1kg/genotypes \
    ${hdfs_root}/user/${USER}/1000kg.adam

# download bulk datafiles
cd /mnt
wget http://compbio.mit.edu/encode-motifs/matches.txt.gz
gunzip matches.txt.gz
~/ephemeral-hdfs/bin/hadoop fs -put \
    matches.txt \
    ${hdfs_root}/user/${USER}/matches.txt
rm matches.txt

wget http://compbio.mit.edu/encode-motifs/motifs.txt

# pull hg19
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit

# pull genes
wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
gunzip Homo_sapiens.GRCh37.75.gtf.gz
awk '{print "chr"$0}' Homo_sapiens.GRCh37.75.gtf | sed 's/chrMT/chrM/g' > hg19.gtf
~/ephemeral-hdfs/bin/hadoop fs -put \
    hg19.gtf \
    ${hdfs_root}/user/${USER}/hg19.gtf
rm *.gtf

cd ~/fig

# run fig
./bin/fig-submit \
    ${hdfs_root}/user/${USER}/1000kg.adam \
    /mnt/hg19.2bit \
    ${hdfs_root}/user/${USER}/hg19.gtf \
    /mnt/motifs.txt \
    ${hdfs_root}/user/${USER}/matches.txt
