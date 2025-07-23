#!/bin/bash

echo "========================================"
pwd
hostname
uname -a
echo $SINGULARITY_CONTAINER
echo $SINGULARITY_NAME
echo $SINGULARITY_ENVIRONMENT
echo $SINGULARITY_BIND
shopt -s expand_aliases
source ~/.bashrc
source /etc/profile
echo "========================================"

echo ""
echo ""
echo ""
echo ""

which setup_mucoll
ls -l /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/mucoll-stack-master-h2ssl2yh2yduqnhsv2i2zcjws74v7mcq/setup.sh

echo "setup_mucoll"
setup_mucoll

echo "mkdir workdir; cd workdir"
mkdir workdir
cd workdir
pwd

echo "tar xf ../package.tar.gz"
tar xf ../package.tar.gz
ls -l

echo "cd mucstudy/"
cd mucstudy/
ls -l

echo "which python3"
which python3

sh bin/gen_lowpt_muon.sh
ls -l
