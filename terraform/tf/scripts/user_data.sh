#! /bin/bash

export DEBIAN_FRONTEND=noninteractive

# directories, logs
mkdir /root/logs
mkdir /root/github

sudo apt-get update 
update-alternatives --remove python /usr/bin/python2 
update-alternatives --install /usr/bin/python python /usr/bin/python3 1 

# Check default python is python3 
sudo apt-get install python3-pip -y

# symlink for pip3 to pip
which pip3
ln -s /usr/bin/pip3 /usr/bin/pip 

#Miscellaneous dependencies
sudo pip install numpy biopython scipy 
sudo apt-get install python-qt4 python-lxml python-six -y
sudo pip install --upgrade ete3 

GITHUB_PATH="/root/github"

cd $GITHUB_PATH

# RAxML
git clone https://github.com/stamatak/standard-RAxML --branch v8.2.12 > /root/logs/raxml_setup.log
cd standard-RAxML; make -f Makefile.SSE3.PTHREADS.gcc; make -f Makefile.gcc; make -f Makefile.SSE3.gcc; make -f Makefile.AVX.gcc; make -f Makefile.PTHREADS.gcc; make-f Makefile.AVX.PTHREADS.gcc; rm *.o

cp raxmlHPC* /usr/bin/

cd $GITHUB_PATH

#Prodigal
git clone https://github.com/hyattpd/Prodigal --branch v2.6.3 > /root/logs/prodigal_setup.log
cd Prodigal; make install >> /root/logs/prodigal_setup.log

cd $GITHUB_PATH
#trimal
git clone https://github.com/scapella/trimal --branch v1.4.1 > /root/logs/trimal_setup.log
cd trimal/source; make >> /root/logs/trimal_setup.log; cp trimal /usr/bin/

#HMMER
apt-get install hmmer -y >> /root/logs/hmmer_setup.log

#MAFFT v7
cd /tmp
wget https://mafft.cbrc.jp/alignment/software/mafft_7.407-1_amd64.deb
dpkg -i mafft_7.407-1_amd64.deb

#papara
wget https://sco.h-its.org/exelixis/resource/download/software/papara_nt-2.5-static_x86_64.tar.gz
tar -xvf papara* --directory=/usr/bin; cd /usr/bin; mv papara* papara

cd $GITHUB_PATH

git clone --single-branch -b autoinstall https://github.com/hallamlab/TreeSAPP.git; cd TreeSAPP;
cd sub_binaries; cp usearch /usr/bin; cd ../
make; make install;

cd ~/

#Move github folder to user home directory
cp -r github /home/${username}

ALIAS="alias treesapp='python /home/${username}/github/TreeSAPP/treesapp.py'"
echo $ALIAS >> /home/${username}/.bashrc

source /home/${username}/.bashrc
