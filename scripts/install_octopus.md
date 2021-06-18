Installing requirements

I will get gcc through conda

```bash
# conda create -n octopus
conda activate octopus
conda install -c anaconda gcc_linux-64
```

It's a bit older one (7.5.0) but hopefully good enought what what I need. The rest of the requirements I will compile.


```bash
CMAKE_VERSION=3.11.4
# PYTHON_VERSION=3.6
HTSLIB_VERSION=master
OCTOPUS_VERSION=master
BOOST_VERSION=1.68.0

N_THREADS=4

DEPENDENCY_DIR="/ceph/users/kjaron/.conda/envs/octopus/dependencies"

mkdir -p $DEPENDENCY_DIR

# conda install python=$PYTHON_VERSION - gets overwritten by the next recepie
conda install -c anaconda cmake=$CMAKE_VERSION
conda install -c conda-forge boost=$BOOST_VERSION

mkdir -p /scratch/$USER/install
cd /scratch/$USER/install

echo "Installing htslib"
HTSLIB_DIR=htslib
git clone -b $HTSLIB_VERSION https://github.com/samtools/htslib.git
cd $HTSLIB_DIR

conda install -c conda-forge autoconf

autoheader
autoconf
./configure --prefix=$DEPENDENCY_DIR

conda install -c anaconda make
conda install -c conda-forge binutils # this upadted gcc, this might cause some mess up with boost, not sure

make -j$N_THREADS && make install
cd ..
rm -rf $HTSLIB_DIR

INCLUDE_PATH=$DEPENDENCY_DIR/include
export PATH=$PATH:$INCLUDE_PATH
export HTSLIB_ROOT=$DEPENDENCY_DIR

echo "Installing Octopus"
git clone -b $OCTOPUS_VERSION https://github.com/luntergroup/octopus.git
cd octopus

conda install -c conda-forge distro
conda install -c conda-forge cxx-compiler # this again upgraded GCC to 9.something
conda install -c anaconda gmp

./scripts/install.py --boost $DEPENDENCY_DIR --threads $N_THREADS

cd ..
```
