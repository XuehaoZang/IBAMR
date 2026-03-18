# Build IBAMR

This is a protocol for setting up environment to run **IBAMR** for aerodynamic simulation.

IBAMR, Immersed Boundary Method Adaptive Mesh Refinement Software Infrastructure, is a distributed-memory parallel implementation of the immersed boundary (IB) method with support for Cartesian grid adaptive mesh refinement (AMR). Support for distributed-memory parallelism is via [MPI](http://www.mcs.anl.gov/research/projects/mpi), the Message Passing Interface.

## References

https://github.com/IBAMR/

https://ibamr.github.io/building

## Dependencies

- System: Windows WSL2 or MacOS XCode

  - All dependencies are downloaded and unpacked in `$HOME/Code` and installed in `$HOME/Applications`

- PETSc: Portable, Extensible Toolkit for Scientific Computation

  - Used as package manager: install MPI HDF5 HYPRE

  - v3.23.3

- Boost, Eigen, muParser, HYPRE, SILO

- MPI: Message Passing Interface
  - Pass data between core processors
  - with Fortran and C compiles `mpicc`

- HDF5: storage format
- libMesh: creating mesh for structures
  - v1.7.8 with debugging + optimized build
- SAMRAI: Structured Adaptive Mesh Refinement Application Infrastructure

## Set up environment

### WSL2 Ubuntu in Windows

#### Create new WSL

``` bash
# PS C:\Users\computer0>
New-Item -ItemType Directory -Path "C:\WSL\IBAMR_Sim”

Invoke-WebRequest -Uri "https://cloud-images.ubuntu.com/releases/22.04/release/ubuntu-22.04-server-cloudimg-amd64-root.tar.xz" -OutFile "C:\WSL\ubuntu2204_rootfs.tar.xz" 

wsl --import IBAMR_Sim "C:\WSL\IBAMR_Sim" "C:\WSL\ubuntu2204_rootfs.tar.xz”

wsl -l -v 
# IBAMR_Sim  Stopped     2 

exit
```

#### Set up user account

```bash
# PS C:\Users\computer0>
wsl -d IBAMR_Sim 

# root@DESKTOP-0:/****mnt****/c/Users/computer0#
adduser sim # Username: sim; Password: CohenLab
usermod -aG sudo sim
groups sim  # sim : sim sudo
exit

# PS C:\Users\computer0>
wsl -d IBAMR_Sim 

sudo tee /etc/wsl.conf <<EOF
[user]
default=sim
EOF

cat /etc/wsl.conf # [user] default=sim
exit

# PS C:\Users\computer0>
wsl –shutdown
```

### **Activate** **conda** **env**

```bash
# PS C:\Users\computer0>
wsl -d IBAMR_Sim		
# sim@DESKTOP-0:/mnt/c/Users/computer0$

# Basics: git gcc cmake
sudo apt update && sudo apt upgrade -y
sudo apt install -y build-essential cmake git gfortran libtool wget m4

# Miniconda	 # PREFIX=/home/sim/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

~/miniconda3/bin/conda init bash 	
source ~/.bashrc	

# (base) sim@DESKTOP-0:/mnt/c/Users/computer0$
conda create -n ibamr python=3.10 cmake make –y
conda activate ibamr
```

### Install `PETSc`

```bash
# PS C:\Users\computer0>
wsl -d IBAMR_Sim		
# sim@DESKTOP-0:/mnt/c/Users/computer0$
conda activate ibamr

# (ibamr) sim@DESKTOP-0:/mnt/c/Users/computer0$
mkdir -p $HOME/Applications	# for software
mkdir -p $HOME/Code		# for code
cd $HOME/Code

# PETSc 3.23.3
wget https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-3.23.3.tar.gz
tar xf petsc-3.23.3.tar.gz
cd petsc-3.23.3
unset PETSC_DIR
unset PETSC_ARCH

# config
# Note: using -O3 -march=native for performance
./configure --force \
  --COPTFLAGS='-O3 -march=native -g' \
  --CXXOPTFLAGS='-O3 -march=native -g' \
  --FOPTFLAGS='-O3 -march=native -g' \
  --download-hdf5=1 \
  --with-shared-libraries=1 \
  --download-hypre=1 \
  --download-openmpi=1 \
  --download-openblas=1 \
  --download-metis=1 \
  --download-parmetis=1 \
  --download-silo=1 \
  --with-x=0 \
  --with-64-bit-indices=0 \
  --with-fortran-bindings=0 \
  --with-mpi=1 \
  --with-debugging=0 \
  --prefix=$HOME/Applications/petsc-3.23.3/

make PETSC_DIR=/home/sim/Code/petsc-3.23.3 PETSC_ARCH=arch-linux-c-opt all 

make PETSC_DIR=/home/sim/Code/petsc-3.23.3 PETSC_ARCH=arch-linux-c-opt install

make PETSC_DIR=/home/sim/Applications/petsc-3.23.3 PETSC_ARCH="" check
```

Installed to `/home/sim/Applications/petsc-3.23.3`

### Install **`libmesh`** (debug+opt build)

```bash
# (ibamr) sim@DESKTOP-0:/mnt/c/Users/computer0$
cd $HOME/Code
wget https://github.com/libMesh/libmesh/releases/download/v1.7.8/libmesh-1.7.8.tar.gz
tar xf libmesh-1.7.8.tar.gz
```

#### Debug build - devel

```bash
cd $HOME/Code/libmesh-1.7.8
mkdir -p build-debug
cd build-debug

# Note: add flag to MPI in PETSc: --with-mpi=$HOME/Applications/petsc-3.23.3 \
# configuration
../configure \
  --prefix=$HOME/Applications/libmesh-1.7.8-debug \
  --with-methods=devel \
  --enable-exodus \
  --enable-triangle \
  --enable-petsc-required \
  --disable-boost \
  --disable-eigen \
  --disable-hdf5 \
  --disable-metaphysicl \
  --disable-openmp \
  --disable-perflog \
  --disable-pthreads \
  --disable-tbb \
  --disable-timestamps \
  --disable-reference-counting \
  --disable-strict-lgpl \
  --disable-glibcxx-debugging \
  --disable-vtk \
  --with-thread-model=none \
  --with-mpi=$HOME/Applications/petsc-3.23.3 \
  PETSC_DIR=$HOME/Applications/petsc-3.23.3 \
  PETSC_ARCH=""

# using 14 cores in parallel
make –j14		
make –j14 install
```

#### Optimized build - opt

```bash
cd $HOME/Code/libmesh-1.7.8
mkdir -p build-opt
cd build-opt

../configure \
  --prefix=$HOME/Applications/libmesh-1.7.8-opt \
  --with-methods=opt \
  --enable-exodus \
  --enable-triangle \
  --enable-petsc-required \
  --disable-boost \
  --disable-eigen \
  --disable-hdf5 \
  --disable-metaphysicl \
  --disable-openmp \
  --disable-perflog \
  --disable-pthreads \
  --disable-tbb \
  --disable-timestamps \
  --disable-reference-counting \
  --disable-strict-lgpl \
  --disable-glibcxx-debugging \
  --disable-vtk \
  --with-thread-model=none \
  --with-mpi=$HOME/Applications/petsc-3.23.3 \
  PETSC_DIR=$HOME/Applications/petsc-3.23.3 \
  PETSC_ARCH=""
  
make -j14
make -j14 install
```

Installed to `/home/sim/Applications/libmesh-1.7.8-debug` and `/home/sim/Applications/ libmesh-1.7.8-opt`

### Install `SAMRAI` v2025.10.29

IBAMR requires a very old version of SAMRAI (version 2.4.4). Developers maintain their own fork with a set of patches which allow that version of SAMRAI to be compiled with modern compilers.

```bash
cd $HOME/Code
# New version of IBAMR requires a 2025.10.29 version of SAMRAI
wget https://github.com/IBAMR/IBSAMRAI2/archive/refs/tags/2025.10.29.tar.gz
tar xf 2025.10.29.tar.gz
```

#### Debug build

```bash
cd $HOME/Code/IBSAMRAI2-2025.10.29
mkdir -p build-debug
cd build-debug

# Notes: needs flags to identify MPI and HDF5 locations
../configure \
  CFLAGS="-fPIC -g -O2" \
  CXXFLAGS="-fPIC -g -O2" \
  FFLAGS="-fPIC -g -O2" \
  --prefix=$HOME/Applications/ibsamrai2-2025.10.29-debug \
  --with-CC="$HOME/Applications/petsc-3.23.3/bin/mpicc" \
  --with-CXX="$HOME/Applications/petsc-3.23.3/bin/mpicxx" \
  --with-F77="$HOME/Applications/petsc-3.23.3/bin/mpif90" \
  --with-hdf5="$HOME/Applications/petsc-3.23.3/" \
  --without-petsc \
  --without-hypre \
  --without-blaslapack \
  --without-cubes \
  --without-eleven \
  --without-petsc \
  --without-sundials \
  --without-x \
  --with-doxygen \
  --enable-debug \
  --disable-opt \
  --enable-implicit-template-instantiation \
  --disable-deprecated
  
make -j14
make -j14 install
```

#### Optimized build

```bash
cd $HOME/Code/IBSAMRAI2-2025.10.29
mkdir -p build-opt
cd build-opt

../configure \
  CFLAGS="-fPIC -O3 -march=native" \
  CXXFLAGS="-fPIC -O3 -march=native" \
  FFLAGS="-fPIC -O3 -march=native" \
  --prefix=$HOME/Applications/ibsamrai2-2025.10.29-opt \
  --with-CC="$HOME/Applications/petsc-3.23.3/bin/mpicc" \
  --with-CXX="$HOME/Applications/petsc-3.23.3/bin/mpicxx" \
  --with-F77="$HOME/Applications/petsc-3.23.3/bin/mpif90" \
  --with-hdf5="$HOME/Applications/petsc-3.23.3/" \
  --without-petsc \
  --without-hypre \
  --without-blaslapack \
  --without-cubes \
  --without-eleven \
  --without-sundials \
  --without-x \
  --enable-opt \
  --enable-implicit-template-instantiation \
  --disable-deprecated
  
make -j14
make -j14 install
```

Installed to `/home/sim/Applications/ibsamrai2-2025.10.29-debug` and `/home/sim/Applications/ibsamrai2-2025.10.29-opt`

Check `ls -l $HOME/Applications/ibsamrai2-2025.10.29-debug/lib/libSAMRAI*`

### Install **`IBAMR`** v0.18.0

```bash
# (ibamr) sim@DESKTOP-0:/mnt/c/Users/computer0$
cd $HOME/Code
wget https://github.com/IBAMR/IBAMR/archive/refs/tags/v0.18.0.tar.gz
tar xf v0.18.0.tar.gz
```

#### Debug build

```bash
cd $HOME/Code/IBAMR-0.18.0
mkdir -p build-debug
cd build-debug

cmake \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_C_COMPILER="$HOME/Applications/petsc-3.23.3/bin/mpicc" \
  -DCMAKE_CXX_COMPILER="$HOME/Applications/petsc-3.23.3/bin/mpicxx" \
  -DCMAKE_Fortran_COMPILER="$HOME/Applications/petsc-3.23.3/bin/mpif90" \
  -DPETSC_ROOT="$HOME/Applications/petsc-3.23.3/" \
  -DHDF5_ROOT="$HOME/Applications/petsc-3.23.3/" \
  -DHYPRE_ROOT="$HOME/Applications/petsc-3.23.3/" \
  -DSILO_ROOT="$HOME/Applications/petsc-3.23.3/" \
  -DSAMRAI_ROOT="$HOME/Applications/ibsamrai2-2025.10.29-debug" \
  -DLIBMESH_ROOT="$HOME/Applications/libmesh-1.7.8-debug" \
  -DLIBMESH_METHOD=DEVEL \
  -DIBAMR_FORCE_BUNDLED_BOOST=ON \
  -DIBAMR_FORCE_BUNDLED_EIGEN3=ON \
  -DIBAMR_FORCE_BUNDLED_MUPARSER=ON \
  -DCMAKE_INSTALL_PREFIX="$HOME/Applications/ibamr-0.18.0-debug" \
  ../
  
make -j14
make -j14 install
```

#### Optimized build

```bash
cd $HOME/Code/IBAMR-0.18.0
mkdir -p build-opt
cd build-opt

# Notes: explicitly add -march=native for performance
cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="$HOME/Applications/petsc-3.23.3/bin/mpicc" \
  -DCMAKE_CXX_COMPILER="$HOME/Applications/petsc-3.23.3/bin/mpicxx" \
  -DCMAKE_Fortran_COMPILER="$HOME/Applications/petsc-3.23.3/bin/mpif90" \
  -DPETSC_ROOT="$HOME/Applications/petsc-3.23.3/" \
  -DHDF5_ROOT="$HOME/Applications/petsc-3.23.3/" \
  -DHYPRE_ROOT="$HOME/Applications/petsc-3.23.3/" \
  -DSILO_ROOT="$HOME/Applications/petsc-3.23.3/" \
  -DSAMRAI_ROOT="$HOME/Applications/ibsamrai2-2025.10.29-opt" \
  -DLIBMESH_ROOT="$HOME/Applications/libmesh-1.7.8-opt" \
  -DLIBMESH_METHOD=OPT \
  -DIBAMR_FORCE_BUNDLED_BOOST=ON \
  -DIBAMR_FORCE_BUNDLED_EIGEN3=ON \
  -DIBAMR_FORCE_BUNDLED_MUPARSER=ON \
  -DCMAKE_INSTALL_PREFIX="$HOME/Applications/ibamr-0.18.0-opt" \
  -DCMAKE_CXX_FLAGS="-O3 -march=native" \
  -DCMAKE_Fortran_FLAGS="-O3 -march=native" \
  ../

make -j14
make -j14 install
```

check `ls -l $HOME/Applications/ibamr-0.18.0-debug/lib`  to see if have `libIBTK3d.a`, `libIBAMR3d.a`

## Git repo (Optional)

1. Fork the https://github.com/IBAMR/IBAMR repo.
2. Clone to local WSL -- using v0.18.0

```bash
wsl -d IBAMR_Sim

cd $HOME/Code
git clone https://github.com/<User>/IBAMR.git IBAMR_git
cd IBAMR_git

# Set remote upstream
git remote add upstream https://github.com/IBAMR/IBAMR.git
git remote -v
# Use v0.18.0
git fetch --all --tags
git reset --hard v0.18.0
```

3. Create branch `sim`

```bash
git checkout -b sim
git branch	# sim
git log -1 --oneline	# v0.18.0
```

4. Remake IBAMR optimized build

```bash
mkdir -p build-opt
cd build-opt

cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="$HOME/Applications/petsc-3.23.3/bin/mpicc" \
  -DCMAKE_CXX_COMPILER="$HOME/Applications/petsc-3.23.3/bin/mpicxx" \
  -DCMAKE_Fortran_COMPILER="$HOME/Applications/petsc-3.23.3/bin/mpif90" \
  -DPETSC_ROOT="$HOME/Applications/petsc-3.23.3/" \
  -DHDF5_ROOT="$HOME/Applications/petsc-3.23.3/" \
  -DHYPRE_ROOT="$HOME/Applications/petsc-3.23.3/" \
  -DSILO_ROOT="$HOME/Applications/petsc-3.23.3/" \
  -DSAMRAI_ROOT="$HOME/Applications/ibsamrai2-2025.10.29-opt" \
  -DLIBMESH_ROOT="$HOME/Applications/libmesh-1.7.8-opt" \
  -DLIBMESH_METHOD=OPT \
  -DIBAMR_FORCE_BUNDLED_BOOST=ON \
  -DIBAMR_FORCE_BUNDLED_EIGEN3=ON \
  -DIBAMR_FORCE_BUNDLED_MUPARSER=ON \
  -DCMAKE_INSTALL_PREFIX="$HOME/Applications/ibamr-sim-opt" \
  -DCMAKE_CXX_FLAGS="-O3 -march=native" \
  -DCMAKE_Fortran_FLAGS="-O3 -march=native" \
  ../
  
  make -j14 install
```

## Ready to go

Now we should have 2 main applications. Both can run official examples and user-defined tasks.

We use the offical read-only version as solver, with changed `main.cpp` files as input.

- **`~/Applications/ibamr-0.18.0-opt`: official stable version, read only.** 
- `~/Applications/ibamr-sim-opt`: customed version based on git repo at v0.18.0. 
  - `~/Code/IBAMR_git` set to `sim` branch to update your own code

------

### 
