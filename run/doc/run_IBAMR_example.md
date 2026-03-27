# Run IBAMR simulation

We will import an official example project and change to the flapping wing setting we need to simulate.

## Test with example

#### Find example files

The example source code can be found in the file named `example.cpp` in each of these examples:

```bash
IBAMR_SRC_DIR
|
├── examples
│   ├── ...
│   ├── IB
│   │   ├── explicit
│   │   └── implicit
│   ├── IBFE
│   │   └── explicit
```

``` bash
# find source dir 
find ~ -name "example.cpp"
```

`IBAMR_SRC_DIR` = `$Home/Code/IBAMR-0.18.0/`

#### Copy example to dir

We will copy over **example 2** (a 3D thin shell in fluid) of IBFE, Immersed Boundary Finite Element.

(ex0 and ex1 are 2d without an `input3d` file.) 

``` bash
# PS C:\Users\computer0> 
wsl -d IBAMR_Sim
# (base) sim@DESKTOP-0:/mnt/c/Users/computer0$ 
conda activate ibamr
# (ibamr) sim@DESKTOP-0:/mnt/c/Users/computer0$
cd $Home		# $Home = /home/sim/

mkdir ~/ex2_test
cd ~/ex2_test
cp /home/sim/Code/IBAMR-0.18.0/examples/IBFE/explicit/ex2/* .
ls
```

#### Configure Cmake

``` bash
cat << 'EOF' > CMakeLists.txt
cmake_minimum_required(VERSION 3.15.0)
project(ex2_reference)

find_package(IBAMR 0.18.0 REQUIRED)

set(SOURCE_FILES example.cpp)

add_executable(main3d ${SOURCE_FILES})

target_link_libraries(main3d IBAMR::IBAMR3d)
set(CMAKE_CXX_FLAGS "${IBAMR_CXX_FLAGS} -O3 -march=native")
EOF
```

Use ` cat CMakeLists.txt` to check if it exists.

#### Cmake

```bash
mkdir build
cd build

cmake \
  -DIBAMR_ROOT=$HOME/Applications/ibamr-0.18.0-opt \
  -DCMAKE_CXX_COMPILER=$HOME/Applications/petsc-3.23.3/bin/mpicxx \
  ../
  
ls -F
```

Use `ls -F` to check if  `main3d*` appears as green executable file.

#### Test run

We will use **4 cores** to test run **~50 steps** to visualize.

``` bash
~/Applications/petsc-3.23.3/bin/mpirun -n 4 ./main3d ../input3d
```

---

## Modify with IBFE ex2

### Copy example to dir

#### Set up dir

``` bash
# PS C:\Users\computer0> 
wsl -d IBAMR_Sim
# (base) sim@DESKTOP-0:/mnt/c/Users/computer0$ 
conda activate ibamr
# (ibamr) sim@DESKTOP-0:/mnt/c/Users/computer0$
cd $Home		# $Home = /home/sim/
mkdir ~/robowing
cd ~/robowing
mkdir run1
```

#### Copy example files

``` bash
cd ~/robowing/run1
cp /home/sim/Code/IBAMR-0.18.0/examples/IBFE/explicit/ex2/example.cpp .
cp /home/sim/Code/IBAMR-0.18.0/examples/IBFE/explicit/ex2/input3d .
```

#### Configure `CMake`

Notes: 

- find_package: look for ibamr 0.18.0 library
- set source file `example.cpp`
- set 3d executable `main3d`
- link to `IBAMR3d` library
- add flag `-O3 -march=native` for performance
- Earlier we used `cmake` to set up the `ibamr-0.18.0-opt`. reference to: https://ibamr.github.io/linking

```bash
cat << 'EOF' > CMakeLists.txt
cmake_minimum_required(VERSION 3.15.0)
project(robowing)

find_package(IBAMR 0.18.0 REQUIRED)

set(SOURCE_FILES example.cpp)

add_executable(main3d ${SOURCE_FILES})

target_link_libraries(main3d IBAMR::IBAMR3d)

set(CMAKE_CXX_FLAGS "${IBAMR_CXX_FLAGS} -O3 -march=native")
EOF
```

Use ` cat CMakeLists.txt` to check if it exists.

```bash
# cmake
mkdir build
cd build

cmake \
  -DIBAMR_ROOT=$HOME/Applications/ibamr-0.18.0-opt \
  -DCMAKE_CXX_COMPILER=$HOME/Applications/petsc-3.23.3/bin/mpicxx \
  ../
  
make -j14
```

Use `ls -F` to check if  `main3d*` appears as green executable file.

#### Test run

``` bash
./main3d ../input3d
```

should see the simulation running:

```bash
IBFEMethod: mesh part 0 is using FIRST order LAGRANGE finite elements.
...
+++++++++++++++++++++++++++++++++++++++++++++++++++
At beginning of timestep # 1
Simulation time is 0.00390625
```

### Define computational domain

Domain is defined as a cubic $[0, L]^3$, $L= 1$

- Geometry in microns

  - wing span $b= 250$, hinge $r= 70$

  - wingtip radius $R=320$

  - max chord $c= 150$, mean chord $\bar c= 75$

  - $L_{real}= 2000\text{um} \approx 26 \bar c \approx 6R$

  - Freqency $f=530$Hz

  - $ U_{tip,max} = 2\pi f \phi_{amp}R = \approx \mathbf{0.74 \text{ m/s}}$

    - $$\omega(t) = \frac{d\phi}{dt} = -\phi_{amp} \cdot (2\pi f) \cdot \sin(2\pi f t)$$
    - $$\omega_{max} = 2\pi f \cdot \phi_{amp}$$

  - $$
    Re = \frac{\rho U_{tip,max} \bar{c}}{\mu} = \frac{ 2\pi f \phi_{amp} R \bar{c}}{\nu}=3.67
    $$

- Non-dimensional and normalized geometry

  - **计算域边长 ($L_{sim}$)**：$1.0$ (对应 $2000\ \mu m$)。

    - $L= 1 \approx 26 \bar c \approx 6R$

  - **特征频率 ($f_{sim}$)**：$1.0$ (对应 $530\text{ Hz}$，即一个周期 $T=1.0$)。

  - **翼尖旋转半径 ($R_{sim}$)**：$320 / 2000 = \mathbf{0.16}$。

  - **平均弦长 ($\bar{c}_{sim}$)**：$75 / 2000 = \mathbf{0.0375}$。

  - **密度 ($\rho$)**：$1.0$

  - Assume $Re=4$, 粘度
    $$
    \mu = \frac{\omega_{max} \cdot R_{tip} \cdot \bar{c}}{Re} =\frac{2\pi f \phi_{amp}R\cdot \bar{c}}{Re}=0.00658
    $$


### Define static wing mesh in `.cpp`

Define non-dimensional and normalized wing mesh

```c++
/* Define normalized wing geometry mesh */
void def_wing_mesh(Mesh& mesh) {
    const double s = 1.0 / 2000.0; // Scaling factor
    
    // Axis center at (0.5, 0.5, 0.5)
    Point p1(0.5, 0.5 + 70.0*s,  0.5); 
    Point p7(0.5, 0.5 + 320.0*s, 0.5);
    Point p4(0.5, 0.5 + 128.0*s, 0.5 - 150.0*s); // Initially pointing -Z

    mesh.add_point(p1);
    mesh.add_point(p4);
    mesh.add_point(p7);

    Elem* elem = mesh.add_elem(new Tri3);
    for (int i = 0; i < 3; ++i) elem->set_node(i) = mesh.node_ptr(i);
    mesh.prepare_for_use();
}
```

in `main`, change structure mesh setting

```c++
int
main(int argc, char* argv[])
{
...
#if (NDIM == 3)
        Mesh mesh(init_info.comm(), 3);
        def_wing_mesh(mesh); /* Call custom wing mesh definition */
#endif
```

### Prescribe wing kinematics (Lagrangian)  in `.cpp`

```c++
/* Prescribed stroke-pitch kinematics: tanh-based pitch and cosine stroke */
void update_lagrangian_structure(Mesh& mesh, double time, double dt, void* ctx){
    const double s = 1.0 / 2000.0;     // Scaling factor (Microns to L=1.0)
    const double f = 530.0;            // Frequency (Hz)
    const double tau = 0.1;            // Pitch transition parameter
    const double C = 1.0 / (M_PI * tau);
    const double xc = 0.5, yc = 0.5, zc = 0.5; // Rotation center

    // Angle calculations (Radians)
    double phi = (40.0 * M_PI / 180.0) * cos(2.0 * M_PI * f * time);
    double inner = C * sin(2.0 * M_PI * f * time);
    double psi = - (45.0 * M_PI / 180.0) / tanh(C) * tanh(inner);

    for (auto& node : mesh.node_ptr_range()) {
        const int id = node->id();
        double r = 0.0, h = 0.0; 

        // Map relative positions from initial mesh
        if (id == 0) { r = 70.0*s;  h = 0.0; }       // P1 (Root)
        if (id == 1) { r = 128.0*s; h = -150.0*s; } // P4 (Chord tip)
        if (id == 2) { r = 320.0*s; h = 0.0; }       // P7 (Wing tip)

        // 1. Pitching (rotation around span axis)
        double px = -h * sin(psi);
        double pz =  h * cos(psi);

        // 2. Stroke rotation (XY plane sweep around center)
        (*node)(0) = xc + px * cos(phi) - r * sin(phi);
        (*node)(1) = yc + px * sin(phi) + r * cos(phi);
        (*node)(2) = zc + pz;
    }
}
```

### Set parameters in `input3d` (SAMRAI script)

#### Physical parameters

using $\mu$ calculated from $Re=4$

```c++
// physical parameters
MU  = 0.00658 /* Re = 3 */ /* Adjust: mu = (rho * U_tip * c_bar) / Re  */
RHO = 1.0
L   = 1.0
```

#### Compulational domain

```c++
CartesianGeometry {
   domain_boxes = [ (0,0,0),(N - 1,N - 1,N - 1) ]
   x_lo = 0,0,0
   x_up = L,L,L
   periodic_dimension = 0,0,0
}
```

- cubic $[0,L]^3$，网格32x32x32

#### Grid spacing

```c++
// grid spacing parameters
MAX_LEVELS = 3                   // maximum number of levels in locally refined grid
REF_RATIO  = 4                   // refinement ratio between levels
N = 32                           // actual    number of grid cells on coarsest grid level
NFINEST = (REF_RATIO^(MAX_LEVELS - 1))*N       // effective number of grid cells on finest grid level
DX0 = L/N                                      // mesh width on coarsest grid level
DX  = L/NFINEST                                // mesh width on finest   grid level
MFAC = 2.0                      // ratio of Lagrangian mesh width to Cartesian mesh width
```

- 弦长分辨率 (Points per Chord)：特征弦长 $\bar{c}$ 上布置 **15 - 20 个网格点** 。

  - 无量纲化的平均弦长 $\bar{c} \approx 0.0375$ 。

  - **计算目标**：为了达到 20 个点的分辨率，最细网格步长 $dx$ 应为：

    $$dx = \frac{0.0375}{20} \approx 0.001875$$

  - **换算为全局网格数**：在 $L=1.0$ 的计算域中，对应的等效最细网格数 $N_{finest}$ 为：

    $$N_{finest} = \frac{1.0}{0.001875} \approx 533$$

- 设定 $N_{finest} = 32 \times 4^{(3-1)} = 32 \times 16 = 512$

  - `N = 32` (Coarsest grid level)：最粗一层的网格数。
  - `MAX_LEVELS = 3`：网格细化的层级数。
  - `REF_RATIO = 4`：每一级的细化比例。
  - `MFAC = 2` (Mesh Factor) ：翅膀三角形（拉格朗日网格）与背景流体网格（欧拉网格）的疏密比例。

#### Gridding algorithm

```c++
GriddingAlgorithm {
   max_levels = MAX_LEVELS
   ratio_to_coarser {
      level_1 = REF_RATIO,REF_RATIO,REF_RATIO
      level_2 = REF_RATIO,REF_RATIO,REF_RATIO
   }
   largest_patch_size {
      level_0 = 128, 128, 128  // all finer levels will use same values as level_0
   }

   smallest_patch_size {
      level_0 =   8,  8,  8  // all finer levels will use same values as level_0
   }

   efficiency_tolerance = 0.85e0  // min % of tag cells in new patch level
   combine_efficiency   = 0.85e0  // chop box if sum of volumes of smaller boxes < efficiency * vol of large box
}
```

- 这里`ratio_to_coarser`的级别数目要和之前的`MAX_LEVELS`对应

#### Solver and timestep

```c++
// solver parameters
CFL_MAX                    = 0.1               // maximum CFL number
DT                         = 0.1*DX           // maximum timestep size
START_TIME                 = 0.0e0             // initial simulation time
END_TIME                   = 100*DT            // final simulation time
GROW_DT                    = 1.0e0             // growth factor for timesteps
NUM_CYCLES                 = 2              // number of cycles of fixed-point iteration

NORMALIZE_PRESSURE         = TRUW             // whether to explicitly force the pressure to have mean zero
ERROR_ON_DT_CHANGE         = TRUE              // whether to emit an error message if the time step size changes
VORTICITY_TAGGING          = TRUE             // whether to tag cells for refinement based on vorticity thresholds
TAG_BUFFER                 = 2                 // size of tag buffer used by grid generation algorithm
REGRID_CFL_INTERVAL        = 0.5               // regrid whenever any material point could have moved 0.5 meshwidths since previous regrid
```

- 低雷诺数下，将 `CFL_MAX` 降至 0.1，并减小时间步长 `DT = 0.1 * DX`， `GROW_DT = 1.0`
- 无量纲时间是0-1，表示一个周期
- 开启 `VORTICITY_TAGGING = TRUE`。当翅膀运动产生涡量时，IBAMR 会自动在涡量集中的区域布置最细层网格

#### IBFE methods

```c++
IBFEMethod {
   ...
   /*Disable check to allow initialization on Level 0 */
   node_outside_patch_check = "NONE"
}
IBHierarchyIntegrator {
   ...
   /* Ensure the grid regrids at the very first timestep */
   regrid_at_initial_time = TRUE 
}
```

#### Boundary conditions

```c++
/* Velocity component: Set all faces to zero-gradient (Extrapolated) */
VelocityBcCoefs_<x> {
   acoef_function_<x> = "0.0"
   bcoef_function_<x> = "1.0"
   gcoef_function_<x> = "0.0"
}
```

- **acoef=0.0 / bcoef=1.0**：定义了 **Neumann** 边界条件。在 IBAMR 中，这通常对应于 `EXTRAPOLATED` 外推条件，表示边界上的物理量梯度为 0 。
- **gcoef=0.0**：这表示边界上没有外部驱动力 。

### Run simulation

Using 4 cores

``` bash
cd ~/robowing/run1/build
./main3d ../input3d
```

should see the simulation running:

```bash
IBFEMethod: mesh part 0 is using FIRST order LAGRANGE finite elements.
...
+++++++++++++++++++++++++++++++++++++++++++++++++++
At beginning of timestep # 1
Simulation time is 0.00390625
```

## Modify with eel3D

ConstraintIB/eel3D (a 3D eel swimming in fluid)

----

- `robowing.cpp`
- `wingKinematics.h`
- `wingKinematics.cpp`
- `wing.vertex`
- `input3d`
- `CMakeLists.txt`

### Define computational domain

Domain is defined as a cubic $[0, L]^3$, $L= 1$

- Geometry in microns

  - wing span $b= 250$, hinge $r= 70$

  - wingtip radius $R=320$

  - max chord $c= 150$, mean chord $\bar c= 75$

  - $L_{real}= 2000\text{um} \approx 26 \bar c \approx 6R$

  - Freqency $f=530$Hz

  - $ U_{tip,max} = 2\pi f \phi_{amp}R = \approx \mathbf{0.74 \text{ m/s}}$

    - $$\omega(t) = \frac{d\phi}{dt} = -\phi_{amp} \cdot (2\pi f) \cdot \sin(2\pi f t)$$
    - $$\omega_{max} = 2\pi f \cdot \phi_{amp}$$

  - $$
    Re = \frac{\rho U_{tip,max} \bar{c}}{\mu} = \frac{ 2\pi f \phi_{amp} R \bar{c}}{\nu}=3.67
    $$

- Non-dimensional and normalized geometry

  - **计算域边长 ($L_{sim}$)**：$1.0$ (对应 $2000\ \mu m$)。

    - $L= 1 \approx 26 \bar c \approx 6R$

  - **特征频率 ($f_{sim}$)**：$1.0$ (对应 $530\text{ Hz}$，即一个周期 $T=1.0$)。

  - **翼尖旋转半径 ($R_{sim}$)**：$320 / 2000 = \mathbf{0.16}$。

  - **平均弦长 ($\bar{c}_{sim}$)**：$75 / 2000 = \mathbf{0.0375}$。

  - **密度 ($\rho$)**：$1.0$

  - Assume $Re=4$, 粘度
    $$
    \mu = \frac{\omega_{max} \cdot R_{tip} \cdot \bar{c}}{Re} =\frac{2\pi f \phi_{amp}R\cdot \bar{c}}{Re}=0.00658
    $$

### Set parameters in `input3d` (SAMRAI script)

#### Physical parameters

using $\mu$ calculated from $Re=4$

```c++
// physical parameters
MU  = 0.00658 /* Re = 4 */ /* Adjust: mu = (rho * U_tip * c_bar) / Re  */
RHO = 1.0
L   = 1.0
```

#### Compulational domain

```c++
CartesianGeometry {
   domain_boxes = [ (0,0,0),(N - 1,N - 1,N - 1) ]
   x_lo = 0,0,0
   x_up = L,L,L
   periodic_dimension = 0,0,0
}
```

- cubic $[0,L]^3$

#### Boundary conditions

```c++
/* Velocity component: Set all faces to zero-gradient (Extrapolated) */
VelocityBcCoefs_<x> {
   acoef_function_<x> = "0.0"
   bcoef_function_<x> = "1.0"
   gcoef_function_<x> = "0.0"
}
```

- **acoef=0.0 / bcoef=1.0**：定义了 **Neumann** 边界条件。在 IBAMR 中，这通常对应于 `EXTRAPOLATED` 外推条件，表示边界上的物理量梯度为 0 。
- **gcoef=0.0**：这表示边界上没有外部驱动力 。

#### Grid spacing

```c++
// grid spacing parameters
MAX_LEVELS = 3                   // maximum number of levels in locally refined grid
REF_RATIO  = 4                   // refinement ratio between levels
N = 32                           // actual    number of grid cells on coarsest grid level
NFINEST = 512       // effective number of grid cells on finest grid level
DX0 = L/N                                      // mesh width on coarsest grid level
DX  = L/NFINEST                                // mesh width on finest   grid level
MFAC = 1.0                      // ratio of Lagrangian mesh width to Cartesian mesh width
```

- 弦长分辨率 (Points per Chord)：特征弦长 $\bar{c}$ 上布置 **15 - 20 个网格点** 。

  - 无量纲化的平均弦长 $\bar{c} \approx 0.0375$ 。

  - **计算目标**：为了达到 20 个点的分辨率，最细网格步长 $dx$ 应为：

    $$dx = \frac{0.0375}{20} \approx 0.001875$$

  - **换算为全局网格数**：在 $L=1.0$ 的计算域中，对应的等效最细网格数 $N_{finest}$ 为：

    $$N_{finest} = \frac{1.0}{0.001875} \approx 533$$

- 设定 $N_{finest} = 32 \times 4^{(3-1)} = 32 \times 16 = 512$

  - `N = 32` (Coarsest grid level)：最粗一层的网格数。
  - `MAX_LEVELS = 3`：网格细化的层级数。
  - `REF_RATIO = 4`：每一级的细化比例。
  - `MFAC = 2` (Mesh Factor) ：翅膀三角形（拉格朗日网格）与背景流体网格（欧拉网格）的疏密比例。

#### Gridding algorithm

```c++
GriddingAlgorithm {
   max_levels = MAX_LEVELS           // Maximum number of levels in hierarchy.
   ratio_to_coarser {
      level_1 = REF_RATIO, REF_RATIO, REF_RATIO      // vector ratio to next coarser level
      level_2 = REF_RATIO, REF_RATIO, REF_RATIO
   }

   largest_patch_size {
      level_0 = 512, 512, 512 // largest patch allowed in hierarchy
                       // all finer levels will use same values as level_0...
   }

   smallest_patch_size {
      level_0 =  8, 8, 8  // smallest patch allowed in hierarchy
                           // all finer levels will use same values as level_0...
   }

   allow_patches_smaller_than_minimum_size_to_prevent_overlaps = TRUE
   efficiency_tolerance   = 0.6e0    // min % of tag cells in new patch level
   combine_efficiency     = 0.8e0    // chop box if sum of volumes of smaller
                                      // boxes < efficiency * vol of large box
}
```

- 这里`ratio_to_coarser`的级别数目要和之前的`MAX_LEVELS`对应

#### Solver and timestep

```c++
// solver parameters
DELTA_FUNCTION       = "IB_4"
START_TIME           = 0.0e0                      // initial simulation time
END_TIME             = 1.0e0                       // final simulation time
MAX_INTEGRATOR_STEPS = 100000                      // Max no of steps.
GROW_DT              = 2.0e0                      // growth factor for timesteps
NUM_CYCLES           = 1                          // number of cycles of fixed-point iteration. For ConstraintIBMethod set it to 1.
CONVECTIVE_OP_TYPE   = "PPM"                      // convective differencing discretization type
CONVECTIVE_FORM      = "ADVECTIVE"                // how to compute the convective terms
NORMALIZE_PRESSURE   = FALSE                      // whether to explicitly force the pressure to have mean zero
CFL_MAX              = 0.1                        // maximum CFL number
DT_MAX               = 0.001                      // maximum timestep size
VORTICITY_TAGGING    = FALSE                       // whether to tag cells for refinement based on vorticity thresholds
TAG_BUFFER           = 2                          // size of tag buffer used by grid generation algorithm
REGRID_CFL_INTERVAL  = 0.3                        // regrid whenever any material point could have moved 0.5 meshwidths since previous regrid
```

- 低雷诺数下，将 `CFL_MAX` 降至 0.1
- 无量纲时间是0-1，表示一个周期
- 先关闭 `VORTICITY_TAGGING = FALSE`：当翅膀运动产生涡量时，IBAMR 会自动在涡量集中的区域布置最细层网格

#### Constraint IB Kinematics

```c++
ConstraintIBKinematics {
     wing {
     structure_names                  = "wing"
     structure_levels                 =  MAX_LEVELS - 1
     calculate_translational_momentum = 1,1,1
     calculate_rotational_momentum    = 0,0,1
     lag_position_update_method       = "CONSTRAINT_POSITION"
     tagged_pt_identifier             = MAX_LEVELS - 1, 0  // level, relative idx of lag point
    }
}
```

### Set wing kinematics in `wingKinematics.cpp/.h`

$$\phi(t) = \phi_{amp} \cdot S(t) \cdot \sin(\omega t)$$

$$\dot{\phi}(t) = \phi_{amp} \left[ \dot{S}(t) \sin(\omega t) + S(t) \omega \cos(\omega t) \right]$$

- 软启动$S(t) = 1 - e^{-t / t_{ramp}},\dot{S}(t) = \frac{1}{t_{ramp}} e^{-t / t_{ramp}}$

$$\psi(t) = \frac{\psi_{amp}}{\tanh(C)} \tanh(C \cos(\omega t))$$

$$\dot{\psi}(t) = \frac{\psi_{amp}}{\tanh(C)} \left[ 1 - \tanh^2(C \cos(\omega t)) \right] \left( -C \omega \sin(\omega t) \right)$$

- $C = \frac{1}{\pi \tau}$，$\tau$代表翅膀翻转时间比例

**软启动：**

$t=0$时，$\phi=\dot \phi=0$，$\psi=45,\dot\psi=0$。

### Set wing initial mesh in `wing.vertex`

1. 展向沿 X 轴，弦向沿 Z 轴；Stroke 绕竖直 Z 轴旋转，Pitch 扭转绕展向 X 轴旋转 。
2. 坐标系中心为$[0.5,0.5,0.5]，$初始格点`wing_flat`生成在$X-Z$平面($Y=0.5$)，长边(LE,span)在$X=Z=0.5$轴上，翅膀挥动时一直$X>0.5$。
3. 在浸入边界法中，对于果蝇这种极薄的柔性薄膜，**“零厚度的二维流形（2D Manifold in 3D space）”是最优、最稳定且最符合标准的做法**。如果添加微小厚度，容易发散。

```bash
# first generate wing_flat.vertex
# compile .cpp to executable
g++ generateWingVertex.cpp -o gen_wing
# run executable
./gen_wing

# rotate to initial wing.vertex
g++ rotateWing.cpp -o rotate_wing
./rotate_wing
```

Vertex info are stored in `./wing.vertex`

### Compile and run simulation

``` bash
cd ~/robowing/run2/build
rm -rf *
ln -s ../wing.vertex .
ln -s ../input3d .
ln -s ../wingKinematics.h .

cmake \
  -DIBAMR_ROOT=$HOME/Applications/ibamr-0.18.0-opt \
  -DCMAKE_CXX_COMPILER=$HOME/Applications/petsc-3.23.3/bin/mpicxx \
  ../
  
make -j14

~/Applications/petsc-3.23.3/bin/mpirun -n 4 ./main3d ../input3d
```

- 如果修改了.cpp，需要重新清空build文件夹重新cmake。

- 如果只修改了input3d文件，不需要make，可以先`rm -rf viz_wing output`清空历史数据后，直接运行。

### Visualize results

Move the files to server to access in Windows

```bash
sudo mkdir -p /mnt/x
sudo mount -t drvfs X: /mnt/x
cd /mnt/x/antenna/sim
mkdir -p run2
```

Reset path in samrai files

```bash
# local
cd ~/robowing/run2
rm -r /mnt/x/antenna/sim/run2
cp -r ~/robowing/run2/* /mnt/x/antenna/sim/run2
cd /mnt/x/antenna/sim/run2
```

####  `ParaView`

-  install via `ParaView-6.1.0-RC1-Windows-Python3.12-msvc2017-AMD64.msi` from official website

-  **File > Open >** `X:\antenna\sim\ex2_viz\viz_IB3d\visit_dump.00064\summary.samrai`
   - Using `VisIt SAMRAI Reader`
-  **Properties > Check** U,P,...
-  **Slice** to visualize field
   - Solid Color > U
-  **File > Open > **`output.ex2`
   - Using `ExodusII Reader`

#### `VisIt` V3.3.3

1. Add plot
   - **Add** > Pseudocolor >  **`U`**（速度大小）、**`P`**（压力）或 **`Omega`**（涡量）
     - Plot 列表中会出现一行绿色的字，代表准备就绪

2. Add Operator
   - 选中Pseudocolor **`U`** > **Operators** > Slicing > Slice
   - 展开 `Pseudocolor` 左边的小加号 `+`，双击出现的 **`Slice`**。
     - 在弹出的属性窗口中，将 **Normal（法向量）** 设为 `0, 1, 0` (Y 轴切面)。
     - 将 **Origin（原点）** 选择为 `Point`，并输入 `0.5 0.5 0.5`。
     - 点击 **Apply** 并关闭窗口。
3. **Draw**
4. Play video
