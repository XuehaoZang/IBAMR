# Run IBAMR simulation

\\wsl.localhost\IBAMR_Sim\home\sim\Code\IBAMR_git\run\IBFEex9\build

Displacement (vector) expressions applied to the mesh
{X_0 - coord(Mesh)[0], X_1 - coord(Mesh)[1], X_2 - coord(Mesh)[2]}

``` bash
cd ~/Code/IBAMR_git/run/IBFEex9/build
rm -rf *
ln -s ../*.cpp .
ln -s ../input3d .
cmake \
  -DIBAMR_ROOT=$HOME/Applications/ibamr-0.18.0-opt \
  -DCMAKE_CXX_COMPILER=$HOME/Applications/petsc-3.23.3/bin/mpicxx \
  ../
  
make -j14

~/Applications/petsc-3.23.3/bin/mpirun -n 12 ./main3d ../input3d
```


``` bash
cd ~/Code/IBAMR_git/run/IBFEex9/build
rm -rf *
ln -s ../*.cpp .
ln -s ../input3d .
ln -s ../*.msh .
cmake \
  -DIBAMR_ROOT=$HOME/Applications/ibamr-0.18.0-opt \
  -DCMAKE_CXX_COMPILER=$HOME/Applications/petsc-3.23.3/bin/mpicxx \
  ../
  
make -j14

~/Applications/petsc-3.23.3/bin/mpirun -n 12 ./main3d ../input3d
```

Total number of systems: 4
--------------------------------------------------
System Index : 0
System Name  : [IB coordinates system]
System Type  : Explicit
Num Variables: 3
  -> Var 0 : X_0
  -> Var 1 : X_1
  -> Var 2 : X_2
--------------------------------------------------
System Index : 1
System Name  : [IB coordinate mapping system]
System Type  : Explicit
Num Variables: 3
  -> Var 0 : dX_0
  -> Var 1 : dX_1
  -> Var 2 : dX_2
--------------------------------------------------
System Index : 2
System Name  : [IB velocity system]
System Type  : Explicit
Num Variables: 3
  -> Var 0 : U_0
  -> Var 1 : U_1
  -> Var 2 : U_2
--------------------------------------------------
System Index : 3
System Name  : [IB force system]
System Type  : Explicit
Num Variables: 3
  -> Var 0 : F_0
  -> Var 1 : F_1
  -> Var 2 : F_2

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
use the U_tip mean instead of u tip max
Domain is defined as a cubic $[-L, L]^3$, $L= 4$

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


### 📐 更新后的无量纲与归一化几何体系

- $L_{scale} = 250\ \mu m$ 
-  **基准长度映射 (Scale Factor)**：$1.0 \text{ (sim)} = 250\ \mu m \text{ (real)}$

- **计算域总边长 ($L_{sim\_total}$)**：$8.0$ (区域是 $[-4, 4]^3$)。
  - 它对应真实的 $8 \times 250 = \mathbf{2000\ \mu m}$（正是你一开始核算的 $L_{real}$）。
  - $L_{sim\_total} = 8 \approx \mathbf{26.6} \bar{c}_{sim} \approx \mathbf{8} R_{sim}$ （边界缓冲依然非常充裕）。

- **特征频率 ($f_{sim}$)**：$1.0$ (对应真实的 $530\text{ Hz}$，周期 $T=1.0$)。

- **最大拍打角 ($\phi_{amp}$)**：$40^\circ = \mathbf{0.6981\text{ rad}}$。

- **翼尖旋转半径 ($R_{sim}$)**：$250 / 250 = \mathbf{1.0}$。

- **平均弦长 ($\bar{c}_{sim}$)**：【**关键变化**】$75 / 250 = \mathbf{0.3}$

- **密度 ($\rho_{sim}$)**：$1.0$

---

### 🧮 核心粘度计算 (基于 $U_{tip,mean}$ 和 $Re=1.3$)

**1. 无量纲平均翼尖速度 ($U_{tip,mean,sim}$)**：
速度公式不变，因为半径 $R_{sim}$ 和频率 $f_{sim}$ 都是 1.0：
$$U_{tip,mean,sim} = 4 \cdot f_{sim} \cdot \phi_{amp} \cdot R_{sim} = 4 \cdot 1.0 \cdot 0.6981 \cdot 1.0 = \mathbf{2.7924}$$

**2. 目标粘度 ($\mu_{sim}$)**：
现在我们将新的弦长 $\bar{c}_{sim} = 0.3$ 代入雷诺数公式反推粘度：
$$
\mu_{sim} = \frac{\rho_{sim} \cdot U_{tip,mean,sim} \cdot \bar{c}_{sim}}{Re} = \frac{1.0 \cdot 2.7924 \cdot 0.3}{1.3} = \frac{0.83772}{1.3} \approx \mathbf{0.6444}
$$


```text
// physical parameters
Re = 1.3
RHO = 1.0
MU = 0.6444       // 对应 R=1.0, c=0.3 的绝对无量纲参数

// grid spacing parameters
L = 4.0           // 计算域半宽 [-4, 4]
```


假设 $T$ 为一个完整的扑动周期（$T = 2\pi/\omega$）。
- $t=0$ (冲程起点)：$\phi=0^\circ, \dot{\phi}=0$；$\psi=45^\circ, \dot{\psi}=0$。此时 $y>0$ 的部分是前缘 (leading edge)，展向与 X 轴重合 ($x>0$)。翅膀准备开始前挥 (Forward stroke/Downstroke)。
- $t=1/4 T$ (第一次冲程反转)：当 stroke 从 $0^\circ$ 增大到 $40^\circ$ 的过程中，pitch 从 $45^\circ$ 逐渐增大到 $90^\circ$。达到临界点时（$\phi=40^\circ, \psi=90^\circ$），展向与 X 轴夹角为 $40^\circ$。此时翅膀完全垂直于运动方向，弦向 (chord) 严格指向 $-Z$ 向。
- $t=1/2 T$ (回桨中点)：进入回挥 (Backstroke)，stroke 从 $40^\circ$ 减小回到初始位置 $\phi=0^\circ$。在此过程中，pitch 继续增大，达到最大反向迎角 $\psi=135^\circ$。此时展向重新与 X 轴重合。前缘依然引导气流，但翅膀此时正在向相反方向拍打。
- $t=3/4 T$ (第二次冲程反转)：stroke 从 $0^\circ$ 继续减小到 $-40^\circ$。pitch 从 $135^\circ$ 逐渐回落到 $90^\circ$。达到临界点时（$\phi=-40^\circ, \psi=90^\circ$），展向与 X 轴夹角为 $-40^\circ$。翅膀再次垂直于运动方向切风，此时弦向 (chord) 严格指向 $+Z$ 向，准备进行下一次翻转。
- $t=T$ (周期结束)：stroke 从 $-40^\circ$ 回到 $0^\circ$，pitch 从 $90^\circ$ 减小回到 $45^\circ$。翅膀完全回到 $t=0$ 时的初始位姿和速度，完成一个完美的“8”字型运动周期。

**参数定义：**
$$S(t) = 1 - e^{-t / t_{ramp}}$$
$$\dot{S}(t) = \frac{1}{t_{ramp}} e^{-t / t_{ramp}}$$
$$C = \frac{1}{\pi \tau}$$
* $\phi_{amp} = 40^\circ$ （拍打振幅）
* $\psi_0 = 90^\circ$ （俯仰中心角/垂直切风角）
* $\psi_{amp} = 45^\circ$ （俯仰振幅）

#### 拍动角 (Stroke, 绕 Z 轴)
$$\phi(t) = \phi_{amp} \cdot S(t) \cdot \sin(\omega t)$$
$$\dot{\phi}(t) = \phi_{amp} \left[ \dot{S}(t) \sin(\omega t) + S(t) \omega \cos(\omega t) \right]$$

#### 俯仰角 (Pitch, 绕展向 X 轴)
$t=0, \cos(0)=1$ 时，$\psi = 90^\circ - 45^\circ = 45^\circ$；当 $t=T/2, \cos(\pi)=-1$ 时，$\psi = 90^\circ - (-45^\circ) = 135^\circ$。
$$\psi(t) = \psi_0 - \frac{\psi_{amp}}{\tanh(C)} \tanh(C \cos(\omega t))$$

$$\dot{\psi}(t) = \frac{\psi_{amp}}{\tanh(C)} \left[ 1 - \tanh^2(C \cos(\omega t)) \right] \left( C \omega \sin(\omega t) \right)$$


