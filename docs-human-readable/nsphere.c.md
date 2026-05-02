# `src/c/nsphere.c` 中文阅读文档

本文档面向维护者，目标不是复述每一行代码，而是把 `nsphere.c` 这个 1.9 万行 C 文件拆成可理解的模块，说明它在做什么、数据如何流动、`main()` 为什么会膨胀到一万多行，以及后续读代码或重构时应该抓住哪些主线。

本文基于当前工作区中的 `src/c/nsphere.c`。该文件总长约 19252 行，其中 `main()` 从第 6578 行开始，到第 16850 行左右结束。

## 一句话总结

`nsphere.c` 是 NSphere 的核心模拟器：它读取命令行参数，生成或读取球对称暗物质 halo 的初始条件，然后按半径排序粒子，用球壳质量近似计算自引力，选择一种时间积分器推进粒子，并可选地加入 SIDM 自相互作用散射、各向异性初始分布、重启、扩展、全粒子快照、轨迹跟踪和后处理输出。

真正的核心思想很简单：

1. 每个粒子只保存半径 `r`、径向速度 `v_r`、角动量 `L` 等低维状态。
2. 每个时间步按半径排序粒子。
3. 排序后的下标近似代表 enclosed mass，即 `M(r)`。
4. 引力只依赖球内质量和当前半径：`a_grav ~ -G M(r) / r^2`。
5. 再加上角动量带来的离心项 `L^2 / r^3`。
6. 用选定积分器推进 `r` 和 `v_r`。
7. 按设定间隔把快照、轨迹、能量、诊断文件写到 `data/`。

复杂度来自外围功能，而不是这个物理主循环本身：NFW/cored/Hernquist 三类初始分布、Eddington 反演、Osipkov-Merritt 各向异性、SIDM Monte Carlo 散射、重启/扩展、双精度滚动缓冲、快照后处理、排序算法基准测试，全都塞进了同一个文件和同一个 `main()`。

## 文件结构地图

大致可以按以下区域阅读：

| 行号范围 | 主要内容 |
| --- | --- |
| 1-185 | 头文件、OpenMP stub、基础结构体、函数前置声明 |
| 191-817 | 全局配置、profile 参数、SIDM 参数、随机种子、粒子块缓冲、轨迹缓冲 |
| 878-1537 | 全局数据清理、双精度快照缓冲、轨迹缓冲、重启文件一致性检查 |
| 1537-2110 | 文件截断、引力/离心力、能量诊断相关工具 |
| 2160-2369 | 命令行参数校验、`--help` 打印、错误退出 |
| 2410-3279 | `all_particle_data` 等二进制快照写入/读取、restart 粒子恢复、chosen particle 文件 |
| 3336-4139 | 时间积分器辅助函数：自适应 leapfrog、Levi-Civita 正则化等 |
| 4139-4558 | 初始条件读写、排序描述、二进制 printf/scanf 工具 |
| 4559-5368 | `--save` 解析、ID 重排、高斯卷积、restart snapshot 检查 |
| 5411-5864 | SIDM step 调度、磁盘空间检查、负 distribution function 修正、文件名解析 |
| 5974-6564 | Hernquist 初始条件生成、profile potential wrapper |
| 6578-16850 | `main()`：命令行、IC、模拟循环、输出、后处理、清理 |
| 16863-17442 | cored/NFW/Hernquist 密度、势能、Eddington integrand、OM integrand |
| 17480-18478 | 粒子排序算法：insertion、quadsort、parallel radix、转置缓冲 |
| 18518-19176 | 3D 向量、总能量、SIDM 串行与并行 graph-coloring 实现 |
| 19193-19252 | Hernquist 解析密度、势能、各向异性分布函数 |

## 粒子数据布局

核心粒子数组是 `double **particles`，分配为 7 行、`npts` 列：

| 行 | 含义 | 备注 |
| --- | --- | --- |
| `particles[0][i]` | 半径 `r` | 单位通常是 kpc |
| `particles[1][i]` | 初始生成阶段是速度大小 `|v|`，进入模拟后是径向速度 `v_r` | 这是最容易踩坑的字段，语义会变 |
| `particles[2][i]` | 带符号角动量 `L` | 生成后转换到 `kpc^2/Myr` |
| `particles[3][i]` | 粒子 ID / final rank ID | 排序后用它重建 `inverse_map` |
| `particles[4][i]` | `mu = v_r / |v|` | 主要用于初始条件方向，restart 后不再重要 |
| `particles[5][i]` | `cos(phi)` | SIDM 需要切向速度方向 |
| `particles[6][i]` | `sin(phi)` | 与 `cos(phi)` 一起持久保存 |

`particles[1]` 的语义变化是理解本文件的第一大关键：初始条件采样时先采样速度大小，后面会乘以 `mu` 并转换单位，变成径向速度。restart 路径直接从文件恢复 `Vrad`，所以不会再经历这个转换。

此外还有几个重要块数组：

| 全局数组 | 作用 |
| --- | --- |
| `L_block` | 一个 snapshot buffer block 中每个粒子的角动量 |
| `Rank_block` | 每个粒子的当前半径排序 rank |
| `R_block` | 半径 |
| `Vrad_block` | 径向速度 |
| `phi_block` | 方位角 `phi`，由 `atan2(sin_phi, cos_phi)` 得到 |
| `scatter_count_block` | SIDM 每粒子散射计数 |
| `ID_block` | 粒子 ID |
| `inverse_map` | `inverse_map[particle_id] = 当前排序数组下标` |

`inverse_map` 在每次排序后都必须重建，否则轨迹跟踪、低角动量粒子跟踪、debug 粒子能量都会找错粒子。

## 全局配置的主要类别

这个文件大量依赖全局变量。重要类别如下：

| 类别 | 代表变量 |
| --- | --- |
| 输出控制 | `g_doDebug`、`g_doDynPsi`、`g_doDynRank`、`g_doAllParticleData` |
| restart/extend | `g_doRestart`、`g_doSimRestart`、`g_doSimExtend`、`skip_file_writes`、`g_restart_mode_active` |
| SIDM | `g_enable_sidm_scattering`、`g_sidm_execution_mode`、`g_sidm_kappa`、`g_sidm_seed`、`g_total_sidm_scatters` |
| 随机种子 | `g_master_seed`、`g_initial_cond_seed`、`g_sidm_seed`、`g_attempt_load_seeds` |
| profile 选择 | `g_use_nfw_profile`、`g_use_hernquist_aniso_profile`、`g_use_hernquist_numerical` |
| profile 参数 | `g_scale_radius_param`、`g_halo_mass_param`、`g_cutoff_factor_param`、`g_falloff_factor_param` |
| 各向异性 | `g_anisotropy_beta`、`g_use_om_profile`、`g_om_anisotropy_radius`、`g_use_numerical_isotropic` |
| 排序 | `g_defaultSortAlg`、`g_sort_columns_buffer`、benchmark 统计变量 |
| 轨迹 | `g_trajectories_buf`、`g_E_buf`、`g_L_buf`、`chosen` |
| 双精度缓冲 | `g_dbl_buf_R`、`g_dbl_buf_Vrad`、`g_dbl_buf_L`、`g_dbl_buf_Rank`、`g_dbl_buf_phi` |

这也是 `main()` 难读的根源之一：很多函数不是显式传参，而是读写全局状态。

## `main()` 的完整流程

`main()` 可以拆成 15 个阶段理解。

### 1. 启动、OpenMP 和 FFTW 初始化

程序一开始打印 NSphere banner。

如果编译时启用了 OpenMP：

1. 获取最大线程数和处理器数。
2. 开启 nested parallelism。
3. 在 macOS/Linux/Windows 上尝试检测 hybrid CPU 的 performance cores。
4. 设置 OpenMP 线程数。
5. 初始化 FFTW 线程支持。

如果没有 OpenMP：

1. 打印性能警告。
2. 如果不是 `--help`，故意等 2.5 秒，让用户看到警告。
3. 使用单线程 stub 函数。

### 2. 设置默认参数

默认参数包括：

| 参数 | 默认值 |
| --- | --- |
| `npts` | 100000 |
| `Ntimes` | 10000 |
| `tfinal_factor` | 5 |
| `nout` | 100 |
| `dtwrite` | 100 |
| `snapshot_block_size` | 100 |
| `tidal_fraction` | 0 |
| `method_select` | 1，即用户层面的默认方法 |
| `display_sort` | 1，即 parallel quadsort |

这里的 `method_select` 先是用户可见编号，后面会映射成内部编号。

### 3. 解析命令行参数

命令行解析是一个很长的 `for` + `if/else if` 链。它显式拒绝 `--option=value`，要求使用空格形式，比如 `--nparticles 100000`。

主要参数类别：

| 类别 | 参数 |
| --- | --- |
| 基础规模 | `--nparticles`、`--ntimesteps`、`--tfinal`、`--nout`、`--dtwrite`、`--snapshot-buffer` |
| 文件标识 | `--tag`、`--methodtag` |
| 保存模式 | `--save raw-data/psi-snaps/full-snaps/debug-energy/all` |
| 初始条件 | `--readinit`、`--writeinit` |
| restart | `--restart [force]`、`--sim-restart [check]`、`--restart-file` |
| extend | `--sim-extend`、`--extend-file` |
| profile | `--profile nfw/cored/hernquist`、`--scale-radius`、`--halo-mass`、`--cutoff-factor`、`--falloff-factor` |
| 各向异性 | `--aniso-beta`、`--aniso-factor`、`--aniso-betascale`、`--lvals-target` |
| SIDM | `--sidm`、`--sidm-mode serial/parallel`、`--sidm-kappa` |
| 随机性 | `--master-seed`、`--init-cond-seed`、`--sidm-seed`、`--load-seeds` |
| 算法选择 | `--method`、`--sort` |

解析过程中会立即检查一些互斥条件，例如：

1. `--readinit` 不能和 `--restart`、`--writeinit` 一起用。
2. `--writeinit` 不能和 `--restart`、`--readinit` 一起用。
3. `--sim-extend` 必须配 `--extend-file`。
4. `--sim-extend` 不能和 `--sim-restart`、`--restart`、`--readinit`、`--ftidal` 一起用。
5. `--aniso-beta` 只允许 Hernquist profile。
6. `--aniso-beta` 不能和 OM 各向异性一起用。

### 4. 方法编号映射

用户看到的 `--method 1..9` 会被映射成内部 `method_select`：

| 用户编号 | 内部编号 | 名称 |
| --- | --- | --- |
| 1 | 5 | Adaptive Leapfrog with Adaptive Levi-Civita |
| 2 | 4 | Full-Step Adaptive Leapfrog + Levi-Civita |
| 3 | 3 | Full-Step Adaptive Leapfrog |
| 4 | 6 | Yoshida 4th-Order |
| 5 | 8 | Adams-Bashforth 3rd-Order |
| 6 | 2 | Leapfrog, velocity half-step |
| 7 | 1 | Leapfrog, position half-step |
| 8 | 7 | Classic RK4 |
| 9 | 0 | Euler |

这张表很重要，因为模拟主循环里判断的是内部编号，而用户和帮助文本看到的是外部编号。

### 5. profile 参数归一化

`--profile` 决定走哪个初始条件路径：

1. `nfw`：`g_use_nfw_profile = 1`
2. `cored`：所有 Hernquist/NFW 标志关闭，走默认 cored branch
3. `hernquist`：
   - 如果使用 OM，则走 numerical Hernquist
   - 否则走 analytical anisotropic Hernquist

NFW 有独立默认值：

| 参数 | NFW 默认 |
| --- | --- |
| scale radius | `RC_NFW_DEFAULT = 1.18` |
| halo mass | `HALO_MASS_NFW = 1.15e9` |
| cutoff factor | 85 |
| falloff factor | 19 |

Cored/Hernquist 默认更多使用通用参数，比如 `HALO_MASS = 1.0e12` 和 `RC = 23`。

`g_active_halo_mass` 最终被设成当前 profile 的 halo mass，后续 `tdyn` 和 N-body force 都用它。

### 6. 输出后缀、目录、lastparams 和 seed

程序创建：

1. `data/`
2. `init/`

然后根据 tag、method tag、`npts`、`Ntimes`、`tfinal_factor` 构建 `g_file_suffix`。多数输出文件会通过 `get_suffixed_filename()` 加上这个后缀。

随后写：

1. `data/lastparams<suffix>.dat`
2. `data/lastparams.dat`，Unix 下是 symlink，Windows 下是 copy

随机种子优先级：

1. 显式 `--init-cond-seed` / `--sidm-seed`
2. `--master-seed` 派生：IC seed = master + 1，SIDM seed = master + 2
3. `--load-seeds` 从上次文件加载
4. 时间 + PID 生成新 seed

最终实际使用的种子会写到：

1. `data/last_initial_seed<suffix>.dat`
2. `data/last_sidm_seed<suffix>.dat`

### 7. 分配粒子数组和初始化 `phi`

`particles` 被分配为 `7 x npts_initial`。如果启用潮汐剥离，`npts_initial` 会比 `npts` 大，保证剥离外层后还剩 `npts` 个。

如果不是 `--readinit` 且不是 restart，则为每个粒子随机初始化 `phi`：

```text
particles[5][i] = cos(phi)
particles[6][i] = sin(phi)
```

SIDM 会用这个角度追踪切向速度方向。

### 8. profile-specific 初始条件生成

这是 `main()` 最大的前半部分。四条路径：

#### NFW 路径

NFW 分布使用带 cutoff 的 NFW-like density。

主要步骤：

1. 选择普通密度导数或 OM augmented density 导数。
2. 在 debug 模式下运行诊断循环，输出不同积分点/样条点下的 `massprofile`、`Psiprofile`、`f_of_E` 等。
3. 计算密度归一化 `nt_nfw`，使总质量等于目标 halo mass。
4. 构造 log-spaced 半径网格。
5. 逐点积分得到 `M(r)`，建立 `splinemass`。
6. 逐点积分得到相对势 `Psi(r)`，建立 `splinePsi`。
7. 建立反函数样条 `r(Psi)`，代码中实际使用 `-Psi` 以满足 GSL 单调性。
8. 用 Eddington 公式计算 `I(E)`。
9. 对 `I(E)` 建插值器，`f(E)` 实际通过 `dI/dE / (sqrt(8) * pi^2)` 得到。
10. 检查负的 `f(E)` 或 `f(Q)`，少量数值 artifact 会插值修正，大量负值会提示用户。
11. 用 `r(M)` 反变换采样半径。
12. 用 rejection sampling 采样速度。
13. 如果启用 OM，把 pseudo velocity `w` 映射成物理 `v_r` 和 `v_t`。
14. 写入 `particles[0..4]` 和粒子 ID。

#### Hernquist analytical anisotropic 路径

非 OM 的 Hernquist 走 `generate_ics_hernquist_anisotropic()`。它支持常数 `beta` 各向异性，核心分布函数在文件末尾 `df_hernquist_aniso()` 中实现，使用超几何函数并处理 `beta = 0.5`、`beta = -0.5` 的特殊情形。

流程上先用 `splines_only=1` 构造质量/势能/分布所需样条，再根据 `readinit`、restart 或正常生成决定填不填 `particles`。

#### Hernquist numerical 路径

当 Hernquist 使用 OM，或者用 `--aniso-factor inf` 触发 numerical isotropic pathway，会走 numerical Hernquist。

它和 NFW 更像：

1. 计算 Hernquist mass normalization。
2. 构造 `M(r)`、`Psi(r)`、`r(Psi)`。
3. 计算 `I(E)` 或 `I(Q)`。
4. 检查负分布函数。
5. 用 `r(M)` 和 rejection sampling 生成粒子。
6. OM 时做 pseudo velocity 到物理 velocity 的映射。

#### Cored Plummer-like 路径

默认非 NFW、非 Hernquist 时走 cored profile。

密度形状大致是：

```text
rho(r) ∝ 1 / (1 + (r/rc)^2)^3
```

它同样先跑诊断循环，再跑主理论计算：

1. 积分得到 mass normalization。
2. 构造 `M(r)` 样条。
3. 构造 `Psi(r)` 样条。
4. 构造 `r(Psi)` 样条。
5. 计算 `I(E)` / `I(Q)`。
6. 建 `f(E)` 插值器。
7. 用 `r(M)` 和 rejection sampling 采样粒子。
8. OM 时做 velocity transformation。

### 9. 初始条件保存、潮汐剥离、单位转换

如果 `--writeinit` 启用，会把剥离前的初始条件写到 `init/<file>`。

如果 `--ftidal > 0`：

1. 先对 `npts_initial` 个粒子按半径排序。
2. 只保留最内侧 `npts` 个。
3. 重新分配 `particles` 到最终大小。
4. 调用 `reassign_orig_ids_with_rank()` 重新分配 ID。

随后，如果不是 restart：

1. `particles[1]` 从速度大小变为径向速度：`v_r = |v| * mu`
2. `v_r` 从 km/s 转换成 kpc/Myr
3. `L` 也转换成 kpc^2/Myr

### 10. 写初始粒子、计算时间步

初始状态写到：

```text
data/particles<suffix>.dat
```

然后计算：

```text
tdyn = 1 / sqrt((VEL_CONV_SQ * G_CONST) * M / r_scale^3)
totaltime = tfinal_factor * tdyn
dt = totaltime / (Ntimes - 1)
```

这里 `r_scale` 根据 profile 选择 NFW scale radius、Hernquist scale radius 或 cored radius。

### 11. 轨迹、能量、低角动量粒子、restart 状态

这一阶段设置多个运行时跟踪系统：

1. 分配轨迹缓冲：`trajectories`、`single_trajectory`、`energy_and_angular_momentum_vs_time`、`lowest_l_trajectories` 或 `chosen_l_trajectories`。
2. 分配双精度滚动缓冲，用于保存最近 4 个 snapshot 的完整 double 精度。
3. 计算所有粒子的初始相对能量 `E_rel = Psi - KE` 和角动量，按 final rank ID 存入 `g_E_init_vals` 和 `g_L_init_vals`。
4. 选择 `nlowest` 个低角动量粒子，或选择最接近 `--lvals-target` 的粒子。
5. 保存/读取 `chosen_particles<suffix>.dat`，用于 restart/extend 保持轨迹连续。
6. 分配 `inverse_map`、SIDM scatter 状态数组和当前 timestep scatter count 数组。
7. 根据 `all_particle_data` 文件状态处理 `--sim-restart`。
8. 根据源文件状态处理 `--sim-extend`。
9. 初始化 `all_particle_data`、`all_particle_phi`、`all_particle_scatter_counts`、`all_particle_ids` 文件。
10. 估算磁盘空间，不足则退出，接近耗尽则询问用户。

### 12. 主模拟循环

核心循环形态：

```text
for (j = 0; j < Ntimes; j++) {
    如果 restart 且 j == 0，跳过重复写入
    根据 method_select 进入某个积分器分支
    每个分支通常：
        排序 particles by radius
        重建 inverse_map
        推进 r 和 v_r
        可选 SIDM scattering
        记录轨迹/能量
        到 dtwrite 间隔时写 block
        打印进度
}
```

每个积分器分支虽然写法不同，但几乎都遵循同样的数据依赖：

1. force 计算需要当前排序 rank。
2. 排序后必须更新 `inverse_map`。
3. 轨迹记录按粒子 ID 找当前位置，因此依赖 `inverse_map`。
4. `dtwrite` 到达时，把当前状态写入 block buffer。
5. block 满时 append 到磁盘。

#### 引力和角动量力

`gravitational_force()` 实际使用预计算常数：

```text
g_precalc_force_const = -(VEL_CONV_SQ * G_CONST) * (M / npts)
force = g_precalc_force_const * current_rank / r^2
```

`current_rank` 近似代表有多少粒子在该半径以内。

角动量项：

```text
effective_angular_force(r, L) = L^2 / r^3
```

因此径向方程大致为：

```text
dv_r/dt = -G M(<r) / r^2 + L^2 / r^3
dr/dt = v_r
```

#### SIDM step

每个积分器分支会调用 `handle_sidm_step()`。如果 `--sidm` 未启用，直接返回。

启用后：

1. parallel 模式且 OpenMP 可用：调用 `perform_sidm_scattering_parallel_graphcolor()`。
2. 否则：调用 `perform_sidm_scattering_serial()`。
3. 本 timestep 散射数累加到 `g_total_sidm_scatters`。

SIDM 的随机数不是普通全局 RNG，而是通过 seed、timestep、particle ID、call index 组合得到 deterministic random，方便 restart 后保持一致。

#### dtwrite 写盘

到达写入步时，代码会：

1. 用 `inverse_map` 按 particle ID 抽取当前 rank、R、Vrad、L。
2. 填 `Rank_block`、`R_block`、`Vrad_block`、`L_block`、`phi_block`、`ID_block`、`scatter_count_block`。
3. block 满时 append 到：
   - `all_particle_data<suffix>.dat`
   - `all_particle_phi<suffix>.dat`
   - `all_particle_scatter_counts<suffix>.dat`
   - `all_particle_ids<suffix>.dat`
4. 同时更新双精度滚动缓冲并刷到：
   - `double_buffer_all_particle_data<suffix>.dat`
   - `double_buffer_all_particle_phi<suffix>.dat`
5. 如果是输出 snapshot 点，计算系统总动能/势能/总能量。
6. 如果轨迹缓冲对应 block 满了，flush 轨迹文件。

### 13. 模拟结束后的直接输出

模拟循环结束后：

1. flush 最后不满一个 block 的粒子数据。
2. 写 `particlesfinal<suffix>.dat`。
3. 根据 profile 写理论曲线：
   - `massprofile<suffix>.dat`
   - `Psiprofile<suffix>.dat`
   - `density_profile<suffix>.dat`
   - `dpsi_dr<suffix>.dat`
   - `drho_dpsi<suffix>.dat`
   - `f_of_E<suffix>.dat`
   - `df_fixed_radius<suffix>.dat`

这些理论输出很多是绘图脚本后续使用的输入。

### 14. 直方图和 snapshot 后处理

程序重新读 `particles<suffix>.dat` 得到初始状态，并从当前 `particles` 得到最终状态，然后生成：

1. `2d_hist_initial<suffix>.dat`
2. `2d_hist_final<suffix>.dat`
3. `combined_histogram<suffix>.dat`

直方图范围不是固定常数，而是用初始+最终数据的 99 分位并乘 1.2 做 padding。

如果 `g_doAllParticleData` 开启，还会对 `all_particle_data` 中指定的 `noutsnaps` 个 snapshot 做并行后处理：

1. 调 `retrieve_all_particle_snapshot()` 读某个 snapshot。
2. 构造 unsorted arrays：按原粒子 ID 顺序。
3. 构造 `PartData` 并按半径排序。
4. 生成 sorted arrays。
5. 计算密度：先保证半径严格单调，必要时过滤和抽稀，再在 log grid 上做平滑。
6. 计算动态势 `PsiA(r)`。
7. 计算 sorted/unsorted energy。
8. 写：
   - `Rank_Mass_Rad_VRad_unsorted_tXXXXX<suffix>.dat`
   - `Rank_Mass_Rad_VRad_sorted_tXXXXX<suffix>.dat`

这部分本质上是后处理，不是时间推进模拟本身。

### 15. 清理

最后释放：

1. 轨迹缓冲
2. 初始/最终临时数组
3. mass/radius/Psi/E arrays
4. GSL workspace 和 RNG
5. `particles`
6. total energy diagnostics
7. debug energy diagnostics
8. GSL splines 和 accelerators
9. global block arrays
10. FFTW threads
11. persistent sort buffer
12. SIDM scatter state

如果使用 adaptive sort benchmark，还会打印排序算法统计。

## 辅助函数模块说明

### 文件名和日志

`get_suffixed_filename()` 根据 `g_file_suffix` 处理输出名：

1. 对 `.dat` 文件，把 suffix 插入扩展名前。
2. 对其他文件，把 suffix 追加到末尾。

`log_message()` 只有 `--log` 开启时才写 `log/nsphere.log`。日志里会带时间、level 和可选 suffix。

### 二进制 I/O

`fprintf_bin()` 和 `fscanf_bin()` 是自定义二进制 printf/scanf：

1. `%d` 写/读 `int`，4 字节。
2. `%f` / `%g` / `%e` 参数在 C varargs 中是 `double`，但写盘时压成 `float`，4 字节。
3. 格式字符串中的空格、换行、注释文本通常只是被忽略，不真的写文本分隔符。

`fprintf_bin_dbl()` 和 `fscanf_bin_dbl()` 类似，但浮点写成 `double`，用于双精度滚动缓冲。

因此很多 `.dat` 文件虽然扩展名像文本，但实际上是二进制。不要用普通文本编辑器或 `numpy.loadtxt()` 直接读。

### Restart/extend

restart 相关有三层：

1. `--restart`：主要用于跳过模拟，直接从已有 `all_particle_data` 做 post-processing。
2. `--sim-restart`：检测不完整模拟，备份并截断到一致 checkpoint，然后继续模拟。
3. `--sim-extend`：复制一个完成的模拟文件到新目标名，并从最后状态继续向后跑。

关键函数：

| 函数 | 作用 |
| --- | --- |
| `find_last_common_complete_snapshot()` | 找出所有相关输出文件共同完整的最后 snapshot |
| `truncate_files_to_snapshot()` | 把主数据、phi、ID、scatter、trajectory 文件截断到一致位置 |
| `load_particles_from_restart()` | 从双精度 buffer 或 float32 主文件恢复 `particles` |
| `parse_nsphere_filename()` | 从标准文件名中解析 N、Ntimes、tfinal |

双精度 buffer 只保存最近 4 个 snapshot。如果 restart 的 snapshot 不在里面，就退回 float32 文件，并打印精度下降警告。

### 初始条件和理论分布

核心函数/概念：

| 函数 | 作用 |
| --- | --- |
| `massintegrand()` | cored profile 的质量积分 integrand |
| `massintegrand_profile_nfwcutoff()` | NFW cutoff profile 的质量积分 integrand |
| `massintegrand_hernquist()` | Hernquist profile 的质量积分 integrand |
| `Psiintegrand()` | 势能积分 integrand，内部调用 profile-specific mass integrand |
| `fEintegrand()` | cored 的 Eddington integrand |
| `fEintegrand_nfw()` | NFW 的 Eddington/OM integrand |
| `fEintegrand_hernquist()` | Hernquist numerical integrand |
| `om_mu_integrand()` | 固定半径输出 `df_fixed_radius` 时，对 OM 模型积分 `mu` |
| `check_and_warn_negative_fQ()` | 检查并修复/提示负分布函数 |

Eddington 反演中，代码常把直接计算的积分叫 `I(E)`，真正的 distribution function 近似通过导数得到：

```text
f(E) = abs(dI/dE) / (sqrt(8) * pi^2)
```

对于 OM 模型，能量变量从 `E` 换成：

```text
Q = E - L^2 / (2 r_a^2)
```

并使用 augmented density `rho_Q(r) = rho(r) * (1 + r^2 / r_a^2)`。

### 时间积分器

基础力学工具：

| 函数 | 作用 |
| --- | --- |
| `gravitational_force()` | 根据排序 rank 和半径计算球壳引力 |
| `effective_angular_force()` | 计算 `L^2/r^3` |
| `doMicroLeapfrog()` | 固定微步数 leapfrog |
| `doAdaptiveFullLeap()` | coarse/fine 对比的自适应 full-step leapfrog |
| `doLeviCivitaLeapfrog()` | `rho = sqrt(r)` 的 Levi-Civita 正则化 |
| `doMicroLeviCivita()` | Levi-Civita 微步 |
| `doSingleTauStepAdaptiveLeviCivita()` | 单个 fictitious time step 的自适应版本 |
| `doAdaptiveFullLeviCivita()` | 完整 adaptive LC step |

Levi-Civita 正则化用于处理接近 `r = 0` 的 stiff orbit。它把 `r` 换成 `rho = sqrt(r)`，并使用 fictitious time `tau`，大致满足：

```text
dt_phys = rho^2 d_tau
```

这样近中心时物理时间步会自动变小。

### 排序

排序的目标是让 `particles` 按半径递增排列。可选算法：

| `--sort` | 内部名称 | 说明 |
| --- | --- | --- |
| 1 | `quadsort_parallel` | 默认，并行 quadsort |
| 2 | `quadsort` | 串行 quadsort |
| 3 | `insertion_parallel` | 并行 insertion sort |
| 4 | `insertion` | 串行 insertion sort |
| 5 | `parallel_radix` | 并行 radix sort |
| 6 | `benchmark_mode` | 每 1000 次排序比较算法并切换最快 |

`sort_particles_with_alg()` 的做法：

1. 把 `particles[component][i]` 转置成 `g_sort_columns_buffer[i][component]`。
2. 对这些 row 指针按第 0 列半径排序。
3. 再转置回原始 `particles` 布局。

这样做方便调用 `quadsort` 等对 row 的排序函数，但也使排序依赖全局持久缓冲。

### SIDM

SIDM 相关分成三层：

1. `handle_sidm_step()`：主循环调用入口，判断是否启用和选择 serial/parallel。
2. `perform_sidm_scattering_serial()`：串行 Monte Carlo 散射。
3. `perform_sidm_scattering_parallel_graphcolor()`：并行 graph coloring，避免同一轮里两个线程同时修改同一个粒子。

粒子真实 3D 速度由以下信息重建：

1. 径向速度 `v_r`
2. 角动量 `L` 给出的切向速度大小 `v_t = L / r`
3. `cos(phi)` 和 `sin(phi)` 给出的切向速度方向

散射时在两粒子的 center-of-mass frame 中随机化相对速度方向，再把新速度投影回 `v_r`、`L`、`phi`。

### 能量和 debug

有两种能量输出：

1. 总系统能量：`total_energy_vs_time<suffix>.dat`
2. 单个 debug 粒子的理论/动态能量对比：`debug_energy_compare<suffix>.dat`

总能量通过 `calculate_system_energies()` 计算。debug 粒子 ID 是硬编码的 `DEBUG_PARTICLE_ID = 4`。

## 主要输出文件

| 文件 | 内容 |
| --- | --- |
| `particles<suffix>.dat` | 初始粒子状态，二进制 |
| `particlesfinal<suffix>.dat` | 最终粒子状态，二进制 |
| `all_particle_data<suffix>.dat` | 主演化数据：rank、R、Vrad、L，每粒子 16 字节 |
| `all_particle_phi<suffix>.dat` | 每粒子 `phi`，float32 |
| `all_particle_ids<suffix>.dat` | 每粒子 ID，int32 |
| `all_particle_scatter_counts<suffix>.dat` | 每粒子 SIDM 散射计数，int32 |
| `double_buffer_all_particle_data<suffix>.dat` | 最近 4 个 snapshot 的 double 精度 R/Vrad/L/rank |
| `double_buffer_all_particle_phi<suffix>.dat` | 最近 4 个 snapshot 的 double 精度 phi |
| `trajectories<suffix>.dat` | 选定粒子半径、速度、mu 的 timestep 轨迹 |
| `single_trajectory<suffix>.dat` | 单粒子轨迹 |
| `energy_and_angular_momentum_vs_time<suffix>.dat` | 选定粒子的 E/L 演化 |
| `lowest_l_trajectories<suffix>.dat` | 最低角动量粒子轨迹 |
| `chosen_l_trajectories<suffix>.dat` | 使用 `--lvals-target` 时的目标 L 粒子轨迹 |
| `chosen_particles<suffix>.dat` | restart/extend 用的轨迹粒子选择 |
| `Rank_Mass_Rad_VRad_sorted_tXXXXX<suffix>.dat` | 后处理 snapshot，按半径排序 |
| `Rank_Mass_Rad_VRad_unsorted_tXXXXX<suffix>.dat` | 后处理 snapshot，按原粒子 ID 顺序 |
| `massprofile<suffix>.dat` | 理论质量剖面 |
| `Psiprofile<suffix>.dat` | 理论势能剖面 |
| `density_profile<suffix>.dat` | 理论密度剖面 |
| `f_of_E<suffix>.dat` | 分布函数 |
| `df_fixed_radius<suffix>.dat` | 固定半径速度分布 |
| `2d_hist_initial/final<suffix>.dat` | 初始/最终相空间二维直方图 |
| `combined_histogram<suffix>.dat` | 初始/最终径向分布对比 |
| `total_energy_vs_time<suffix>.dat` | 总能量诊断 |
| `debug_energy_compare<suffix>.dat` | debug 粒子能量诊断 |

## 最容易误解或出错的点

1. `particles[1]` 会从速度大小变成径向速度。读代码时必须看当前阶段。
2. 很多 `.dat` 是二进制，不是文本。
3. `main()` 里的 `method_select` 是内部编号，不等于用户传入的 `--method`。
4. 排序后数组下标表示半径 rank，粒子身份要通过 `particles[3]` 和 `inverse_map` 找回。
5. `normalization` 是全局变量，不同 profile/诊断循环都会使用，读写顺序很关键。
6. restart、sim-restart、sim-extend 是三套不同语义，不要混用。
7. `skip_file_writes` 和 `skip_simulation` 不是同一件事：一个控制写文件，一个控制模拟循环。
8. `g_doDebug` 默认是 1，所以 debug 相关理论文件和检查默认会发生。
9. OM 模型会把采样变量从 `E` 换成 `Q`，速度采样后还要做 pseudo velocity transformation。
10. SIDM 会修改速度和角动量，Adams-Bashforth 等多步法还要处理散射后的历史状态。
11. 排序算法用持久全局缓冲，不适合在多个独立模拟上下文中复用。
12. OpenMP 下有不少 `single`、`barrier`、`parallel for`，任何移动代码都要重新检查同步关系。

## 为什么 `main()` 会变得这么长

这个 `main()` 至少承担了以下职责：

1. CLI parser
2. 参数校验
3. profile 参数归一化
4. OpenMP/FFTW runtime 初始化
5. 文件后缀和目录管理
6. seed 管理
7. 初始条件理论计算
8. 初始条件采样
9. restart/extend 文件校验与复制
10. 轨迹粒子选择
11. 模拟主循环
12. 积分器 dispatch
13. SIDM dispatch
14. block buffer 写盘
15. energy diagnostics
16. 理论 profile 输出
17. histogram 输出
18. snapshot 后处理
19. 内存和 GSL/FFTW 清理
20. benchmark summary

这些职责彼此耦合在局部变量和全局变量上，所以 `main()` 只能不断向下膨胀。

## 建议的重构边界

如果后续要拆，不建议一上来改算法。先按“无行为变化”的边界抽函数/模块。

### 第一阶段：只拆配置和 I/O

可以先引入：

```c
typedef struct {
    int npts;
    int npts_initial;
    int Ntimes;
    int tfinal_factor;
    int nout;
    int noutsnaps;
    int dtwrite;
    int snapshot_block_size;
    double tidal_fraction;
    int method_display;
    int method_internal;
    int sort_display;
    char method_str[32];
    char custom_tag[256];
    char file_suffix[256];
} RunConfig;
```

然后拆出：

1. `parse_args(argc, argv, &cfg)`
2. `validate_config(&cfg)`
3. `derive_profile_config(&cfg, &profile)`
4. `init_output_dirs()`
5. `write_lastparams(&cfg)`
6. `init_seeds(&cfg)`

这一步风险最低，因为不碰数值算法。

### 第二阶段：封装粒子状态

把裸 `double **particles` 包成：

```c
typedef struct {
    int n;
    double *r;
    double *vrad;
    double *L;
    double *id;
    double *mu;
    double *cos_phi;
    double *sin_phi;
} ParticleSet;
```

短期可以仍然保留底层 7 行数组，只提供访问宏或 helper，避免一次性改太多排序代码。

### 第三阶段：拆 profile context

把 NFW/cored/Hernquist 的样条和数组归到一起：

```c
typedef struct {
    ProfileType type;
    double halo_mass;
    double scale_radius;
    double cutoff_factor;
    double rmax;
    double Psimin;
    double Psimax;
    double normalization;
    int num_points;
    double *radius;
    double *mass;
    double *Psi;
    double *Evalues;
    double *Ivalues;
    gsl_spline *mass_spline;
    gsl_spline *Psi_spline;
    gsl_spline *r_of_Psi_spline;
    gsl_interp *fofE_interp;
    gsl_interp_accel *mass_acc;
    gsl_interp_accel *Psi_acc;
    gsl_interp_accel *r_of_Psi_acc;
    gsl_interp_accel *fofE_acc;
} ProfileContext;
```

这样可以把 `normalization` 从危险全局变量变成 profile 私有字段。

### 第四阶段：积分器函数表

把主循环里的 `if (method_select == ...)` 换成函数表：

```c
typedef void (*IntegratorStepFn)(SimulationState *sim, const RunConfig *cfg);
```

每个积分器只负责一件事：推进一个 timestep。主循环负责公共流程：SIDM、轨迹记录、写盘、进度。

这能显著减少重复的“记录轨迹”和“更新 inverse_map”代码。

### 第五阶段：分离后处理

`Rank_Mass_Rad_VRad_*`、density smoothing、histogram、理论 profile 输出都可以搬到独立模块。它们依赖 `all_particle_data`，但不应该阻塞主模拟器的可读性。

## 推荐阅读顺序

如果只是想理解当前代码，不建议从第 1 行读到第 19252 行。更高效的顺序：

1. 先读本文档的“粒子数据布局”和“`main()` 的完整流程”。
2. 看 `main()` 开头的参数解析和方法映射。
3. 只选一个 profile 路径读，比如默认 NFW 或 cored。
4. 跳到主模拟循环，先理解默认 method 1 对应的内部 method 5。
5. 看 `gravitational_force()` 和 `effective_angular_force()`。
6. 看 `sort_particles_with_alg()`，理解排序如何改变粒子顺序。
7. 看 `append_all_particle_data_chunk_to_file()` 和 `retrieve_all_particle_snapshot()`，理解数据文件格式。
8. 如果关心 SIDM，再读 `handle_sidm_step()` 和两个 `perform_sidm_scattering_*()`。
9. 最后再读 snapshot 后处理和 density smoothing。

## 维护时的基本原则

1. 改 profile 生成逻辑时，必须同时检查理论 profile 输出、`df_fixed_radius` 和负 `f(E)` 检查。
2. 改粒子数组字段时，必须检查初始条件、restart、排序、SIDM、trajectory、post-processing 全链路。
3. 改排序时，必须保证 `particles[3]` 和 `inverse_map` 的语义不变。
4. 改写盘格式时，必须同步 Python 绘图脚本和 restart/extend 读取逻辑。
5. 改 timestep/restart 逻辑时，必须验证 `dtwrite`、`nwrite_total`、`restart_initial_nwrite`、`g_restart_initial_timestep` 的关系。
6. 改 SIDM 时，必须保持 deterministic random 与 timestep/particle ID 的关系，否则 restart 后不可复现。
7. 改 OpenMP 区域时，必须重新检查 `single`、`barrier`、critical I/O 和共享变量写入。

## 总结

`nsphere.c` 的物理核心并不难：球对称 N-body shell 模拟，按半径排序，用 rank 估计 enclosed mass，再积分径向运动。难点是这个文件同时承载了整个程序生命周期和大量研究功能，导致状态散落在全局变量、`main()` 局部变量和 profile-specific 分支里。

理解它的关键不是追每一行，而是抓住三条数据流：

1. `particles` 如何从初始条件生成，变成模拟态，再被排序和写盘。
2. `all_particle_data` 如何作为 restart、extend 和后处理的中心文件。
3. profile 理论对象 `M(r)`、`Psi(r)`、`f(E)/f(Q)` 如何服务于初始采样和诊断输出。

后续如果要让它变得可维护，首要目标应当是把 `main()` 从“所有事情的执行脚本”拆成“调度器”，并把配置、profile、粒子态、积分器、写盘、后处理分别封装起来。
