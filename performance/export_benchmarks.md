## Benchmarking the exporting

Here are the benchmarking results for analyses performed both without and with the `export` keyword for `gorder` v1.1 and v1.2.

### Benchmarked trajectory
- Atomistic, 256 POPC lipids, 64 500 atoms, 10 000 trajectory frames

### System configuration
- **CPU:** 8-core Intel Core i7-11700  
- **SSD:** Samsung 870 EVO  
- **OS:**  GNU/Linux Mint 20.2  

### Compiler versions
- `rustc` 1.89.0, `gcc` 9.4.0

### Important analysis options
- 'Leaflets' analyses: global leaflet classification
- 'Normals' analyses: dynamic membrane normals estimation

### Other information
- Benchmarked using `hyperfine` with cold cache.
- Each analysis was run 5 times.

---

### 'Leaflets' analyses

#### Frequency: every

| `gorder` version | Threads | Export | Analysis time [s]   |
|:----------------:|:-------:|:------:|:-------------------:|
|       1.1        |    1    |   no   |   28.039 ± 0.022    |
|       1.2        |    1    |   no   |   27.996 ± 0.040    |
|       1.2        |    1    |  yes   |   28.018 ± 0.019    |
|       1.1        |    8    |   no   |    9.605 ± 0.030    |
|       1.2        |    8    |   no   |    9.683 ± 0.137    |
|       1.2        |    8    |  yes   |    9.685 ± 0.073    |

---

#### Frequency: every 10th

| `gorder` version | Threads | Export | Analysis time [s]   |
|:----------------:|:-------:|:------:|:-------------------:|
|       1.1        |    1    |   no   |   18.268 ± 0.023    |
|       1.2        |    1    |   no   |   18.263 ± 0.008    |
|       1.2        |    1    |  yes   |   18.282 ± 0.013    |
|       1.1        |    8    |   no   |    7.330 ± 0.030    |
|       1.2        |    8    |   no   |    7.314 ± 0.014    |
|       1.2        |    8    |  yes   |    7.323 ± 0.013    |

---

#### Frequency: once

| `gorder` version | Threads | Export | Analysis time [s]   |
|:----------------:|:-------:|:------:|:-------------------:|
|       1.1        |    1    |   no   |   17.082 ± 0.017    |
|       1.2        |    1    |   no   |   17.075 ± 0.028    |
|       1.2        |    1    |  yes   |   17.058 ± 0.014    |
|       1.1        |    8    |   no   |    6.125 ± 0.019    |
|       1.2        |    8    |   no   |    6.110 ± 0.031    |
|       1.2        |    8    |  yes   |    6.125 ± 0.017    |

---

### 'Normals' analyses

| `gorder` version | Threads | Export | Analysis time [s]   |
|:----------------:|:-------:|:------:|:-------------------:|
|       1.1        |    1    |   no   |   24.114 ± 0.012    |
|       1.2        |    1    |   no   |   24.010 ± 0.014    |
|       1.2        |    1    |  yes   | **25.403 ± 0.149**  |
|       1.1        |    8    |   no   |    6.624 ± 0.036    |
|       1.2        |    8    |   no   |    6.679 ± 0.108    |
|       1.2        |    8    |  yes   | **8.168 ± 0.085**   |
