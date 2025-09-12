## Benchmarking leaflet classification methods

Here are the benchmarking results for various leaflet classification methods using `gorder` v1.2. Classification was performed for every analyzed frame.

### Benchmarked trajectories
- **Atomistic (CHARMM36):** 256 POPC lipids, 64 500 atoms, 10 000 trajectory frames  
- **Coarse-grained (Martini 3):** 512 POPC lipids, 16 800 beads, 10 000 trajectory frames  

### System configuration
- **CPU:** 8-core Intel Core i7-11700  
- **SSD:** Samsung 870 EVO  
- **OS:**  GNU/Linux Mint 20.2  

### Compiler versions
- `rustc` 1.89.0, `gcc` 9.4.0

### Important analysis options
- Radius of the cylinder for local leaflet classification was 2.0 nm.

### Other information
- Benchmarked using `hyperfine` with cold cache.
- Analyses for 'No assignment', 'Global', 'Individual', and 'FromFile' were run 5 times each.
- Analyses for 'Local' and 'Clustering' were run once.

---

### Benchmarking results

#### Atomistic system

| Method          | Threads | Analysis time [s] | Rel. to No |
|:---------------:|:-------:|:-----------------:|:----------:|
| No assignment   |    1    | 16.251 ± 0.032    |   100%     |
| Global          |    1    | 27.996 ± 0.040    |   172%     |
| Local           |    1    | ~1138             |  7000%     |
| Individual      |    1    | 17.250 ± 0.014    |   106%     |
| Clustering      |    1    | ~87               |   535%     |
| FromFile        |    1    | 17.679 ± 0.022    |   109%     |
|
| No assignment   |    8    | 5.858 ± 0.032     |   100%     |
| Global          |    8    | 9.683 ± 0.137     |   165%     |
| Local           |    8    | ~507              |  8650%     |
| Individual      |    8    | 6.287 ± 0.018     |   107%     |
| Clustering      |    8    | ~16               |   273%     |
| FromFile        |    8    | 6.738 ± 0.026     |   115%     |

---

#### Coarse-grained system

| Method         | Threads | Analysis time [s]  | Rel. to No |
|:--------------:|:-------:|:------------------:|:----------:|
| No assignment  |    1    | 4.720 ± 0.012      |   100%     |
| Global         |    1    | 7.146 ± 0.029      |   151%     |
| Local          |    1    | ~233               |  4940%     |
| Individual     |    1    | 4.979 ± 0.008      |   105%     |
| FromFile       |    1    | 5.886 ± 0.003      |   125%     |
|                                                            |
| No assignment  |    8    | 1.893 ± 0.010      |   100%     |
| Global         |    8    | 2.468 ± 0.007      |   130%     |
| Local          |    8    | ~39                |  2060%     |
| Individual     |    8    | 2.214 ± 0.005      |   117%     |
| FromFile       |    8    | 3.175 ± 0.011      |   168%     |