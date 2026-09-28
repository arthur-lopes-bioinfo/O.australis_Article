# Molecular Docking Results — Octopamine Receptor vs OaHLEα16

## 1. General Analysis Information

Molecular docking between the **octopamine receptor** and the enterotoxin **OaHLEα16** was performed using the **HADDOCK 2.4** server.

The run was successfully completed, including the *post-processing* stage.

The probable extracellular interaction region of the octopamine receptor was defined using the following residue ranges:

- **Pro78–Cys107**
- **Phe156–Thr170**
- **Pro194–Gly204**

These residues were used to guide docking toward the proposed extracellular interaction region of the receptor.

### Clustering Summary

- **Clustered structures:** 22
- **Number of clusters:** 4
- **Percentage of water-refined models included in clusters:** 11%
- **Maximum number of models considered for clustering by HADDOCK:** 200

HADDOCK ranks clusters primarily according to the **HADDOCK score**, with more negative values generally associated with more favorable solutions. The **Z-score** indicates how many standard deviations a cluster score differs from the mean score of all clusters, with more negative values indicating a relatively more favorable cluster.

## 2. Cluster Results

### Cluster 2 — Best Cluster

**Cluster 2** showed the best overall result, with a **HADDOCK score of -105.3 ± 9.2** and a **Z-score of -1.7**.

| Parameter                                     | Result            |
| --------------------------------------------- | ----------------- |
| HADDOCK score                                 | **-105.3 ± 9.2**  |
| Cluster size                                  | 6                 |
| RMSD from the overall lowest-energy structure | 29.2 ± 0.1 Å      |
| Van der Waals energy                          | -71.9 ± 5.4       |
| Electrostatic energy                          | -160.1 ± 18.3     |
| Desolvation energy                            | -32.3 ± 3.2       |
| Restraint violation energy                    | 309.1 ± 54.0      |
| Buried surface area (BSA)                     | 2328.2 ± 105.2 Å² |
| Z-score                                       | **-1.7**          |

Cluster 2 contains **six structures** and simultaneously exhibited the most negative HADDOCK score and the lowest Z-score among all identified clusters.

In addition, it showed a van der Waals energy of **-71.9 ± 5.4**, an electrostatic energy of **-160.1 ± 18.3**, and a desolvation energy of **-32.3 ± 3.2**.

The buried surface area at the interface was **2328.2 ± 105.2 Å²**, which was also the largest value observed among the four clusters.

Based on these criteria, **Cluster 2 was considered the most favorable docking solution for the octopamine receptor–OaHLEα16 complex**.

### Cluster 4

**Cluster 4** showed a HADDOCK score of **-64.0 ± 9.0** and a Z-score of **0.3**.

| Parameter                                     | Result           |
| --------------------------------------------- | ---------------- |
| HADDOCK score                                 | -64.0 ± 9.0      |
| Cluster size                                  | 5                |
| RMSD from the overall lowest-energy structure | 35.3 ± 0.1 Å     |
| Van der Waals energy                          | -63.2 ± 7.1      |
| Electrostatic energy                          | -114.7 ± 25.3    |
| Desolvation energy                            | -30.9 ± 3.9      |
| Restraint violation energy                    | 530.8 ± 32.7     |
| Buried surface area (BSA)                     | 2158.8 ± 97.5 Å² |
| Z-score                                       | 0.3              |

The cluster contains **five structures**. Although it showed favorable van der Waals, electrostatic, and desolvation energy contributions, its HADDOCK score was considerably less favorable than that observed for Cluster 2.

The restraint violation energy was also higher, reaching **530.8 ± 32.7**.

### Cluster 3

**Cluster 3** showed a HADDOCK score of **-62.9 ± 11.0** and a Z-score of **0.3**.

| Parameter                                     | Result            |
| --------------------------------------------- | ----------------- |
| HADDOCK score                                 | -62.9 ± 11.0      |
| Cluster size                                  | 5                 |
| RMSD from the overall lowest-energy structure | 37.3 ± 0.1 Å      |
| Van der Waals energy                          | -64.5 ± 8.1       |
| Electrostatic energy                          | -142.0 ± 24.1     |
| Desolvation energy                            | -21.7 ± 4.0       |
| Restraint violation energy                    | 517.4 ± 104.2     |
| Buried surface area (BSA)                     | 2203.0 ± 137.8 Å² |
| Z-score                                       | 0.3               |

The cluster contains **five structures** and showed an electrostatic energy of **-142.0 ± 24.1**, indicating an important contribution of electrostatic interactions to complex formation.

However, both the HADDOCK score and the Z-score were less favorable compared with Cluster 2.

### Cluster 1

**Cluster 1** exhibited the least favorable HADDOCK score among the four clusters, corresponding to **-47.5 ± 5.1**, with a Z-score of **1.0**.

| Parameter                                     | Result            |
| --------------------------------------------- | ----------------- |
| HADDOCK score                                 | -47.5 ± 5.1       |
| Cluster size                                  | 6                 |
| RMSD from the overall lowest-energy structure | 35.8 ± 0.7 Å      |
| Van der Waals energy                          | -44.7 ± 8.6       |
| Electrostatic energy                          | -163.9 ± 44.2     |
| Desolvation energy                            | -23.9 ± 6.6       |
| Restraint violation energy                    | 538.5 ± 25.9      |
| Buried surface area (BSA)                     | 1723.9 ± 188.6 Å² |
| Z-score                                       | 1.0               |

Despite showing a highly favorable electrostatic contribution (**-163.9 ± 44.2**), this cluster exhibited a weaker van der Waals contribution, a higher restraint violation energy, and a smaller buried surface area compared with Cluster 2.

## 3. Comparison Among Clusters

| Cluster | HADDOCK score    | Z-score  | VdW             | Electrostatic     | Desolvation     | Restraints       | BSA (Å²)           | Size |
| ------- | ---------------- | -------- | --------------- | ----------------- | --------------- | ---------------- | ------------------ | ---- |
| **2**   | **-105.3 ± 9.2** | **-1.7** | **-71.9 ± 5.4** | -160.1 ± 18.3     | **-32.3 ± 3.2** | **309.1 ± 54.0** | **2328.2 ± 105.2** | 6    |
| 4       | -64.0 ± 9.0      | 0.3      | -63.2 ± 7.1     | -114.7 ± 25.3     | -30.9 ± 3.9     | 530.8 ± 32.7     | 2158.8 ± 97.5      | 5    |
| 3       | -62.9 ± 11.0     | 0.3      | -64.5 ± 8.1     | -142.0 ± 24.1     | -21.7 ± 4.0     | 517.4 ± 104.2    | 2203.0 ± 137.8     | 5    |
| 1       | -47.5 ± 5.1      | 1.0      | -44.7 ± 8.6     | **-163.9 ± 44.2** | -23.9 ± 6.6     | 538.5 ± 25.9     | 1723.9 ± 188.6     | 6    |

## 4. Selection of the Best Conformation

Based on the parameters provided by HADDOCK, **Cluster 2** was selected as the most representative cluster for the interaction between the octopamine receptor and the OaHLEα16 enterotoxin.

This cluster exhibited:

- the **lowest HADDOCK score** (-105.3 ± 9.2);
- the **best Z-score** (-1.7);
- the **most favorable van der Waals energy** (-71.9 ± 5.4);
- the **most favorable desolvation energy** (-32.3 ± 3.2);
- the **lowest restraint violation energy** (309.1 ± 54.0);
- the **largest buried surface area** (2328.2 ± 105.2 Å²).

Therefore, the highest-ranked structure belonging to **Cluster 2** was considered the primary candidate conformation for subsequent structural analyses of the receptor–toxin interface.
