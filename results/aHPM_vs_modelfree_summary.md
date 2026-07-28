# aHPM vs model-free comparison

Old run: `ml_results_aHPM_lag_K/` (log-lin μ paired with aHPM exit_lag_rate / N_max).
New run: `ml_results_modelfree/` (log-lin μ, Buchanan lag, q95 N_max — all Kinbiont).

## Parameter: `gr`

Top 5 old:

| compound | imp_old | imp_new | rank_old | rank_new | delta_rank |
|---|---|---|---|---|---|
| Glucose_(mM) | 0.7071 | 0.7071 | 1 | 1 | 0 |
| Na3C6H5O7/2H2O_(mM) | 0.07487 | 0.07487 | 2 | 2 | 0 |
| FeSO4/7H2O_(mM) | 0.03773 | 0.03773 | 3 | 3 | 0 |
| Isoleucine_(mM) | 0.03068 | 0.03068 | 4 | 4 | 0 |
| K2HPO4_(mM) | 0.029 | 0.029 | 5 | 5 | 0 |

Top 5 new:

| compound | imp_old | imp_new | rank_old | rank_new | delta_rank |
|---|---|---|---|---|---|
| Glucose_(mM) | 0.7071 | 0.7071 | 1 | 1 | 0 |
| Na3C6H5O7/2H2O_(mM) | 0.07487 | 0.07487 | 2 | 2 | 0 |
| FeSO4/7H2O_(mM) | 0.03773 | 0.03773 | 3 | 3 | 0 |
| Isoleucine_(mM) | 0.03068 | 0.03068 | 4 | 4 | 0 |
| K2HPO4_(mM) | 0.029 | 0.029 | 5 | 5 | 0 |

Biggest rank shifts (|Δrank| ≥ 5):

_None._

## Parameter: `lag`

Top 5 old:

| compound | imp_old | imp_new | rank_old | rank_new | delta_rank |
|---|---|---|---|---|---|
| KH2PO4_(mM) | 0.1197 | 0.1841 | 1 | 2 | 1 |
| K2HPO4_(mM) | 0.1175 | 0.04749 | 2 | 5 | 3 |
| Isoleucine_(mM) | 0.115 | 0.01546 | 3 | 10 | 7 |
| MgSO4/7H2O_(mM) | 0.1032 | 0.009111 | 4 | 18 | 14 |
| CaSO4/2H2O_(mM) | 0.06655 | 2.573e-05 | 5 | 41 | 36 |

Top 5 new:

| compound | imp_old | imp_new | rank_old | rank_new | delta_rank |
|---|---|---|---|---|---|
| Serine_(mM) | 0.01354 | 0.3496 | 15 | 1 | -14 |
| KH2PO4_(mM) | 0.1197 | 0.1841 | 1 | 2 | 1 |
| Leucine_(mM) | 0.006993 | 0.07766 | 23 | 3 | -20 |
| ZuSO4/7H2O_(mM) | 0.03545 | 0.04915 | 10 | 4 | -6 |
| K2HPO4_(mM) | 0.1175 | 0.04749 | 2 | 5 | 3 |

Biggest rank shifts (|Δrank| ≥ 5):

| compound | imp_old | imp_new | rank_old | rank_new | delta_rank |
|---|---|---|---|---|---|
| Thiamine/HCl_(mM) | 0.001893 | 0.01794 | 40 | 8 | -32 |
| KCl_(mM) | 0.003468 | 0.01931 | 33 | 7 | -26 |
| MgCl2/6H2O_(mM) | 0.001314 | 0.009542 | 41 | 17 | -24 |
| CaCl2/2H2O_(mM) | 0.002777 | 0.01141 | 36 | 13 | -23 |
| Leucine_(mM) | 0.006993 | 0.07766 | 23 | 3 | -20 |
| Asparagine/H2O_(mM) | 0.004872 | 0.01427 | 30 | 11 | -19 |
| Pyridoxine_(mM) | 0.001063 | 0.006957 | 43 | 24 | -19 |
| Tryptophan_(mM) | 0.002291 | 0.007424 | 37 | 21 | -16 |
| NaCl_(mM) | 0.000681 | 0.005594 | 44 | 29 | -15 |
| Serine_(mM) | 0.01354 | 0.3496 | 15 | 1 | -14 |
| Proline_(mM) | 0.006193 | 0.009972 | 26 | 16 | -10 |
| Valine_(mM) | 0.006302 | 0.01037 | 24 | 15 | -9 |
| ZuSO4/7H2O_(mM) | 0.03545 | 0.04915 | 10 | 4 | -6 |
| Alanine_(mM) | 0.002017 | 0.00253 | 39 | 33 | -6 |
| Cystine/HCl/H2O_(mM) | 0.01085 | 0.01155 | 18 | 12 | -6 |
| AminobenzoicAcid_(mM) | 0.01921 | 0.03741 | 11 | 6 | -5 |
| Glucose_(mM) | 0.006177 | 0.002843 | 27 | 32 | 5 |
| Arginine/HCl_(mM) | 0.01497 | 0.00794 | 14 | 20 | 6 |
| Isoleucine_(mM) | 0.115 | 0.01546 | 3 | 10 | 7 |
| Na2S2O3/5H2O_(mM) | 0.003278 | 5.477e-06 | 34 | 44 | 10 |
| Glutamine_(mM) | 0.004363 | 2.118e-05 | 31 | 42 | 11 |
| Na2MoO4/2H2O_(mM) | 0.008733 | 0.002463 | 21 | 34 | 13 |
| MgSO4/7H2O_(mM) | 0.1032 | 0.009111 | 4 | 18 | 14 |
| FeSO4/7H2O_(mM) | 0.01593 | 0.006401 | 12 | 26 | 14 |
| (NH4)2SO4_(mM) | 0.04754 | 0.00701 | 9 | 23 | 14 |
| CuSO4/5H2O_(mM) | 0.01129 | 0.002413 | 17 | 35 | 18 |
| Histidine/HCl/H2O_(mM) | 0.0548 | 0.005852 | 7 | 27 | 20 |
| Phenylalanine_(mM) | 0.05304 | 0.005695 | 8 | 28 | 20 |
| NH4Cl_(mM) | 0.009512 | 7.022e-06 | 20 | 43 | 23 |
| Na2HPO4_(mM) | 0.0573 | 0.001276 | 6 | 37 | 31 |
| CaSO4/2H2O_(mM) | 0.06655 | 2.573e-05 | 5 | 41 | 36 |

## Parameter: `N_max`

Top 5 old:

| compound | imp_old | imp_new | rank_old | rank_new | delta_rank |
|---|---|---|---|---|---|
| Valine_(mM) | 0.219 | 0.01344 | 1 | 5 | 4 |
| MgCl2/6H2O_(mM) | 0.2008 | 0.00432 | 2 | 8 | 6 |
| Riboflavin_(mM) | 0.08951 | 0.0001477 | 3 | 32 | 29 |
| MgSO4/7H2O_(mM) | 0.06448 | 0.001655 | 4 | 17 | 13 |
| Na3C6H5O7/2H2O_(mM) | 0.05954 | 0.03269 | 5 | 3 | -2 |

Top 5 new:

| compound | imp_old | imp_new | rank_old | rank_new | delta_rank |
|---|---|---|---|---|---|
| Glucose_(mM) | 0.03149 | 0.8499 | 8 | 1 | -7 |
| Leucine_(mM) | 0.009092 | 0.03393 | 17 | 2 | -15 |
| Na3C6H5O7/2H2O_(mM) | 0.05954 | 0.03269 | 5 | 3 | -2 |
| FeSO4/7H2O_(mM) | 0.01401 | 0.019 | 13 | 4 | -9 |
| Valine_(mM) | 0.219 | 0.01344 | 1 | 5 | 4 |

Biggest rank shifts (|Δrank| ≥ 5):

| compound | imp_old | imp_new | rank_old | rank_new | delta_rank |
|---|---|---|---|---|---|
| Lysine/HCl_(mM) | 0.002801 | 0.01232 | 37 | 6 | -31 |
| Threonine_(mM) | 0.001097 | 0.004147 | 39 | 9 | -30 |
| Arginine/HCl_(mM) | 0.0004489 | 0.002019 | 44 | 14 | -30 |
| Tyrosine_(mM) | 0.002878 | 0.0044 | 36 | 7 | -29 |
| Methionine_(mM) | 0.0006772 | 0.0009899 | 43 | 20 | -23 |
| ZuSO4/7H2O_(mM) | 0.003354 | 0.001466 | 33 | 18 | -15 |
| Leucine_(mM) | 0.009092 | 0.03393 | 17 | 2 | -15 |
| Na2MoO4/2H2O_(mM) | 0.004161 | 0.001867 | 28 | 15 | -13 |
| Glutamine_(mM) | 0.001635 | 0.0002164 | 38 | 28 | -10 |
| FeSO4/7H2O_(mM) | 0.01401 | 0.019 | 13 | 4 | -9 |
| Na2S2O3/5H2O_(mM) | 0.007677 | 0.00211 | 20 | 13 | -7 |
| KCl_(mM) | 0.000974 | 9.693e-05 | 41 | 34 | -7 |
| Cystine/HCl/H2O_(mM) | 0.001002 | 9.721e-05 | 40 | 33 | -7 |
| Glucose_(mM) | 0.03149 | 0.8499 | 8 | 1 | -7 |
| AminobenzoicAcid_(mM) | 0.005383 | 0.001757 | 22 | 16 | -6 |
| CuSO4/5H2O_(mM) | 0.004176 | 0.0006929 | 27 | 22 | -5 |
| Tryptophan_(mM) | 0.003127 | 4.26e-05 | 34 | 40 | 6 |
| MgCl2/6H2O_(mM) | 0.2008 | 0.00432 | 2 | 8 | 6 |
| KH2PO4_(mM) | 0.009997 | 0.0006686 | 16 | 23 | 7 |
| AsparticAcid_(mM) | 0.003594 | 5.131e-05 | 32 | 39 | 7 |
| H3BO3_(mM) | 0.004142 | 7.136e-05 | 29 | 37 | 8 |
| Phenylalanine_(mM) | 0.003064 | 8.058e-06 | 35 | 43 | 8 |
| Na2HPO4_(mM) | 0.01265 | 0.0006025 | 14 | 24 | 10 |
| (NH4)2SO4_(mM) | 0.01937 | 0.0007817 | 10 | 21 | 11 |
| CaCl2/2H2O_(mM) | 0.007729 | 0.0001833 | 19 | 30 | 11 |
| Thiamine/HCl_(mM) | 0.00488 | 8.299e-05 | 23 | 35 | 12 |
| MgSO4/7H2O_(mM) | 0.06448 | 0.001655 | 4 | 17 | 13 |
| CaSO4/2H2O_(mM) | 0.01055 | 0.0001516 | 15 | 31 | 16 |
| Glycine_(mM) | 0.004685 | 7.388e-06 | 25 | 44 | 19 |
| Asparagine/H2O_(mM) | 0.008096 | 1.604e-05 | 18 | 42 | 24 |
| FolicAcid_(mM) | 0.02084 | 8.016e-05 | 9 | 36 | 27 |
| Riboflavin_(mM) | 0.08951 | 0.0001477 | 3 | 32 | 29 |
| Histidine/HCl/H2O_(mM) | 0.05419 | 5.474e-05 | 6 | 38 | 32 |
