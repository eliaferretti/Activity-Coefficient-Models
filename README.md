# 1964 - WILSON Model

$$\log\left(\gamma_i\right) = -\log\left[ \sum_{j} x_j \Lambda_{ij} \right] + 1 - \sum_k\left[ \frac{ x_k \Lambda_{ki} }{ \sum_{j} x_j \Lambda_{kj} } \right]$$

$$ \Lambda_{ij} = \frac{v_j^L}{v_i^L}\exp\left(-\frac{a_{ij}-a_{ii}}{RT}\right) = \frac{v_j^L}{v_i^L}\exp\left(-\frac{\lambda_{ij}}{RT}\right)$$

$$ \Lambda_{ii} = 1 $$


# 1968 - NRTL Model (NonRandom Two-Liquid)

$$ \log\left(\gamma_i\right) = \frac{ \sum_{j} \tau_{ji} G_{ji} x_j}{\sum_{m} G_{mi} x_m } + \sum_{j}\left[ \frac{x_jG_{ij}}{\sum_{m}x_mG_{mj}} \left( \tau_{ij}-\frac{\sum_{n}x_n\tau_{nj}G_{nj}}{\sum_{m} x_mG_{mj}}\right) \right]     $$

$$ \tau_{ij} = \frac{g_{ij}-g_{jj}}{RT} $$

$$ G_{ij} = \exp \left( -\alpha_{ij}\tau_{ij}\right) $$

$$ \tau_{ij} = 0,     G_{ij} = 1$$

# 1975 - UNIQUAC Model (Universal QuasiChemical)

$$ \log\left(\gamma_i\right) = \log\left(\gamma_i^C\right) + \log\left(\gamma_i^R\right) $$

$$ \log\left(\gamma_i^C\right) = \log\left(\frac{q_i}{\sum_jx_jq_j}\right) + \frac{z}{2}\log\left(\frac{q_i\sum_jx_jr_j}{r_i\sum_jx_jq_j}\right) + I_i - \frac{q_i}{\sum_jx_jq_j}\sum_jx_jI_j $$

$$ \log\left(\gamma_i^R\right) = q_i \left[1-\log\left(\sum_j\theta_j\tau_{ji}\right) - \sum_j\left(\right)\right]$$
