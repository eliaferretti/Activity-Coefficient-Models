# 1964 - WILSON Model

$$\log(\gamma_i) = -\log\left[ \sum_{j} x_j \Lambda_{ij} \right] + 1 - \sum_k\left[ \frac{ x_k \Lambda_{ki} }{ \sum_{j} x_j \Lambda_{kj} } \right]$$

$$ \Lambda_{ij} = \frac{v_j^L}{v_i^L}\exp\left(-\frac{a_{ij}-a_{ii}}{RT}\right) = \frac{v_j^L}{v_i^L}\exp\left(-\frac{\lambda_{ij}}{RT}\right)$$

$$ \Lambda_{ii} = 1 $$


# 1968 - NRTL Model (NonRandom Two-Liquid)

$$ \log(\gamma_i) = \frac{ \sum_{j} \tau_{ji} G_{ji} x_j}{\sum_{m} G_{mi} x_m } + \sum_{j}\left( \frac{x_jG_{ij}}{\sum_{m}x_mG_{mj}} \left( \tau_{ij}-\frac{\sum_{n}x_n\tau_{nj}G_{nj}}{\sum_{m} x_mG_{mj}}\right) \right)     $$

$$ \tau_{ij} = \frac{g_{ij}-g_{jj}}{RT} $$

$$ G_{ij} = \exp \left( -\alpha_{ij}\tau_{ij}\right) $$

$$ \tau_{ij} = 0,     G_{ij} = 1$$

# 1975 - UNIQUAC Model (Universal QuasiChemical)

$$ \log(\gamma_i) = \log(\gamma_i^C) + \log(\gamma_i^R) $$

$$ \log(\gamma_i^C) = \log(\frac{q_i}{\sum_jq_jx_j}) $$

$$ \log(\gamma_i^R) = $$
