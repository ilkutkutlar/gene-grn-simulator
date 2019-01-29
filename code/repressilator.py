from math import *

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Transcription-related values
ps_active = 0.5  # Promoter strength (active)
ps_repr = 5 * (10 ** -4)  # Promoter strength (repressed)
transcription_rate_active = ps_active * 60
transcription_rate_repr = ps_repr * 60

# Decay-related values
mRNA_half_life = 2
protein_half_life = 10
mRNA_decay_rate = log(2, e) / mRNA_half_life
protein_decay_rate = log(2, e) / protein_half_life

# Translation-related values
translation_efficiency = 19.97
translation_rate = translation_efficiency * mRNA_decay_rate

# Other values
n = 2  # Hill coefficient
Km = 40  # Activation coefficient

# Derived values

# Active transcription rate (rescaled?)
alpha = transcription_rate_active * translation_efficiency \
        * protein_half_life / (log(2, e) * Km)
# Repressed transcription rate (rescaled?)
alpha0 = transcription_rate_repr * translation_efficiency \
         * protein_half_life / (log(2, e) * Km)
# Translation rate
beta = protein_decay_rate / mRNA_decay_rate


def dy_dt(y, t):
    # Retrieve the current values of mRNAs

    # y represents the state of the network at time t.
    # It includes values for the concentrations of mRNA and protein.
    m_lacIi = y[0]
    m_tetRi = y[1]
    m_cli = y[2]

    p_lacIi = y[3]
    p_tetRi = y[4]
    p_cli = y[5]

    # Calculate changes in mRNA and protein using the equations
    # given in the original paper.
    m_lacI = -m_lacIi + (alpha / (1 + pow(p_cli, n))) + alpha0
    m_tetR = -m_tetRi + (alpha / (1 + pow(p_lacIi, n))) + alpha0
    m_cl = -m_cli + (alpha / (1 + pow(p_tetRi, n))) + alpha0

    p_lacI = beta * (m_lacI - p_lacIi)
    p_tetR = beta * (m_tetR - p_tetRi)
    p_cl = beta * (m_cl - p_cli)

    return [m_lacI, m_tetR, m_cl, p_lacI, p_tetR, p_cl]


# Initial values for proteins and mRNAs
# Values taken from the original paper
p_lacI0 = 10
p_tetR0 = 10
p_cl0 = 10

m_lacI0 = 100
m_tetR0 = 80
m_cl0 = 50

# Initial state
y0 = [m_lacI0, m_tetR0, m_cl0, p_lacI0, p_tetR0, p_cl0]

# time grid -> The time space for which a graph will be drawn
t = np.linspace(0, 40, 1000)

# solve the ODEs
soln = odeint(dy_dt, y0, t)

m_lacI = soln[:, 0]
m_tetR = soln[:, 1]
m_cl = soln[:, 2]
p_lacI = soln[:, 3]
p_tetR = soln[:, 4]
p_cl = soln[:, 5]

# plot results
plt.figure()

plt.plot(t, m_lacI, label='mLacI')
plt.plot(t, m_tetR, label='mTetR')
plt.plot(t, m_cl, label='mCl')

plt.xlabel('Time')
plt.ylabel('mRNA')
plt.title('mRNA amounts')
plt.legend(loc=0)

plt.draw()
plt.show()
