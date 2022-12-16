import numpy as np
import matplotlib.pyplot as plt

signal=np.load('signal.npy')
plt.plot(signal['time'], signal['value'])
plt.xlim(14950,15600)
plt.xlabel('Time (ns)')
plt.ylabel('Amplitude (V)')