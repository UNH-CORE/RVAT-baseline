from __future__ import division, print_function
from . import processing
import numpy as np
from nose.tools import assert_almost_equal


def test_find_amp_phase(plot=False):
    # Create a time series with a known amplitude and phase
    amp = 2.4
    phase = 1.5
    mean = 0.5
    npeaks = 3
    places = 2 # Decimal places the values should match
    noise_std = 0.01
    # Angle will start at zero first
    angle = np.linspace(0, 4*np.pi, num=300)
    angle_deg = np.rad2deg(angle)
    data = amp*np.cos(npeaks*(angle - phase)) + mean
    noise = np.random.normal(0, noise_std, len(data))
    data += noise
    amp_fit, phase_fit = processing.find_amp_phase(angle_deg, data, npeaks)
    if not plot:
        assert_almost_equal(amp, amp_fit, places=places)
        assert_almost_equal(phase, phase_fit, places=places)
    # Now check that we can fit a time series that doesn't start at zero
    data_old = data.copy()
    angle_old = angle.copy()
    angle_deg_old = angle_deg.copy()
    angle_deg += 247.05
    angle = np.deg2rad(angle_deg)
    data = amp*np.cos(npeaks*(angle - phase)) + mean + noise
    if plot:
        import matplotlib.pyplot as plt
        plt.plot(angle_deg_old, data_old, marker="o")
        plt.plot(angle_deg_old, amp_fit*np.cos(npeaks*(angle_old - phase_fit))
                 + mean, linewidth=4, alpha=0.6, color="red")
        plt.plot(angle_deg, data, marker="^", color="g")
    amp_fit, phase_fit = processing.find_amp_phase(angle_deg, data, npeaks)
    if not plot:
        assert_almost_equal(amp, amp_fit, places=places)
        assert_almost_equal(phase, phase_fit, places=places)
    else:
        plt.plot(angle_deg, amp_fit*np.cos(npeaks*(angle - phase_fit))
                 + mean, linewidth=4, alpha=0.6, color="purple")
        plt.show()
