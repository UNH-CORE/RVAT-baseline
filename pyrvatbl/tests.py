from __future__ import division, print_function
from . import processing
import numpy as np
from nose.tools import assert_almost_equal


def test_find_amp_phase():
    # Create a time series with a known amplitude and phase
    amp = 2.4
    phase = 1.5
    mean = 0.0
    npeaks = 3
    places = 2 # Decimal places the values should match
    noise_std = 0.01
    # Angle will start at zero first
    angle = np.linspace(0, 4*np.pi, num=200)
    angle_deg = np.rad2deg(angle)
    data = mean + amp*np.cos(npeaks*(angle - phase))
    noise = np.random.normal(0, noise_std, len(data))
    data += noise
    amp_fit, phase_fit = processing.find_amp_and_phase(angle_deg, data, npeaks)
    assert_almost_equal(amp, amp_fit, places=places)
    assert_almost_equal(phase, phase_fit, places=places)
