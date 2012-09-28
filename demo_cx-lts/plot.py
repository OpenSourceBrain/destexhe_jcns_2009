"""
Script to plot data from demo_cx05_N=500b_LTS simulation.

Usage: python plot.py <spike_times_file> <spike_counts_file> <vm_file> <neuron_index>

Andrew Davison, 2012
License: Modified BSD (see LICENSE.txt)
"""

import sys
import numpy
import matplotlib
import matplotlib.pyplot as plt

TSTOP = 5000 # ms
matplotlib.rcParams.update({
    'font.size': 9,
    'xtick.direction': 'out',
    'ytick.direction': 'out',})

def get_version():
    from mercurial import ui, hg
    from binascii import hexlify
    repo = hg.repository(ui.ui(), "..")
    ctx = repo.parents()[0]
    return hexlify(ctx.node()[:6])
__version__ = get_version()


spike_times_file, spike_counts_file, vm_file, neuron_id = sys.argv[1:]

fig = plt.figure(figsize=(8, 3))
fig.dpi = 120

# Plot spike times
with open(spike_times_file) as fp:
    data = numpy.loadtxt(fp)
ids, times = data.T

ax = fig.add_axes((0.1, 0.12, 0.6, 0.55), frameon=False)
ax.set_xlim([0, TSTOP])
ax.plot(times, ids, 'b.', markersize=0.2)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.set_xlabel("Time (ms)")
ax.set_ylabel("Cell number")

# Plot firing rate histogram
with open(spike_counts_file) as fp:
    fp.readline() # first line is metadata
    n_spikes = numpy.loadtxt(fp) * 1000.0/TSTOP  # calculate firing rate
    
bins = numpy.arange(0, 100, 5.0)
ax = fig.add_axes((0.78, 0.2, 0.2, 0.5), frameon=False)
plt.hist(n_spikes, bins)
# add the left and bottom axis lines back in
xmin, xmax = ax.get_xaxis().get_view_interval()
ymin, ymax = ax.get_yaxis().get_view_interval()
ax.add_artist(plt.Line2D((xmin, xmin), (ymin, ymax), color='black', linewidth=1))
ax.add_artist(plt.Line2D((xmin, xmax), (ymin, ymin), color='black', linewidth=1))
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.set_xlabel("Firing rate (Hz)")
ax.set_ylabel("Number of cells")

# Plot sample membrane potential trace
data = numpy.loadtxt(vm_file)
id, t, v = data[data[:, 0] == int(neuron_id)].T

ax = fig.add_axes((0.1, 0.73, 0.6, 0.25), frameon=False)
ax.set_xlim([0, TSTOP])
ax.plot(t, v, 'r', linewidth=0.8)
ax.xaxis.set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.set_xlabel("Time (ms)")
ax.set_ylabel("Vm (mV)")

fig.text(0.85, 0.95, __version__)
plt.savefig("demo_cx05_N=500b_LTS_%s.png" % __version__)
