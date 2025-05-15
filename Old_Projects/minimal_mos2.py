import tbplas as tb
import numpy as np
import matplotlib.pyplot as plt

# Get cell fr MoS2
mos2_cell = tb.make_mos2_soc()

# Generate k-path
k_points = np.array([
    [0.0, 0.0, 0.0],    # Gamma
    [1./2, 0.0, 0.0],   # M
    [2./3, 1./3, 0.0],  # K
    [0.0, 0.0, 0.0],    # Gamma
])
k_label = ["G", "M", "K", "G"]
k_path, k_idx = tb.gen_kpath(k_points, [40, 40, 40])


# Calculate bands
k_len, bands = mos2_cell.calc_bands(k_path)


# Visualise bands
# vis = tb.Visualizer()
# vis.plot_bands(k_len, bands, k_idx, k_label)

num_bands = bands.shape[1]
for i in range(num_bands):
    plt.plot(k_len, bands[:, i], color="r", linewidth=1.0)
for idx in k_idx:
    plt.axvline(k_len[idx], color='k', linewidth=1.0)
plt.xlim((0, np.amax(k_len)))
plt.xticks(k_len[k_idx], k_label)
plt.xlabel("k (1/nm)")
plt.ylabel("Energy (eV)")
plt.savefig("mos2_bands.png", dpi=300)

