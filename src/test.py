import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.interpolate import griddata

import json
from conformer import Conformer
from clustering import calc_pca, obtain_markers_from_cluster


confs = json.load(open("files/checkpoint.json"))
ensemble = [Conformer.load_raw(confs[i]) for i in confs]

pca_scores, clusters, colors, numbers, z = calc_pca(ensemble, ncluster=5)


fig = plt.figure()

gs = gridspec.GridSpec(2, 1, height_ratios=[50, 1])

ax = fig.add_subplot(gs[0])
color_axis = fig.add_subplot(gs[1])


x_ = pca_scores[:, 0]
y_ = pca_scores[:, 1]

resolution = 1500  # Risoluzione della griglia
xi = np.linspace(min(x_), max(x_), resolution)
yi = np.linspace(min(y_), max(y_), resolution)
xi, yi = np.meshgrid(xi, yi)
zi = griddata((x_, y_), z, (xi, yi), method="cubic", rescale=True)
# zi[zi <= 0] = 0.01

im = ax.pcolormesh(xi, yi, zi, shading="auto", cmap="coolwarm", alpha=0.75)
ax.contour(xi, yi, zi, "--", levels=10, colors="grey", linewidths=0.5, alpha=0.6)


# ax2 = plt.subplot(gs[1])
cbar = plt.colorbar(im, cax=color_axis, orientation="horizontal")
cbar.set_label("Potenziale")


for x, y, m, c, n in zip(
    pca_scores[:, 0],
    pca_scores[:, 1],
    np.array(list(map(obtain_markers_from_cluster, clusters))),
    colors,
    numbers,
):
    ax.scatter(x, y, marker=m, label=f"CONF {n}", c=c)
ax.set_xlabel("Principal Component 1")
ax.set_ylabel("Principal Component 2")
ax.set_title("Prova")


y = ax.get_ylim()
x = ax.get_xlim()

ax.vlines(0, y[0], y[1], "#353535", "--", alpha=0.2)
ax.hlines(0, x[0], x[1], "#353535", "--", alpha=0.2)

ax.set_xlim(x)
ax.set_ylim(y)
ax.legend(
    loc="upper left",
    bbox_to_anchor=(1.05, 1.0),
    fancybox=True,
    shadow=True,
    ncol=2,
    title="Conformers",
)

plt.tight_layout()


plt.show()
