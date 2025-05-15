

import numpy as np
import tbplas as tb # Import the tbplas library
import matplotlib.pyplot as plt
from numpy.linalg import norm
from math import exp, sqrt

#Using class SK_Custom_SK1954_fixed, made by Max Sinner for evaluating SK integrals
from sk_custom_class import SK_Custom_SK1954_fixed

#--------------------
#0. Literature Values
#--------------------

# Slater-Koster parameters from Table II
sk_params_literature = {
    "Vpdσ": -2.619, "Vpdπ": -1.396,
    "Vddσ": -0.933, "Vddπ": -0.478, "Vddδ": -0.442,
    "Vppσ": 0.696, "Vppπ": 0.278,
}

sk_kwargs = {
    "v_sss": 0.0, "v_sps": 0.0,
    "v_pps": sk_params_literature.get("Vppσ", 0.0),
    "v_ppp": sk_params_literature.get("Vppπ", 0.0),
    "v_dds": sk_params_literature.get("Vddσ", 0.0),
    "v_ddp": sk_params_literature.get("Vddπ", 0.0),
    "v_ddd": sk_params_literature.get("Vddδ", 0.0),
    "v_pds": sk_params_literature.get("Vpdσ", 0.0),
    "v_pdp": sk_params_literature.get("Vpdπ", 0.0),
    "v_sds": 0.0
}

# Orbital energies - CORRECTED based on Cappelluti et al. 2013, Table II
e_mo_dz2 = -1.512    # Δ0
e_mo_d1 = -3.025     # Δ1 undetermined, using Δ2 as approximation
e_mo_d2 = -3.025     # Δ2
e_s_p = -1.276       # Δp (Used for px, py, AND pz)
e_s_s = -8.236       # Δz (Interpreted as S 's' energy)

orbital_energy_map = {}
# Mo
orbital_energy_map["Mo:dz2"] = e_mo_dz2
orbital_energy_map["Mo:dxz"] = e_mo_d1
orbital_energy_map["Mo:dyz"] = e_mo_d1
orbital_energy_map["Mo:dx2-y2"] = e_mo_d2
orbital_energy_map["Mo:dxy"] = e_mo_d2
# S top
orbital_energy_map["S_top:s"] = e_s_s
orbital_energy_map["S_top:px"] = e_s_p
orbital_energy_map["S_top:py"] = e_s_p
orbital_energy_map["S_top:pz"] = e_s_p # Corrected
# S bottom
orbital_energy_map["S_bot:s"] = e_s_s
orbital_energy_map["S_bot:px"] = e_s_p
orbital_energy_map["S_bot:py"] = e_s_p
orbital_energy_map["S_bot:pz"] = e_s_p # Corrected




#--------------------------
#1. Building MoS2 in tbplas
#--------------------------

# Lattice vectors for MoS2 (in Angstroms) - from Cappelluti et al. 2013 paper, Fig 1
vectors = np.array([
    [3.16, 0.0, 0.0],
    [-1.58, 2.736, 0.0],  # a/2 * (-1, √3, 0) where a=3.16
    [0.0, 0.0, 20.0],     # vacuum spacing
])
# Atomic fractional coordinates - from paper
coord_mo = np.array([
    [1/3, 2/3, 0.5],      # Mo atom
])
coord_s = np.array([
    [1/3, 2/3, 0.43],     # S atom (top, z=0.5 - 0.07)
    [1/3, 2/3, 0.57],     # S atom (bottom, z=0.5 + 0.07)
])

# Define orbitals - Mo with 5 d orbitals, S with 1 s and 3 p orbitals (13 total)
mo_orbital_types = ["dz2", "dxz", "dyz", "dx2-y2", "dxy"]
s_orbital_types = ["s", "px", "py", "pz"] # 4 orbitals per S atom

# Generate orbital coordinates and labels with atom identifiers
mo_orbital_coord = [row for row in coord_mo for _ in range(len(mo_orbital_types))]
mo_orbital_label = [f"Mo:{orb}" for orb in mo_orbital_types] * len(coord_mo)

s_top_orbital_coord = [coord_s[0] for _ in range(len(s_orbital_types))]
s_top_orbital_label = [f"S_top:{orb}" for orb in s_orbital_types]

s_bot_orbital_coord = [coord_s[1] for _ in range(len(s_orbital_types))]
s_bot_orbital_label = [f"S_bot:{orb}" for orb in s_orbital_types]

orbital_coord = mo_orbital_coord + s_top_orbital_coord + s_bot_orbital_coord
orbital_label = mo_orbital_label + s_top_orbital_label + s_bot_orbital_label

print("\nOrbitals added:")
for i in range(len(orbital_label)):
    print("    ",orbital_label[i], ":", orbital_coord[i])

orbital_energy_list = [orbital_energy_map[lbl] for lbl in orbital_label]



###CLASS:####

class mos2exp:
    def __init__(self,name):
        self.name = name
        #Creating Primitive Cell
        self.cell = tb.PrimitiveCell(lat_vec=vectors, unit=tb.ANG)
        for i in range(len(orbital_label)):
            self.cell.add_orbital_cart(orbital_coord[i], unit=tb.ANG, energy=orbital_energy_list[i], label=orbital_label[i])

    def addhoppings(self, cutoff_distance, amax, bmax, minhop):
        # Neighbor search
        neighbors = tb.find_neighbors(self.cell, a_max=amax, b_max=bmax, max_distance=cutoff_distance)


        print("Starting to add hoppings")
        # Add hopping terms using the custom SK class eval
        hopping_count = 0
        sk_instance = SK_Custom_SK1954_fixed() # Use the fixed custom class
        for term_idx, term in enumerate(neighbors):
            i, j = term.pair
            full_label_i = self.cell.get_orbital(i).label
            full_label_j = self.cell.get_orbital(j).label
            rij_array = np.array(term.rij)

            # Extract lm part of label (e.g., "dz2", "px", "dxz") - needed by custom class
            lm_i = full_label_i.split(":")[1]
            lm_j = full_label_j.split(":")[1]

            try:
                # Call the eval method of the CUSTOM SK instance
                hop = sk_instance.eval(r=rij_array, label_i=lm_i, label_j=lm_j, **sk_kwargs)

                if np.abs(hop) > minhop:
                    self.cell.add_hopping(term.rn, i, j, hop)
                    hopping_count += 1
            except Exception as e:
                print(f"\nERROR calculating/adding hopping for term {term_idx}:")
                print(f"  Pair: ({full_label_i}, {full_label_j}), Indices: ({i}, {j})")
                print(f"  Vector (rij): {term.rij}, Cell offset (rn): {term.rn}")
                print(f"  LM labels passed to eval: ({lm_i}, {lm_j})")
                print(f"  Error: {e}")

        print("Hoppings added")

    def calcbandstructure(self):
        # Calculate band structure
        k_points = np.array([
            [0.0, 0.0, 0.0],         # Γ
            [0.5, 0.0, 0.0],         # M
            [1/3, 1/3, 0.0],         # K
            [0.0, 0.0, 0.0]          # Γ
            ])
        k_label = ["Γ", "M", "K", "Γ"]
        k_path, k_idx = tb.gen_kpath(k_points, [50, 50, 50])
        k_len, bands = self.cell.calc_bands(k_path)
        plt.figure(figsize=(8, 6))
        vis = tb.Visualizer()
        figname = str(self.name) + "_bandstructure.png"
        vis.plot_bands(k_len, bands, k_idx=k_idx, k_label=k_label, fig_name=figname)





#########
# Main
# #######

"""
cutoff_distance = 3.5
amax=4
bmax = 4
minhop = 1e-6

classtest = mos2exp("classtest")
classtest.addhoppings(cutoff_distance,amax,bmax,minhop)
classtest.calcbandstructure()
"""