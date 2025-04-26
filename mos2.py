# Start of script
#Building MoS2 in TBPLas
#Version of Max Sinner

print("DEBUG: Script starting...")

import numpy as np
import tbplas as tb # Import the tbplas library
import matplotlib.pyplot as plt
from numpy.linalg import norm
from math import exp, sqrt

# --- Definitions (Constants, Coordinates, Orbitals, Energies, SK Params) ---
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

print(f"Mo orbitals{len(mo_orbital_coord)}")

s_top_orbital_coord = [coord_s[0] for _ in range(len(s_orbital_types))]
s_top_orbital_label = [f"S_top:{orb}" for orb in s_orbital_types]

s_bot_orbital_coord = [coord_s[1] for _ in range(len(s_orbital_types))]
s_bot_orbital_label = [f"S_bot:{orb}" for orb in s_orbital_types]

orbital_coord = mo_orbital_coord + s_top_orbital_coord + s_bot_orbital_coord

print(f"length or orbital vectors:{len(orbital_coord)}")
orbital_label = mo_orbital_label + s_top_orbital_label + s_bot_orbital_label
num_orbitals = len(orbital_label)
print(f"DEBUG: Orbital setup complete. Total orbitals: {num_orbitals}")

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

orbital_energy_list = [orbital_energy_map[lbl] for lbl in orbital_label]
print(f"DEBUG: Corrected on-site energies assigned (S pz uses Δp={e_s_p} eV, S s uses Δz={e_s_s} eV).")

# Slater-Koster parameters from Table II
sk_params_literature = {
    "Vpdσ": -2.619, "Vpdπ": -1.396,
    "Vddσ": -0.933, "Vddπ": -0.478, "Vddδ": -0.442,
    "Vppσ": 0.696, "Vppπ": 0.278,
}
print("DEBUG: Hopping parameters defined.")


# --- Custom SK Class Definition with fixed d-d evaluation ---
class SK_Custom_SK1954_fixed:
    """
    Custom Slater-Koster class implementing s-s, s-p, p-p, s-d, p-d, d-d interactions
    based on Slater & Koster, Phys. Rev. 94, 1498 (1954), Table I.
    Includes fixes for d-d combination handling.
    Accepts standard orbital labels ('s', 'px', 'py', 'pz', 'dxy', 'dyz', 'dxz', 'dx2-y2', 'dz2').
    """
    def __init__(self):
        pass

    # Map python labels to S&K Table I labels internally
    sk_label_map = {
        's': 's', 'px': 'x', 'py': 'y', 'pz': 'z',
        'dxy': 'xy', 'dyz': 'yz', 'dxz': 'zx', # Convert dxz to zx for internal use
        'dx2-y2': 'x2-y2', 'dz2': '3z2-r2'
    }
    # Map S&K integrals to python kwargs
    sk_integral_map = {
        '(ssσ)': 'v_sss', '(spσ)': 'v_sps', '(sdσ)': 'v_sds',
        '(ppσ)': 'v_pps', '(ppπ)': 'v_ppp',
        '(pdσ)': 'v_pds', '(pdπ)': 'v_pdp',
        '(ddσ)': 'v_dds', '(ddπ)': 'v_ddp', '(ddδ)': 'v_ddd'
    }

    def eval(self, r, label_i, label_j, **kwargs):
        """
        Evaluate hopping integral using Slater-Koster scheme.
        Handles orbital order and uses formulas from S&K Table I.
        """
        r_vec = np.array(r); r_norm = norm(r_vec)
        if r_norm < 1e-10: return 0.0

        sk_i = self.sk_label_map.get(label_i)
        sk_j = self.sk_label_map.get(label_j)
        if sk_i is None or sk_j is None:
             print(f"DEBUG WARNING: Unknown label passed to custom eval: {label_i} or {label_j}")
             return 0.0

        type_i = label_i[0]; type_j = label_j[0]
        l, m, n = r_vec / r_norm
        v = {k: kwargs.get(v, 0.0) for k, v in self.sk_integral_map.items()}

        # --- Dispatch based on type combination ---
        if type_i == 's' and type_j == 's': return v['(ssσ)']
        elif type_i == 's' and type_j == 'p':
            if sk_j == 'x': return l * v['(spσ)']
            elif sk_j == 'y': return m * v['(spσ)']
            else: return n * v['(spσ)']
        elif type_i == 'p' and type_j == 's':
            if sk_i == 'x': return (-l) * v['(spσ)']
            elif sk_i == 'y': return (-m) * v['(spσ)']
            else: return (-n) * v['(spσ)']
        elif type_i == 'p' and type_j == 'p':
            if sk_i == 'x' and sk_j == 'x': return l**2 * v['(ppσ)'] + (1 - l**2) * v['(ppπ)']
            elif sk_i == 'y' and sk_j == 'y': return m**2 * v['(ppσ)'] + (1 - m**2) * v['(ppπ)']
            elif sk_i == 'z' and sk_j == 'z': return n**2 * v['(ppσ)'] + (1 - n**2) * v['(ppπ)']
            elif (sk_i == 'x' and sk_j == 'y') or (sk_i == 'y' and sk_j == 'x'):
                return l * m * v['(ppσ)'] - l * m * v['(ppπ)']
            elif (sk_i == 'y' and sk_j == 'z') or (sk_i == 'z' and sk_j == 'y'):
                return m * n * v['(ppσ)'] - m * n * v['(ppπ)']
            elif (sk_i == 'z' and sk_j == 'x') or (sk_i == 'x' and sk_j == 'z'):
                return n * l * v['(ppσ)'] - n * l * v['(ppπ)']
            else: print(f"DEBUG WARNING: Unhandled p-p: {label_i}({sk_i}), {label_j}({sk_j})"); return 0.0
        elif type_i == 's' and type_j == 'd': return self._eval_sd_sk1954(l, m, n, sk_i, sk_j, v)
        elif type_i == 'd' and type_j == 's': return self._eval_sd_sk1954(-l, -m, -n, sk_j, sk_i, v)
        elif type_i == 'p' and type_j == 'd': return self._eval_pd_sk1954(l, m, n, sk_i, sk_j, v)
        elif type_i == 'd' and type_j == 'p': return self._eval_pd_sk1954(-l, -m, -n, sk_j, sk_i, v)
        elif type_i == 'd' and type_j == 'd':
             if sk_i > sk_j: sk_i_calc, sk_j_calc = sk_j, sk_i
             else: sk_i_calc, sk_j_calc = sk_i, sk_j
             return self._eval_dd_sk1954(l, m, n, sk_i_calc, sk_j_calc, v)
        else: print(f"DEBUG WARNING: Unhandled type combo: {type_i}, {type_j}"); return 0.0

    def _eval_sd_sk1954(self, l, m, n, sk_s, sk_d, v):
        """ Helper for s-d terms from S&K 1954 Table I. Assumes sk_s = 's' """
        if abs(v['(sdσ)']) < 1e-10: return 0.0
        if sk_d == 'xy': return sqrt(3) * l * m * v['(sdσ)']
        elif sk_d == 'yz': return sqrt(3) * m * n * v['(sdσ)']
        elif sk_d == 'zx': return sqrt(3) * n * l * v['(sdσ)']
        elif sk_d == 'x2-y2': return sqrt(3)/2 * (l**2 - m**2) * v['(sdσ)']
        elif sk_d == '3z2-r2': return (n**2 - (l**2 + m**2)/2) * v['(sdσ)']
        else: print(f"DEBUG WARNING: Unknown s-d combo: {sk_s}, {sk_d}"); return 0.0

    def _eval_pd_sk1954(self, l, m, n, sk_p, sk_d, v):
        """ Helper for p-d terms from S&K 1954 Table I. Assumes sk_p is p, sk_d is d"""
        if abs(v['(pdσ)']) < 1e-10 and abs(v['(pdπ)']) < 1e-10: return 0.0
        if sk_d == 'xy':
            if sk_p == 'x': return sqrt(3)*l**2*m * v['(pdσ)'] + m*(1-2*l**2) * v['(pdπ)']
            elif sk_p == 'y': return sqrt(3)*l*m**2 * v['(pdσ)'] + l*(1-2*m**2) * v['(pdπ)']
            elif sk_p == 'z': return sqrt(3)*l*m*n * v['(pdσ)'] - 2*l*m*n * v['(pdπ)']
        elif sk_d == 'yz':
            if sk_p == 'x': return sqrt(3)*l*m*n * v['(pdσ)'] - 2*l*m*n * v['(pdπ)']
            elif sk_p == 'y': return sqrt(3)*m**2*n * v['(pdσ)'] + n*(1-2*m**2) * v['(pdπ)']
            elif sk_p == 'z': return sqrt(3)*m*n**2 * v['(pdσ)'] + m*(1-2*n**2) * v['(pdπ)']
        elif sk_d == 'zx':
            if sk_p == 'x': return sqrt(3)*l**2*n * v['(pdσ)'] + n*(1-2*l**2) * v['(pdπ)']
            elif sk_p == 'y': return sqrt(3)*l*m*n * v['(pdσ)'] - 2*l*m*n * v['(pdπ)']
            elif sk_p == 'z': return sqrt(3)*l*n**2 * v['(pdσ)'] + l*(1-2*n**2) * v['(pdπ)']
        elif sk_d == 'x2-y2':
            if sk_p == 'x': return sqrt(3)/2*l*(l**2-m**2) * v['(pdσ)'] + l*(1-(l**2-m**2)) * v['(pdπ)']
            elif sk_p == 'y': return sqrt(3)/2*m*(l**2-m**2) * v['(pdσ)'] - m*(1+(l**2-m**2)) * v['(pdπ)']
            elif sk_p == 'z': return sqrt(3)/2*n*(l**2-m**2) * v['(pdσ)'] - n*(l**2-m**2) * v['(pdπ)']
        elif sk_d == '3z2-r2':
            if sk_p == 'x': return l*(n**2 - (l**2+m**2)/2) * v['(pdσ)'] - sqrt(3)*l*n**2 * v['(pdπ)']
            elif sk_p == 'y': return m*(n**2 - (l**2+m**2)/2) * v['(pdσ)'] - sqrt(3)*m*n**2 * v['(pdπ)']
            elif sk_p == 'z': return n*(n**2 - (l**2+m**2)/2) * v['(pdσ)'] + sqrt(3)*n*(l**2+m**2) * v['(pdπ)']
        else: print(f"DEBUG WARNING: Unknown p-d combo: {sk_p}, {sk_d}"); return 0.0

    def _eval_dd_sk1954(self, l, m, n, sk_i, sk_j, v):
        """ Fully implemented d-d terms from S&K 1954 Table I. Assumes sk_i <= sk_j """
        if abs(v['(ddσ)']) < 1e-10 and abs(v['(ddπ)']) < 1e-10 and abs(v['(ddδ)']) < 1e-10: return 0.0

        # --- d-d: t2g - t2g ---
        if sk_i == 'xy' and sk_j == 'xy':
             # Using formula from user's original script for delta term, S&K paper unclear
             return 3*l**2*m**2*v['(ddσ)'] + (l**2+m**2-4*l**2*m**2)*v['(ddπ)'] + (n**2 + (l**2 - m**2)**2)*v['(ddδ)']
        elif sk_i == 'xy' and sk_j == 'yz': return 3*l*m**2*n*v['(ddσ)'] + l*n*(1-4*m**2)*v['(ddπ)'] + l*n*(m**2-1)*v['(ddδ)']
        elif sk_i == 'xy' and sk_j == 'zx': return 3*l**2*m*n*v['(ddσ)'] + m*n*(1-4*l**2)*v['(ddπ)'] + m*n*(l**2-1)*v['(ddδ)']
        elif sk_i == 'yz' and sk_j == 'yz': return 3*m**2*n**2*v['(ddσ)'] + (m**2+n**2-4*m**2*n**2)*v['(ddπ)'] + (l**2 + (m**2 - n**2)**2)*v['(ddδ)']
        elif sk_i == 'yz' and sk_j == 'zx': return 3*m*n**2*l*v['(ddσ)'] + m*l*(1-4*n**2)*v['(ddπ)'] + m*l*(n**2-1)*v['(ddδ)']
        elif sk_i == 'zx' and sk_j == 'zx': return 3*n**2*l**2*v['(ddσ)'] + (n**2+l**2-4*n**2*l**2)*v['(ddπ)'] + (m**2 + (n**2 - l**2)**2)*v['(ddδ)']

        # --- d-d: eg - eg ---
        elif sk_i == 'x2-y2' and sk_j == 'x2-y2': return 3/4*(l**2-m**2)**2*v['(ddσ)'] + (l**2+m**2-(l**2-m**2)**2)*v['(ddπ)'] + (n**2+1/4*(l**2-m**2)**2)*v['(ddδ)']
        elif sk_i == '3z2-r2' and sk_j == '3z2-r2': return (n**2-(l**2+m**2)/2)**2*v['(ddσ)'] + 3*n**2*(l**2+m**2)*v['(ddπ)'] + 3/4*(l**2+m**2)**2*v['(ddδ)']
        # Mixed eg-eg: sorted sk_i = '3z2-r2', sk_j = 'x2-y2'
        elif sk_i == '3z2-r2' and sk_j == 'x2-y2': return sqrt(3)/2*(l**2-m**2)*(n**2-(l**2+m**2)/2)*v['(ddσ)'] + sqrt(3)*n**2*(m**2-l**2)*v['(ddπ)'] + sqrt(3)/4*(1+n**2)*(l**2-m**2)*v['(ddδ)']

        # --- d-d: t2g - eg --- (Sorted: eg comes first alphabetically: '3z2...', 'x2...')
        # Case 1: sk_i = '3z2-r2'
        elif sk_i == '3z2-r2' and sk_j == 'xy': return sqrt(3)*l*m*(n**2-(l**2+m**2)/2)*v['(ddσ)'] - 2*sqrt(3)*l*m*n**2*v['(ddπ)'] + sqrt(3)/2*l*m*(1+n**2)*v['(ddδ)']
        elif sk_i == '3z2-r2' and sk_j == 'yz': return sqrt(3)*m*n*(n**2-(l**2+m**2)/2)*v['(ddσ)'] + sqrt(3)*m*n*(l**2+m**2-n**2)*v['(ddπ)'] - sqrt(3)/2*m*n*(l**2+m**2)*v['(ddδ)']
        elif sk_i == '3z2-r2' and sk_j == 'zx': return sqrt(3)*n*l*(n**2-(l**2+m**2)/2)*v['(ddσ)'] + sqrt(3)*n*l*(l**2+m**2-n**2)*v['(ddπ)'] - sqrt(3)/2*n*l*(l**2+m**2)*v['(ddδ)']
        # Case 2: sk_i = 'x2-y2'
        elif sk_i == 'x2-y2' and sk_j == 'xy': return 3/2*l*m*(l**2-m**2)*v['(ddσ)'] + 2*l*m*(m**2-l**2)*v['(ddπ)'] + 1/2*l*m*(l**2-m**2)*v['(ddδ)']
        elif sk_i == 'x2-y2' and sk_j == 'yz': return 3/2*m*n*(l**2-m**2)*v['(ddσ)'] - m*n*(1+2*(l**2-m**2))*v['(ddπ)'] + m*n*(1+1/2*(l**2-m**2))*v['(ddδ)']
        elif sk_i == 'x2-y2' and sk_j == 'zx': return 3/2*n*l*(l**2-m**2)*v['(ddσ)'] + n*l*(1-2*(l**2-m**2))*v['(ddπ)'] - n*l*(1-1/2*(l**2-m**2))*v['(ddδ)']

        # --- Fallback ---
        else:
            # This should now only catch combinations not listed in S&K Table I or logic errors
            print(f"DEBUG WARNING: Fallback in _eval_dd: {sk_i}, {sk_j}")
            return 0.0

print("DEBUG: SK_Custom_SK1954_fixed class defined.")

# --- Main Script ---

# Create primitive cell
print("DEBUG: Creating PrimitiveCell...")
cell = tb.PrimitiveCell(lat_vec=vectors, unit=tb.ANG)
print("DEBUG: PrimitiveCell created.")

# Add orbitals with energies
print("DEBUG: Adding orbitals...")
for i, label in enumerate(orbital_label):
    coord = orbital_coord[i]
    energy = orbital_energy_list[i]
    cell.add_orbital(coord, energy=energy, label=label)
print(f"DEBUG: Finished adding {len(cell.orbitals)} orbitals.")

# Hopping calculation setup using CUSTOM SK class
sk_instance = SK_Custom_SK1954_fixed() # Use the fixed custom class
print(f"DEBUG: Using CUSTOM SK_Custom_SK1954_fixed instance.")

# Prepare the keyword arguments dictionary
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
print(f"DEBUG: Prepared SK kwargs: {sk_kwargs}")

# Neighbor search
cutoff_distance = 3.5 # Angstroms
print(f"DEBUG: Finding neighbors within {cutoff_distance} Å...")
neighbors = tb.find_neighbors(cell, a_max=2, b_max=2, max_distance=cutoff_distance)
print(f"DEBUG: Found {len(neighbors)} neighbor pairs.")

# Add hopping terms using the custom SK class eval
hopping_count = 0
print("DEBUG: Adding hopping terms using custom SK_Custom_SK1954_fixed.eval...")
for term_idx, term in enumerate(neighbors):
    i, j = term.pair
    full_label_i = cell.get_orbital(i).label
    full_label_j = cell.get_orbital(j).label
    rij_array = np.array(term.rij)

    # Extract lm part of label (e.g., "dz2", "px", "dxz") - needed by custom class
    lm_i = full_label_i.split(":")[1]
    lm_j = full_label_j.split(":")[1]

    try:
        # Call the eval method of the CUSTOM SK instance
        hop = sk_instance.eval(r=rij_array, label_i=lm_i, label_j=lm_j, **sk_kwargs)

        if np.abs(hop) > 1e-6:
            cell.add_hopping(term.rn, i, j, hop)
            hopping_count += 1
    except Exception as e:
        print(f"\nERROR calculating/adding hopping for term {term_idx}:")
        print(f"  Pair: ({full_label_i}, {full_label_j}), Indices: ({i}, {j})")
        print(f"  Vector (rij): {term.rij}, Cell offset (rn): {term.rn}")
        print(f"  LM labels passed to eval: ({lm_i}, {lm_j})")
        print(f"  Error: {e}")
        # break # Uncomment to stop on first error

print(f"DEBUG: Finished adding {hopping_count} non-zero hopping terms.") # Expecting more non-zero terms now

# Check if any hoppings were added
if hopping_count == 0:
    print("\n" + "="*30)
    print(" WARNING: No non-zero hopping terms were added.")
    print(f" Check cutoff distance ({cutoff_distance} Å) or SK parameters/formulas.")
    print("="*30 + "\n")
else:
    # Calculate band structure
    k_points = np.array([
        [0.0, 0.0, 0.0],         # Γ
        [0.5, 0.0, 0.0],         # M
        [1/3, 1/3, 0.0],         # K
        [0.0, 0.0, 0.0]          # Γ
    ])
    k_label = ["Γ", "M", "K", "Γ"]
    k_path, k_idx = tb.gen_kpath(k_points, [50, 50, 50])
    print("DEBUG: Calculating band structure...")
    try:
        k_len, bands = cell.calc_bands(k_path)
        print("DEBUG: Calculation complete.")

        # Plot band structure
        print("DEBUG: Plotting...")
        plt.figure(figsize=(8, 6))
        vis = tb.Visualizer()
        vis.plot_bands(k_len, bands, k_idx=k_idx, k_label=k_label, fig_name="newest.png")

        plt.title(f"MoS2 Band Structure (Corrected Params, Fixed Custom SK1954, cutoff={cutoff_distance}Å)")
        plt.ylim(-8.5, 4.5) # Focus near the gap
        plt.ylabel("Energy (eV)")
        plt.axhline(0, color='grey', linestyle='--', linewidth=0.5)
        plt.grid(True, axis='y', linestyle=':', linewidth=0.5)
        output_filename = f"mos2_bands_corrected_params_custom_sk_v3.png" # New filename
        plt.savefig(output_filename, dpi=300, bbox_inches='tight')
        print(f"DEBUG: Band structure plot saved to {output_filename}")
        # plt.show()
        plt.close()

    except Exception as e:
        print(f"ERROR during band structure calculation or plotting: {e}")

print("\nDEBUG: Script finished.")
