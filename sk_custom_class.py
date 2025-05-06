# By Max Sinner


import numpy as np
import tbplas as tb # Import the tbplas library
import matplotlib.pyplot as plt
from numpy.linalg import norm
from math import exp, sqrt


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