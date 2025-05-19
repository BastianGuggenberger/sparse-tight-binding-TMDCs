import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh, svd
import argparse

class BandStructure:
    def __init__(self, stc_file=None):
        self.norb = 22  # need to fix this to read from beloew line where atoms is found [LINE 6 in stencil format]
        if stc_file:
            self.hoppings = self.parse_stc(stc_file)
        else:
            self.hoppings = {}
        
    def parse_stc(self, filename):
        """Parse .stc file into hopping dictionary"""
        hoppings = {}
        current_R = None
        
        with open(filename, 'r') as f:
            lines = f.readlines()
            
        for line in lines:
            if line.strip():
                parts = line.split()
                # Check if it's an R vector line
                if len(parts) == 3 and all(p.replace('-','').isdigit() for p in parts):
                    current_R = (int(parts[0]), int(parts[1]))
                    hoppings[current_R] = np.zeros((self.norb, self.norb), dtype=complex)
                # Check if it's a hopping line
                elif current_R is not None and '(' in line:
                    parts = line.split('(')
                    indices = parts[0].split()
                    value = float(parts[1].split(',')[0])
                    i, j = int(indices[0]), int(indices[1])
                    hoppings[current_R][i,j] = value
                    
        return hoppings
    
    def construct_hk(self, k):
        """Construct Hamiltonian at k-point"""
        Hk = np.zeros((self.norb, self.norb), dtype=complex)
        
        for R, t_matrix in self.hoppings.items():
            phase = np.exp(2j * np.pi * (k[0]*R[0] + k[1]*R[1]))
            Hk += t_matrix * phase
            
        # Make Hermitian
        Hk = 0.5 * (Hk + Hk.conj().T)
        
        return Hk
    
    def get_bands(self, k_path):
        """Calculate bands along k-path"""
        bands = []
        for k in k_path:
            Hk = self.construct_hk(k)
            eigenvals = eigh(Hk, eigvals_only=True)
            bands.append(eigenvals)
            
        return np.array(bands)

def normalize_hoppings(hoppings):
    """Normalize hopping values for consistent comparison"""
    # Convert to real if all values are essentially real
    for R, mat in hoppings.items():
        if np.all(np.abs(np.imag(mat)) < 1e-10):
            hoppings[R] = np.real(mat)
    return hoppings

def create_norm_based_reduction(bs_orig, n_components):
    """Implement norm-based reduction method"""
    # Calculate norms and find important vectors
    vector_norms = {R: np.linalg.norm(hopping) 
                   for R, hopping in bs_orig.hoppings.items()}
    
    important_Rs = sorted(vector_norms.keys(), 
                         key=lambda R: vector_norms[R], 
                         reverse=True)[:n_components]
    
    # Create reduced hopping set
    hoppings_reduced = {}
    for R in important_Rs:
        hoppings_reduced[R] = bs_orig.hoppings[R].copy()
    
    # Create reduced Band Structure object
    bs_red = BandStructure()
    bs_red.norb = bs_orig.norb  # copy number of orbitals
    bs_red.hoppings = normalize_hoppings(hoppings_reduced)
    
    return bs_red, vector_norms

def create_svd_based_reduction(bs_orig, n_components):
    """Implement SVD-based reduction method"""
    # Convert to matrix for SVD
    R_vectors = list(bs_orig.hoppings.keys())
    hoppings_list = list(bs_orig.hoppings.values())
    
    # Create the matrix T (R_vectors × (norb*norb))
    T_matrix = np.array([h.flatten() for h in hoppings_list])
    
    # Perform SVD
    U, s, Vh = svd(T_matrix, full_matrices=False)
    
    # Truncate to keep only top N singular values
    U_N = U[:, :n_components]
    s_N = s[:n_components]
    Vh_N = Vh[:n_components, :]
    
    # Reconstruct the reduced matrix
    T_reduced = U_N @ np.diag(s_N) @ Vh_N
    
    # Convert back to hopping dictionary format
    hoppings_svd = {}
    for i, R in enumerate(R_vectors):
        hoppings_svd[R] = T_reduced[i].reshape(bs_orig.norb, bs_orig.norb)
    
    # Enforce Hermiticity
    for R, t_matrix in list(hoppings_svd.items()):
        R_minus = (-R[0], -R[1])
        if R_minus not in hoppings_svd:
            hoppings_svd[R_minus] = t_matrix.T.conj()
        else:
            avg = 0.5 * (t_matrix + hoppings_svd[R_minus].T.conj())
            hoppings_svd[R] = avg
            hoppings_svd[R_minus] = avg.T.conj()
    
    # Create band structure object
    bs_svd = BandStructure()
    bs_svd.norb = bs_orig.norb  # copy number of orbitals
    bs_svd.hoppings = normalize_hoppings(hoppings_svd)
    
    return bs_svd, s

def compare_methods(stc_file, n_components=7, output_file='comprehensive_comparison.png'):
    """Compare all three methods in a single plot with difference panels"""
    # Load original bands
    bs_orig = BandStructure(stc_file)
    
    # Create norm-based reduction
    bs_norm, vector_norms = create_norm_based_reduction(bs_orig, n_components)
    
    # Create SVD-based reduction
    bs_svd, singular_values = create_svd_based_reduction(bs_orig, n_components)
    
    # Set up k-path
    k_points = np.array([
        [0, 0],      # Γ
        [1/3, 1/3],  # K
        [0.5, 0],    # M
        [0, 0]       # Γ
    ])
    k_labels = ['Γ', 'K', 'M', 'Γ']
    
    # Generate k-path
    nk = 50
    k_path = []
    k_indices = []
    
    for i in range(len(k_points)-1):
        k_start = k_points[i]
        k_end = k_points[i+1]
        for j in range(nk):
            t = j/nk
            k = k_start + t*(k_end - k_start)
            k_path.append(k)
        k_indices.append(i*nk)
    k_path.append(k_points[-1])
    k_indices.append(len(k_path)-1)
    k_path = np.array(k_path)
    
    # Calculate bands
    bands_orig = bs_orig.get_bands(k_path)
    bands_norm = bs_norm.get_bands(k_path)
    bands_svd = bs_svd.get_bands(k_path)
    
    # Calculate errors
    norm_error = np.mean(np.abs(bands_orig - bands_norm))
    norm_max_error = np.max(np.abs(bands_orig - bands_norm))
    svd_error = np.mean(np.abs(bands_orig - bands_svd))
    svd_max_error = np.max(np.abs(bands_orig - bands_svd))
    
    # Create figure with 2 rows, 2 columns
    fig = plt.figure(figsize=(12, 10))
    
    # Main band structure plot (top)
    ax1 = plt.subplot2grid((2, 2), (0, 0), colspan=2)
    x = np.arange(len(k_path))
    
    # Plot all three band structures in one graph
    for i, band in enumerate(bands_orig.T):
        if i == 0:
            ax1.plot(x, band, 'k-', linewidth=2, alpha=0.8, label='Original')
        else:
            ax1.plot(x, band, 'k-', linewidth=2, alpha=0.8)
            
    for i, band in enumerate(bands_norm.T):
        if i == 0:
            ax1.plot(x, band, 'b--', linewidth=1.5, alpha=0.7, label='Norm-Based')
        else:
            ax1.plot(x, band, 'b--', linewidth=1.5, alpha=0.7)
            
    for i, band in enumerate(bands_svd.T):
        if i == 0:
            ax1.plot(x, band, 'r:', linewidth=1.5, alpha=0.7, label='SVD-Based')
        else:
            ax1.plot(x, band, 'r:', linewidth=1.5, alpha=0.7)
    
    ax1.set_title(f'Band Structure Comparison - NbSe₂\nOriginal ({len(bs_orig.hoppings)} R vectors) vs Reduced (N={n_components})', fontsize=14)
    ax1.set_xlabel('k-path')
    ax1.set_ylabel('Energy (eV)')
    ax1.set_xticks(k_indices)
    ax1.set_xticklabels(k_labels)
    ax1.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper right')
    
    # Error plot (bottom left)
    ax2 = plt.subplot2grid((2, 2), (1, 0))
    
    # Calculate band-by-band error
    norm_band_errors = np.mean(np.abs(bands_orig - bands_norm), axis=0)
    svd_band_errors = np.mean(np.abs(bands_orig - bands_svd), axis=0)
    
    band_indices = np.arange(len(norm_band_errors))
    width = 0.35
    
    ax2.bar(band_indices - width/2, norm_band_errors, width, label='Norm-Based', color='blue', alpha=0.7)
    ax2.bar(band_indices + width/2, svd_band_errors, width, label='SVD-Based', color='red', alpha=0.7)
    
    ax2.set_title('Mean Error by Band')
    ax2.set_xlabel('Band Index')
    ax2.set_ylabel('Mean Error (eV)')
    ax2.set_xticks(band_indices)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Add text box with overall error metrics
    error_text = (
        f"Norm-Based Reduction:\n"
        f"Mean Error: {norm_error:.4f} eV\n"
        f"Max Error: {norm_max_error:.4f} eV\n\n"
        f"SVD-Based Reduction:\n"
        f"Mean Error: {svd_error:.4f} eV\n"
        f"Max Error: {svd_max_error:.4f} eV"
    )
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax2.text(1.05, 0.95, error_text, transform=ax2.transAxes, fontsize=10,
            verticalalignment='top', bbox=props)
    
    # Singular value / Norm plot (bottom right)
    ax3 = plt.subplot2grid((2, 2), (1, 1))
    
    # Plot singular values
    sorted_norms = sorted(vector_norms.items(), key=lambda x: x[1], reverse=True)
    top_norms = [norm for _, norm in sorted_norms[:n_components*2]]
    
    ax3.semilogy(np.arange(len(singular_values[:n_components*2])), singular_values[:n_components*2], 'ro-', label='Singular Values')
    ax3.semilogy(np.arange(len(top_norms)), top_norms, 'bs-', label='Vector Norms')
    
    ax3.axvline(x=n_components-0.5, color='k', linestyle='--', alpha=0.5, label=f'N={n_components} Cutoff')
    ax3.set_title('Singular Values vs Vector Norms')
    ax3.set_xlabel('Index')
    ax3.set_ylabel('Magnitude (log scale)')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.show()
    
    # Print analysis
    print(f"\nComparison with N={n_components} components:")
    print("\nNorm-Based Reduction:")
    print(f"Mean Error: {norm_error:.6f} eV")
    print(f"Max Error: {norm_max_error:.6f} eV")
    
    print("\nSVD-Based Reduction:")
    print(f"Mean Error: {svd_error:.6f} eV")
    print(f"Max Error: {svd_max_error:.6f} eV")
    
    # Print R-vector details
    print(f"\nOriginal R-vectors: {len(bs_orig.hoppings)}")
    print(f"Reduced R-vectors: {n_components} ({n_components/len(bs_orig.hoppings)*100:.1f}%)")
    
    return {
        'norm_error': norm_error,
        'norm_max_error': norm_max_error,
        'svd_error': svd_error,
        'svd_max_error': svd_max_error,
        'original_vectors': len(bs_orig.hoppings),
        'reduced_vectors': n_components
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare band structure reduction methods.')
    parser.add_argument('--stc_file', type=str, required=True, help='Path to the .stc file')
    parser.add_argument('--n_components', type=int, default=20, help='Number of components for reduction')
    parser.add_argument('--output', type=str, default='comprehensive_comparison.png', help='Path to output plot file')
    args = parser.parse_args()
    
    compare_methods(args.stc_file, n_components=args.n_components, output_file=args.output)
