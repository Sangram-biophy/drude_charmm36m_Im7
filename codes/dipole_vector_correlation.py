import numpy as np
from tqdm import tqdm


def compute_neighbor_dipole_cosines(dipoles):
    """
    Computes average cos(θ) between neighboring dipole vectors (i and i+1) per frame.
    
    Parameters:
        dipoles: np.ndarray of shape (n_frames, n_residues, 3)
        
    Returns:
        avg_cosines: np.ndarray of shape (n_residues - 1,)
                     Average cos(θ) between dipoles i and i+1 over time
    """
    n_frames, n_residues, _ = dipoles.shape

    # Normalize dipole vectors
    norm_dipoles = dipoles / np.linalg.norm(dipoles, axis=2, keepdims=True)

    # Compute cosine between adjacent residues for each frame
    cosines = np.sum(
        norm_dipoles[:, :-1, :] * norm_dipoles[:, 1:, :], axis=2
    )  # shape: (n_frames, n_residues - 1)

    # Average over time
    avg_cosines = np.mean(cosines, axis=0)  # shape: (n_residues - 1,)
    sd_cosine = np.std(cosines, axis=0)
    return avg_cosines,sd_cosine

def compute_dipole_correlation_stats(dipoles):
    """
    Computes the mean and standard deviation of dipole-dipole cosines
    between all residue pairs over time.
    
    Parameters:
        dipoles: np.ndarray of shape (n_frames, n_residues, 3)
    
    Returns:
        mean_corr: np.ndarray of shape (n_residues, n_residues)
                   Mean cosine of angle between dipole vectors of residues i and j
        std_corr: np.ndarray of shape (n_residues, n_residues)
                  Standard deviation of cosine over time
    """
    n_frames, n_residues, _ = dipoles.shape

    # Normalize dipoles to unit vectors
    norm_dipoles = dipoles / np.linalg.norm(dipoles, axis=2, keepdims=True)

    # Preallocate array to store time series of correlations
    corr_time_series = np.empty((n_frames, n_residues, n_residues))

    for f in range(n_frames):
        v = norm_dipoles[f]  # shape: (n_residues, 3)
        corr_time_series[f] = np.dot(v, v.T)  # cosine similarity matrix at frame f

    # Mean and std over time
    mean_corr = np.mean(corr_time_series, axis=0)
    std_corr = np.std(corr_time_series, axis=0)

    return mean_corr, std_corr

dipole_vector=np.load("residuewise_backbone_dipole_vector_drude_LP_skip2.npy")
avg_cosines,sd_cosine=compute_neighbor_dipole_cosines(dipole_vector)
resid=np.arange(1,np.shape(dipole_vector)[1])
np.save("cosine_dipole_vector_neighbour.npy",np.array([resid,avg_cosines,sd_cosine]).T,allow_pickle=True)
mean_corr, std_corr=compute_dipole_correlation_stats(dipole_vector)
np.save("mean_correlation_dipole_vector.npy",mean_corr,allow_pickle=True)
np.save("std_correlation_dipole_vector.npy",std_corr,allow_pickle=True)
