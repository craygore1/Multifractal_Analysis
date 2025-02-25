import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve
from scipy.stats import linregress
import cv2
from typing import Tuple, List
import time

start_time = time.perf_counter()

def load_and_preprocess_image(image_path: str, threshold: int = 128) -> np.ndarray:
    """
    Loads an image in grayscale and applies binary thresholding.
    """
    image = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)

    if image is None:
        raise FileNotFoundError(f"Error: Could not read the image at {image_path}")
    
    # Apply binary thresholding
    _, binary_image = cv2.threshold(image, threshold, 255, cv2.THRESH_BINARY)
    
    # Convert pixel values from 255 to 1 for easier processing
    binary_image[binary_image == 255] = 1
    return (binary_image == 1).astype(int)

def create_ellipse_mask(image_shape: Tuple[int, int]) -> np.ndarray:
    """
    Creates an elliptical mask for the given image dimensions.
    """
    height, width = image_shape
    center_x, center_y = width // 2, height // 2
    radius_x, radius_y = center_x, center_y

    y, x = np.meshgrid(np.arange(height), np.arange(width), indexing="ij")
    return ((y - center_y) ** 2 / radius_y**2 + (x - center_x) ** 2 / radius_x**2) <= 1

def compute_box_counting(mylog: np.ndarray, sz: Tuple[int, int], maxboxes: int = 5) -> Tuple[np.ndarray, np.ndarray]:
    """
    Computes elliptical box-counting for multifractal analysis.

    Parameters:
        mylog (np.ndarray): Binary mask indicating object presence.
        sz (Tuple[int, int]): Size of the image (height, width).
        maxboxes (int): Maximum number of subdivisions (default: 5).

    Returns:
        Tuple[np.ndarray, np.ndarray]: Normalized probabilities (prbM), log-scaled box sizes (X).
    """
    # Generate coordinate vectors
    xvec = np.round(np.linspace(-sz[1] // 2, sz[1] // 2, sz[1]))
    yvec = np.round(np.linspace(-sz[0] // 2, sz[0] // 2, sz[0]))

    # Create mesh grids for spatial coordinates
    Xim, Yim = np.meshgrid(xvec, yvec)

    # Convert Cartesian coordinates to polar coordinates
    theta, rho = np.arctan2(Yim, Xim), np.hypot(Xim, Yim)

    # Extract valid indices where the object exists
    myind = np.where(mylog)
    thetaind, rhoind = theta[myind], rho[myind]
    Xind, Yind = Xim[myind], Yim[myind]

    # Calculate total area and aspect ratio
    totalarea = np.pi * max(xvec) * max(yvec)
    aspect = max(max(xvec), max(yvec)) / min(max(xvec), max(yvec))

    # Initialize storage matrix
    M = np.zeros((2**(2*maxboxes), maxboxes + 1))  # Resizing for proper 2D bin counts

    # Box-counting loop
    for i in range(maxboxes + 1):
        # Construct linear system for subdivision
        d1 = np.pi * np.ones(2**i)
        d2 = -np.pi * np.ones(2**i - 1)
        A = np.diag(d1) + np.diag(d2, -1)

        # Compute areas and solve for major/minor axis lengths
        areas = (totalarea / (2**i)) * np.ones(2**i)
        ab = solve(A, areas)
        minorrange = np.sqrt(ab / aspect)
        majorrange = ab / minorrange

        # Insert zero at the beginning to avoid indexing errors
        minorrange = np.insert(minorrange, 0, 0)
        majorrange = np.insert(majorrange, 0, 0)

        # Avoid division by zero
        majorrange = np.clip(np.nan_to_num(majorrange, nan=1e-10), 1e-10, None)
        minorrange = np.clip(np.nan_to_num(minorrange, nan=1e-10), 1e-10, None)

        # Define angular bins
        thetarange = np.linspace(-np.pi, np.pi, (2**i) + 1)

        # **Vectorized Assignment of Bins**
        theta_bins = np.digitize(thetaind, thetarange) - 1  # Assign angular bins

        # Compute normalized squared distances for all points at once
        norm_sq = (Xind[:, None]**2 / majorrange[None, :-1]**2) + (Yind[:, None]**2 / minorrange[None, :-1]**2)

        # Assign each point to the correct elliptical ring (vectorized)
        radial_bin_indices = np.sum(norm_sq >= 1, axis=1) - 1  # Find the highest valid bin

        # **Ensure all indices are within bounds**
        theta_bins = np.clip(theta_bins, 0, len(thetarange) - 2)
        radial_bin_indices = np.clip(radial_bin_indices, 0, len(majorrange) - 2)

        # Use a 2D histogram to properly partition across **both** dimensions
        histogram, _, _ = np.histogram2d(theta_bins, radial_bin_indices, 
                                         bins=[len(thetarange) - 1, len(majorrange) - 1])

        # Flatten histogram and store
        M[:histogram.size, i] = histogram.flatten()

    # Normalize probabilities
    prbM = M / (np.sum(M, axis=0) + 1e-10)

    # Compute log-scaled box sizes
    truesz = (2 * np.pi) / (2 ** np.arange(maxboxes + 1))
    X = np.log2(truesz)

    return prbM, X

def compute_multifractal_spectrum(prbM: np.ndarray, qvals: np.ndarray, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Computes the multifractal spectrum (Dq, tauq, myalpha, falpha) from probability measures.

    Parameters:
        prbM (np.ndarray): Probability matrix of shape (n_regions, n_scales).
        qvals (np.ndarray): Array of q-values for multifractal analysis.
        X (np.ndarray): Log-scaled box sizes.

    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        - Dq (np.ndarray): Generalized fractal dimensions.
        - tauq (np.ndarray): Mass exponent function.
        - myalpha (np.ndarray): Singularity strengths.
        - falpha (np.ndarray): Multifractal spectrum.
    """
    num_qvals, num_scales = len(qvals), len(prbM[0])

    # Initialize output arrays
    yD, yalph, yf = np.zeros((num_qvals, num_scales)), np.zeros((num_qvals, num_scales)), np.zeros((num_qvals, num_scales))
    Dq, tauq, myalpha, falpha = np.zeros(num_qvals), np.zeros(num_qvals), np.zeros(num_qvals), np.zeros(num_qvals)

    # Compute multifractal spectrum
    for idx, k in enumerate(qvals):
        for a in range(num_scales):
            nonzero_prbM = prbM[:, a][prbM[:, a] > 0]  # Exclude zero probabilities

            if k == 1:
                yD[idx, a] = np.sum(nonzero_prbM * np.log2(nonzero_prbM))
            else:
                sum_powered = np.maximum(np.sum(nonzero_prbM**k), 1e-30)
                yD[idx, a] = np.log2(sum_powered)

                mu = (nonzero_prbM**k) / sum_powered
                yalph[idx, a] = np.sum(mu * np.log2(nonzero_prbM))
                yf[idx, a] = np.sum(mu * np.log2(mu))

    # Compute fractal dimensions and multifractal spectrum
    for idx_q, q in enumerate(qvals):
        slope, _, _, _, _ = linregress(X, yD[idx_q])

        if q == 1:
            Dq[idx_q] = abs(slope)
        else:
            tauq[idx_q] = slope
            Dq[idx_q] = tauq[idx_q] / (q - 1)
            myalpha[idx_q] = linregress(X, yalph[idx_q])[0]
            falpha[idx_q] = linregress(X, yf[idx_q])[0]

    return Dq, tauq, myalpha, falpha


def plot_multifractal_spectrum(qvals: np.ndarray, Dq: np.ndarray, myalpha: np.ndarray, falpha: np.ndarray, plots: bool = True) -> None:
    """
    Plots the multifractal spectrum, including D(q) vs. q and f(α) vs. α.

    Parameters:
        qvals (np.ndarray): Array of q-values.
        Dq (np.ndarray): Generalized fractal dimensions D(q).
        myalpha (np.ndarray): Singularity strengths α.
        falpha (np.ndarray): Multifractal spectrum f(α).
        plots (bool): Whether to generate plots (default: True).
    """
    if not plots:
        return  # Skip plotting if disabled

    # Plot D(q) vs. q
    plt.figure()
    plt.plot(qvals, Dq, color='#0F6FC6', linewidth=1.25)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.xlabel('$q$', fontsize=16)
    plt.ylabel('$D(q)$', fontsize=16)
    plt.show()

    # Plot f(α) vs. α
    plt.figure()
    plt.plot(myalpha, falpha, color='#0F6FC6', marker='.')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.xlabel(r'$\alpha$', fontsize=16)
    plt.ylabel(r'$f(\alpha)$', fontsize=16)
    plt.xlim(0, 5)
    plt.ylim(0, 2)
    plt.show()

def main(image_path: str, qvals: np.ndarray, maxboxes: int = 5, threshold: int = 128, plots: bool = True) -> None:
    """
    Function to perform multifractal analysis on an image.

    Parameters:
        image_path (str): Path to the input image.
        qvals (np.ndarray): Array of q-values for multifractal analysis.
        maxboxes (int): Maximum number of subdivisions (default: 5).
        threshold (int): Threshold for binary image processing (default: 128).
        plots (bool): Whether to generate plots (default: True).
    """
    # Load and preprocess the image
    binary_image = load_and_preprocess_image(image_path, threshold)
    plt.imshow(binary_image)

    # Create and apply elliptical mask
    ellipse_mask = create_ellipse_mask(binary_image.shape)
    mylog = binary_image * ellipse_mask
    #plt.imshow(mylog)

    # Compute the box-counting probabilities and log-scaled box sizes
    prbM, X = compute_box_counting(mylog, binary_image.shape, maxboxes)

    # Compute multifractal spectrum
    Dq, tauq, myalpha, falpha = compute_multifractal_spectrum(prbM, qvals, X)

    # Plot results
    plot_multifractal_spectrum(qvals, Dq, myalpha, falpha, plots)
    
    return prbM

if __name__ == "__main__":
    image_path = r"C:\Users\woods\OneDrive\Documents\Research\Full Branching\Naga\Naga1 CM.png"
    qvals = np.linspace(-10,10, 100)  # Initialize q-values
    prbM_test = main(image_path, qvals, maxboxes=6, threshold=128, plots=True)


end_time = time.perf_counter()
execution_time = end_time - start_time
print(f"Execution time: {execution_time} seconds")

