import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve
from typing import Tuple
import cv2

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


def visualize_elliptical_partitions(mylog: np.ndarray, sz: Tuple[int, int], maxboxes: int = 5):
    """
    Visualizes the elliptical partitions used in the box-counting method.

    Parameters:
        mylog (np.ndarray): Binary mask indicating object presence.
        sz (Tuple[int, int]): Size of the image (height, width).
        maxboxes (int): Maximum number of subdivisions (default: 5).
    """
    # Generate coordinate vectors
    xvec = np.linspace(-sz[1] // 2, sz[1] // 2, sz[1])
    yvec = np.linspace(-sz[0] // 2, sz[0] // 2, sz[0])

    # Create mesh grids for spatial coordinates
    Xim, Yim = np.meshgrid(xvec, yvec)

    # Convert Cartesian coordinates to polar coordinates
    theta, rho = np.arctan2(Yim, Xim), np.hypot(Xim, Yim)

    # Calculate total area and aspect ratio
    totalarea = np.pi * max(xvec) * max(yvec)
    aspect = max(max(xvec), max(yvec)) / min(max(xvec), max(yvec))

    # Plot the original binary mask
    plt.figure(figsize=(8, 8))
    plt.imshow(mylog, cmap='gray', origin='lower', extent=[xvec[0], xvec[-1], yvec[0], yvec[-1]])

    # Box-counting loop for drawing partitions
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

        # **Draw angular partitions (radial lines)**
        for theta_val in thetarange:
            x_line = [majorrange[-1] * np.cos(theta_val), -majorrange[-1] * np.cos(theta_val)]
            y_line = [majorrange[-1] * np.sin(theta_val), -majorrange[-1] * np.sin(theta_val)]
            plt.plot(x_line, y_line, 'r-', alpha=0.5, linewidth=1)

        # **Draw radial partitions (elliptical rings)**
        theta_grid = np.linspace(0, 2 * np.pi, 200)
        for j in range(1, len(majorrange)):  # Skip 0 to avoid a point at the center
            x_ellipse = majorrange[j] * np.cos(theta_grid)
            y_ellipse = minorrange[j] * np.sin(theta_grid)
            plt.plot(x_ellipse, y_ellipse, 'b-', alpha=0.5, linewidth=1)

    # Labels and display
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.title("Elliptical Partitions for Box Counting")
    plt.axis("equal")  # Keep aspect ratio consistent
    plt.show()

image_path = r"file.png"
threshold = 128
   
binary_image = load_and_preprocess_image(image_path, threshold)
plt.imshow(binary_image)

# Create and apply elliptical mask
ellipse_mask = create_ellipse_mask(binary_image.shape)
mylog = binary_image * ellipse_mask

sz = mylog.shape

visualize_elliptical_partitions(mylog, sz, maxboxes=5)