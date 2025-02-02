import cv2

# Read the image in grayscale
image = cv2.imread(r"C:\Users\woods\OneDrive\Documents\Research\Full Branching\Naga\Naga1 CM.png", cv2.IMREAD_GRAYSCALE)

# Check if image is loaded
if image is None:
    print("Error: Could not read the image.")
    exit()

# Apply binary thresholding
_, binary_image = cv2.threshold(image, 128, 255, cv2.THRESH_BINARY)

# Save and display the binary image
cv2.imwrite("binary_image.jpg", binary_image)
cv2.imshow("Binarized Image", binary_image)
cv2.waitKey(0)
cv2.destroyAllWindows()
