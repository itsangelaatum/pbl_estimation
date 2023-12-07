import matplotlib.pyplot as plt
import numpy as np

# Sample latitude and longitude data (replace with your own data)
latitudes = np.random.uniform(-90, 90, 100)
longitudes = np.random.uniform(-180, 180, 100)

# Create a 2D histogram
heatmap, xedges, yedges = np.histogram2d(latitudes, longitudes, bins=100)

# Create a figure and axis
fig, ax = plt.subplots()

# Plot the heatmap
cax = ax.imshow(heatmap.T, cmap='YlOrRd', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])

# Add a colorbar
cbar = fig.colorbar(cax, label='Frequency')

# Set labels and title
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title('Latitude-Longitude Heatmap')

# Show the plot
plt.show()
