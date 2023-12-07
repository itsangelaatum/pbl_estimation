import matplotlib.pyplot as plt
import math

def plot_colormap(longitude, latitude, density, title='', f_size=(10,5)):
    plt.figure(figsize=f_size)
    if title != '':
        plt.title(title, fontsize=20)
    plt.ylabel("Latitude", fontsize=15)
    plt.xlabel("Longitude", fontsize=15)

    lons, lats = np.meshgrid(longitude, latitude)

    plt.contourf(lons, lats, density, cmap='viridis')
    plt.colorbar()

    min_lon = math.floor(min(longitude))
    max_lon = math.ceil(max(longitude))+1
    min_lat = math.floor(min(latitude))
    max_lat = math.ceil(max(latitude))+1
    plt.xticks(range(min_lon, max_lon, 30))
    plt.yticks(range(min_lat, max_lat, 30))

    plt.xlim(min_lon, max_lon)
    plt.ylim(min_lat, max_lat)

    plt.show()
