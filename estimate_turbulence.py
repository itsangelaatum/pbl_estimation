import numpy as np
from gmi_daily_v8 import GMIdaily
from example_usage import *

def read_data(filename='f35_20230918v8.2.gz'):
    dataset = GMIdaily(filename, missing=np.nan)
    if not dataset.variables: sys.exit('problem reading file')
    return dataset

def estimate_turbulence_struct(diff):
    # D_IWV_hat(s) = <[IWV(x)-IWV(x+s)]^2>
    diff = np.array(diff)
    d_iwv_hat = np.mean(np.square(diff))

    # TODO: get and subtract noise bias
    # noise_bias = 0
    # d_iwv = d_iwv_hat - noise_bias

    return d_iwv_hat
    
def estimate_beta(s_diwv_dict):
    # log(D_IWV(s)) = beta log(s) + log(C_IWV)
    # linear_regression
    s_s = [s for s in s_diwv_dict.keys()]
    diwv_s = []
    [diwv_s.extend(s_diwv_dict[s]) for s in s_s]

    y_s = np.log(diwv_s)

    x_s = np.log(s_s*np.ones(np.size(diwv_s)))
    A = np.vstack([x_s, np.ones(len(x_s))]).T
    beta, log_ciwv = np.linalg.lstsq(A, y_s, rcond=None)[0]
    return (beta, log_ciwv)


def get_adjacent_differences(dataset, lat, lon, k, iasc = 0, lat_max = 720, lon_max = 1440):
    iwv_x = dataset[iasc, lat, lon]
    differences = []
    if not np.isnan(iwv_x):
        if lat - k >= 0:
            iwv_xs = dataset[iasc, lat - k, lon]
            differences.append(iwv_x-iwv_xs) if not np.isnan(iwv_xs) else 0
        if lat + k < lat_max:
            iwv_xs = dataset[iasc, lat + k, lon]
            differences.append(iwv_x-iwv_xs) if not np.isnan(iwv_xs) else 0
        if lon - k >= 0:
            iwv_xs = dataset[iasc, lat, lon - k]
            differences.append(iwv_x-iwv_xs) if not np.isnan(iwv_xs) else 0
        if lon + k < lon_max:
            iwv_xs = dataset[iasc, lat, lon + k]
            differences.append(iwv_x-iwv_xs) if not np.isnan(iwv_xs) else 0
    return differences

def get_diag_differences(dataset, lat, lon, k, iasc = 0, lat_max = 720, lon_max = 1440):
    iwv_x = dataset[iasc, lat, lon]
    differences = []
    if not np.isnan(iwv_x):
        if lat - k >= 0 and lon - k >= 0:
            iwv_xs = dataset[iasc, lat - k, lon]
            differences.append(iwv_x-iwv_xs) if not np.isnan(iwv_xs) else 0
        if lat + k < lat_max and lon - k >= 0:
            iwv_xs = dataset[iasc, lat + k, lon]
            differences.append(iwv_x-iwv_xs) if not np.isnan(iwv_xs) else 0
        if lat - k >= 0 and lon + k < lon_max:
            iwv_xs = dataset[iasc, lat - k, lon - k]
            differences.append(iwv_x-iwv_xs) if not np.isnan(iwv_xs) else 0
        if lat + k < lat_max and lon + k < lon_max:
            iwv_xs = dataset[iasc, lat + k, lon + k]
            differences.append(iwv_x-iwv_xs) if not np.isnan(iwv_xs) else 0
    return differences

def estimate_beta_for_area(dataset, lat_min, lon_min, lat_max, lon_max): 
    # print(f"latitude range = [{lat_min},{lat_max}]")
    # print(f"longitude range = [{lon_min},{lon_max}]")

    s_diwv_dict = {}
    for k in range(1,4):
        # 0.25 deg between each measurement and each 0.25deg = 27.75km
        s_adjacent, s_diag = round(k * 0.25 * 27.75, 4), round(k * 0.25 * np.sqrt(2) * 27.75, 4)
        adjacent_differences, diag_differences = [], []
        s_diwv_dict[s_adjacent] , s_diwv_dict[s_diag] = [], []
        for lat in range(lat_min, lat_max):
            for lon in range(lon_min, lon_max):
                adjacent_differences.extend(get_adjacent_differences(dataset, lat, lon, k))
                diag_differences.extend(get_diag_differences(dataset, lat, lon, k))

        d_iwv_hat = estimate_turbulence_struct(adjacent_differences)
        s_diwv_dict[s_adjacent].append(d_iwv_hat)
        # print(f"s = {round(s_adjacent, 4)}, and d_iwv(s) = {d_iwv_hat}")

        d_iwv_hat = estimate_turbulence_struct(diag_differences)
        s_diwv_dict[s_diag].append(d_iwv_hat)
        # print(f"s = {round(s_diag, 4)}, and d_iwv(s) = {d_iwv_hat}")
    
    beta, c_iwv = estimate_beta(s_diwv_dict)
    if not np.isnan(beta):
        print(f"latitude range = [{lat_min},{lat_max}]")
        print(f"longitude range = [{lon_min},{lon_max}]")
        print(f"beta = {round(beta, 4)}, and C_iwv(s) = {c_iwv}")
    return beta

# def plot_lat_lon_beta(latlon_beta_dict):


def run_algorithm(dataset, area = 20):
    latlon_beta_dict = {}
    for i in range(0,720//area-1):
        for j in range(0,1440//area-1):
            beta = estimate_beta_for_area(dataset, i*area, (i+1)*area, j*area, (j+1)*area)
            latlon_beta_dict[(i*area, j*area)] = beta
    print(latlon_beta_dict)
    


if __name__ == '__main__':
    import sys
    dataset = read_data()
    iwv = dataset.variables['vapor']
    show_dimensions(dataset)
    # show_variables(dataset)
    # show_validrange(dataset)
    # show_somedata(dataset)
    run_algorithm(iwv)
