from itertools import combinations, product

class TurbulenceStructureEstimator():
    def __init__(self, dataset, area, ilat_min=0, ilon_min=0, ilat_max=720, ilon_max=1440, verbose=False):
        self.dataset = dataset
        self.area = area
        self.ilat_min = ilat_min
        self.ilon_min = ilon_min
        self.ilat_max = ilat_max
        self.ilon_max = ilon_max
        self.area_sdiwv_dict = {}
        self.area_sdiffs_dict = {}
        self.verbose = verbose

    def estimate_turbulence_struct(self, diff):
        # D_IWV_hat(s) = <[IWV(x)-IWV(x+s)]^2>
        diff = np.array(diff)
        d_iwv_hat = np.mean(np.square(diff))

        # TODO: get and subtract noise bias
        # noise_bias = 0
        # d_iwv = d_iwv_hat - noise_bias

        return d_iwv_hat

    def estimate_beta(self, s_diwv_dict):
        # log(D_IWV(s)) = beta log(s) + log(C_IWV)
        # linear_regression
        s_s = list(s_diwv_dict.keys())
        diwv_s = list(s_diwv_dict.values())

        y_s = np.log(diwv_s)

        x_s = np.log(s_s*np.ones(np.size(diwv_s)))
        A = np.vstack([x_s, np.ones(len(x_s))]).T
        beta, log_ciwv = np.linalg.lstsq(A, y_s, rcond=None)[0]
        ciwv = np.exp(log_ciwv)
        return (beta, ciwv)

    def get_adjacent_differences(self, lat, lon, k):
        iwv_x = self.dataset[lat, lon]
        differences = []
        if not np.isnan(iwv_x):
            if lat - k >= self.ilat_min:
                iwv_xs = self.dataset[lat - k, lon]
                differences.append(iwv_x-iwv_xs) if not np.isnan(iwv_xs) else 0
            if lat + k < self.ilat_max:
                iwv_xs = self.dataset[lat + k, lon]
                differences.append(iwv_x-iwv_xs) if not np.isnan(iwv_xs) else 0
            if lon - k >= self.ilon_min:
                iwv_xs = self.dataset[lat, lon - k]
                differences.append(iwv_x-iwv_xs) if not np.isnan(iwv_xs) else 0
            if lon + k < self.ilon_max:
                iwv_xs = self.dataset[lat, lon + k]
                differences.append(iwv_x-iwv_xs) if not np.isnan(iwv_xs) else 0
        return differences

    def get_diag_differences(self, lat, lon, k):
        iwv_x = self.dataset[lat, lon]
        differences = []
        if not np.isnan(iwv_x):
            if lat - k >= self.ilat_min and lon - k >= self.ilon_min:
                iwv_xs = self.dataset[lat - k, lon]
                differences.append(iwv_x-iwv_xs) if not np.isnan(iwv_xs) else 0
            if lat + k < self.ilat_max and lon - k >= self.ilon_min:
                iwv_xs = self.dataset[lat + k, lon]
                differences.append(iwv_x-iwv_xs) if not np.isnan(iwv_xs) else 0
            if lat - k >= self.ilat_min and lon + k < self.ilon_max:
                iwv_xs = self.dataset[lat - k, lon - k]
                differences.append(iwv_x-iwv_xs) if not np.isnan(iwv_xs) else 0
            if lat + k < self.ilat_max and lon + k < self.ilon_max:
                iwv_xs = self.dataset[lat + k, lon + k]
                differences.append(iwv_x-iwv_xs) if not np.isnan(iwv_xs) else 0
        return differences

    def estimate_beta_for_area(self, lat_min, lat_max, lon_min, lon_max):
        s_diwv_dict = {}
        s_diffs_dict = {}
        for k in range(1,4):
            # 0.25 deg between each measurement and each 0.25deg = 27.75km
            s_adjacent, s_diag = round(k * 0.25 * 27.75, 4), round(k * 0.25 * np.sqrt(2) * 27.75, 4)
            adjacent_differences, diag_differences = [], []
            s_diwv_dict[s_adjacent] , s_diwv_dict[s_diag] = [], []
            for lat in range(lat_min, lat_max):
                for lon in range(lon_min, lon_max):
                    adjacent_differences.extend(self.get_adjacent_differences(lat, lon, k))
                    diag_differences.extend(self.get_diag_differences(lat, lon, k))

            d_iwv_hat = self.estimate_turbulence_struct(adjacent_differences)
            s_diwv_dict[s_adjacent] = d_iwv_hat
            # print(f"s = {round(s_adjacent, 4)}, and d_iwv(s) = {d_iwv_hat}")

            d_iwv_hat = self.estimate_turbulence_struct(diag_differences)
            s_diwv_dict[s_diag] = d_iwv_hat
            # print(f"s = {round(s_diag, 4)}, and d_iwv(s) = {d_iwv_hat}")

            s_diffs_dict[s_adjacent] = adjacent_differences
            s_diffs_dict[s_diag] = diag_differences
        self.area_sdiffs_dict[(lat_min, lon_min)] = s_diffs_dict
        self.area_sdiwv_dict[(lat_min, lon_min)] = s_diwv_dict

        beta, c_iwv = self.estimate_beta(s_diwv_dict)
        if self.verbose and not np.isnan(beta):
            print(f"latitude range = [{lat_idx_to_deg(lat_min)},{lat_idx_to_deg(lat_max)}]")
            print(f"longitude range = [{lon_idx_to_deg(lon_min)},{lon_idx_to_deg(lon_max)}]")
            print(f"beta = {round(beta, 4)}, and C_iwv(s) = {c_iwv}")
        return beta, c_iwv

    def run_algorithm(self):
        latlon_betaciwv_dict = {}

        for i in range((self.ilat_max-self.ilat_min)//self.area-1):
            for j in range((self.ilon_max-self.ilon_min)//self.area-1):
                beta, c_iwv = self.estimate_beta_for_area(self.ilat_min+i*self.area, \
                                                   self.ilat_min+(i+1)*self.area,\
                                                   self.ilon_min+j*self.area, \
                                                   self.ilon_min+(j+1)*self.area)
                latlon_betaciwv_dict[(self.ilat_min+i*self.area, self.ilon_min+j*self.area)] = (beta, c_iwv)
        self.ilatilon_betaciwv_dict = latlon_betaciwv_dict
        return latlon_betaciwv_dict

    def plot_s_vs_diwv_and_divw_hat(self, lat_d, lon_d):
        lat = lat_deg_to_idx(lat_d)
        lon = lon_deg_to_idx(lon_d)
        plt.figure()
        plt.title(r's vs $D_{IWV}(s)$ and $\hat{D}_{IWV}(s)$ for '+str(lat_d)+', '+str(lon_d))
        plt.ylabel(r'$D_{IWV}(s)$')
        plt.xlabel("s")
        s = list(self.area_sdiwv_dict[lat, lon].keys())
        plt.plot(s, self.area_sdiwv_dict[lat, lon].values(), marker='o')

        beta, ciwv = self.ilatilon_betaciwv_dict[(lat, lon)]
        plt.plot(s, ciwv*(s**beta))
        plt.legend([r'$D_{iwv}(s)$', r'$C_{iwv}s^{\beta}$'])
        plt.show()

    def plot_s_vs_diwv(self, lat_d, lon_d):
        lat = lat_deg_to_idx(lat_d)
        lon = lon_deg_to_idx(lon_d)
        plt.figure()
        plt.title(r's vs $D_{IWV}(s)$ for '+str(lat_d)+', '+str(lon_d))
        plt.ylabel(r'$D_{IWV}(s)$')
        plt.xlabel("s")
        s = list(self.area_sdiwv_dict[lat, lon].keys())
        plt.plot(s, self.area_sdiwv_dict[lat, lon].values(), marker='o')
        plt.show()

    def plot_s_vs_diffs_diwv(self, lat_d, lon_d):
        lat = lat_deg_to_idx(lat_d)
        lon = lon_deg_to_idx(lon_d)
        plt.figure()
        plt.title(r's (km) vs $IWV(x)-IWV(x+s)$ for '+str(lat_d)+', '+str(lon_d))
        plt.ylabel(r'$IWV(x)-IWV(x+s)$')
        plt.xlabel("s")

        s_s = list(self.area_sdiffs_dict[lat, lon].keys())
        s_pts = []
        [s_pts.extend([s]*len(self.area_sdiffs_dict[lat, lon][s])) for s in s_s]
        diff_pts = []
        [diff_pts.extend(self.area_sdiffs_dict[lat,lon][s]) for s in s_s]

        plt.scatter(s_pts, diff_pts, marker='+')
        plt.show()