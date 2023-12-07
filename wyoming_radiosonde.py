"""
# What it is

**pysoundings** is just a simple little module designed to pull atmospheric sounding data from the [University of Wyoming, College of Engineering, Department of Atmospheric Science's website](http://weather.uwyo.edu/upperair/sounding.html). All you need is the station number and the date and this module will return a [pandas.DataFrame](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html) containing the sounding data. You may still want to go to the website to determine what stations are available, but there are data from a lot of stations in the United States.
"""

__docformat__ = "numpy"


from bs4 import BeautifulSoup
import requests
import pandas
import io
import warnings
import certifi
from datetime import datetime

def build_param_dict(stnm, dt):
    """
    Builds a dictionary containing the station number, year, month, and day/hour for the desired atmospheric sounding data.
    
    Parameters
    ----------
    stnm : string
        String of the station identifier, e.g. '70261' for PAFA, Fairbanks Int'l Airport.
    year : string
        String of the year, e.g. '2021'.
    month : string
        String of the month, e.g. '01'.
    day_hour : string
        String of the combined day and hour, e.g. '0100' for the first day of the month as '01', and for the beginning of the day in UTC as '00'.

    Returns
    -------
    param_dict : dict
        A dictionary containing the station number, year, month, and day/hour of the desired date and location.
    """
    param_dict = {'id': stnm, 'datetime': dt.strftime('%Y-%m-%d%%20%H:00:00')}

    return param_dict


def build_url(param_dict):
    """
    Builds the URL needed to query [University of Wyoming, College of Engineering, Department of Atmospheric Science's website](http://weather.uwyo.edu/upperair/sounding.html) to get the proper sounding data.
    
    Parameters
    ----------
    param_dict : dict
        A dictionary containing the station number, year, month, and day/hour of the desired date and location.
    
    Returns
    -------
    full_url : string
        String of the query URL with the proper date and location of the desired atmospheric sounding.
    """
    
    base_url = 'https://weather.uwyo.edu/cgi-bin/bufrraob.py?type=TEXT:LIST'

    full_url = base_url
    for key in param_dict.keys():
        full_url += '&' + key + '=' + param_dict[key]

    return full_url


def format_data(html_string):
    """
    Takes a string containing the html container that has the fixed width formatted atmospheric sounding data and stores it in a [pandas.DataFrame](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html), as well as extracts the units for the data.
    
    Parameters
    ----------
    html_string : string
        A string of the html container that holds the atmospheric sounding data.

    Returns
    -------
    data : DataFrame
        A [pandas.DataFrame](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html) containing the atmospheric data with a labeled header.
    units : string
        String of the units of the pandas DataFrame columns.
    """
    
    fwf = html_string.split('\n')[5:-2]
    header = html_string.split('\n')[2].split()
    units = html_string.split('\n')[3]

    data = pandas.DataFrame([x.split() for x in fwf], columns=header)

    return data, units


def pull_data(url):
    """
    This function makes an http request with the given URL to retrieve the atmospheric sounding data from the [University of Wyoming, College of Engineering, Department of Atmospheric Science's website](http://weather.uwyo.edu/upperair/sounding.html). It then puts the data into a [pandas.DataFrame](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html) and the units into a string.

    Parameters
    ----------
    url : string
        String of the url containing the information necessary to query the sounding database.

    Returns
    -------
    data : DataFrame
        A [pandas.DataFrame](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html) containing the atmospheric data with a labeled header.
    units : string
        String of the units of the pandas DataFrame columns.
    """
    # Path to certifi's CA bundle
    ca_bundle_path = certifi.where()
    print(url)

    page = requests.get(url, verify=False)
#     page = requests.get(url)
    soup = BeautifulSoup(page.content, 'lxml')
    data_container = soup.pre

    if data_container == None:
        warnings.warn('No data found for the given parameters')
        data, units = None, None

    else: 
        data, units = format_data(data_container.prettify())

    return data, units


def get_data(stnm, datetime):
    """
    This function calls all the necessary functions to build the url with the given date and station number, and then extracts the data from the loaded webpage from the [University of Wyoming, College of Engineering, Department of Atmospheric Science's website](http://weather.uwyo.edu/upperair/sounding.html) that contains the atmospheric sounding data.

    Parameters
    ----------
    stnm : string
        String of the station identifier, e.g. '70261' for PAFA, Fairbanks Int'l Airport.
    year : string
        String of the year, e.g. '2021'.
    month : string
        String of the month, e.g. '01'.
    day_hour : string
        String of the combined day and hour, e.g. '0100' for the first day of the month as '01', and for the beginning of the day in UTC as '00'.

    Returns
    -------
    data : DataFrame
        A [pandas.DataFrame](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html) containing the atmospheric data with a labeled header.
    units : string
        String of the units of the pandas DataFrame columns.

    Examples
    --------

    >>> stnm, year, month, day_hour = '70261', '2021', '01', '0100'
    >>> data, units = get_data(stnm, year, month, day_hour)
    >>> data
           PRES   HGHT  TEMP  DWPT  RELH  MIXR   DRCT  SKNT   THTA   THTE   THTV
    0    1000.0     93   NaN   NaN   NaN   NaN    NaN   NaN    NaN    NaN    NaN
    1     994.0    134 -23.3 -25.7  81.0  0.48   85.0   1.0  250.3  251.6  250.3
    2     992.0    149 -19.3 -20.9  87.0  0.73   80.0   1.0  254.4  256.5  254.6
    3     990.0    164 -18.9 -20.9  84.0  0.73   75.0   2.0  255.0  257.0  255.1
    4     988.0    179 -15.7 -17.6  85.0  0.98   70.0   2.0  258.3  261.1  258.5
    ..      ...    ...   ...   ...   ...   ...    ...   ...    ...    ...    ...
    153    10.0  30830 -59.7 -89.7   1.0  0.01  345.0  94.0  795.6  795.8  795.7
    154     9.9  30893 -60.1 -90.1   1.0  0.01  345.0  95.0  796.5  796.6  796.5
    155     9.6  31090 -59.8 -90.1   1.0  0.01  345.0  98.0  804.9  805.0  804.9
    156     9.1  31394 -59.3 -90.0   1.0  0.01  345.0  94.0  818.0  818.2  818.0
    157     8.0  32229 -57.9 -89.9   1.0  0.02    NaN   NaN  855.2  855.4  855.2
        [158 rows x 11 columns]
    >>> units
    '    hPa     m      C      C      %    g/kg    deg   knot     K      K      K '
    """
    
    param_dict = build_param_dict(stnm, datetime)
    url = build_url(param_dict)
    data, units = pull_data(url)

    return data, units


def save_data(csv_path, data, units):
    """
    Saves the DataFrame containing the atmospheric sounding data, and the accompanying units, and saves the data in a csv file.

    Parameters
    ----------
    csv_path : string
        String of the file path to save the data at.
    data : DataFrame
        A [pandas.DataFrame](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html) containing the atmospheric data with a labeled header.
    units : string
        String of the units of the pandas DataFrame columns.
    """
    
    with open(csv_path, 'a') as file:
        file.write('# units:' + units + '\n')
        data.to_csv(file)