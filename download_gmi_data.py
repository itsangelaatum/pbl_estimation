from ftplib import FTP
from time import sleep

def get_gmi_data(year='23', mon=range(1,12)):
    # Replace these with your FTP server details
    ftp_host = 'ftp.remss.com'
    ftp_user = 'awangs@umich.edu'
    ftp_passwd = 'awangs@umich.edu'

    # Create an FTP instance and connect to the server
    ftp = FTP(ftp_host)
    ftp.login(user=ftp_user, passwd=ftp_passwd)

    # Change to the directory where you want to pull data from
    year = '23'
    for m in mon:
        m = '0'+str(m) if m<10 else str(m)
        directory = '/gmi/bmaps_v08.2/y20'+year+'/m'+m+'/'
        ftp.cwd(directory)

        # List files and directories in the current directory
        file_list = ftp.nlst()

        for file_name in file_list:
            if 'd3d' in file_name or len(file_name)<18:
                continue
            # Download a file from the FTP server
            file_to_download = file_name  # Replace with the name of the file you want to download
            local_filename = file_name  # Local file name to save the downloaded file
            print(file_name)

            with open(local_filename, 'wb') as local_file:
                ftp.retrbinary('RETR ' + file_to_download, local_file.write)

    # Close the FTP connection
    ftp.quit()