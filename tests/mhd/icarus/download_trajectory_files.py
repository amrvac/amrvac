# Usage of this python script:
#
# python3 download_trajectory_files.py -o 'path_to_the_output_folder' --date YYYY-MM-DD
# or
# python3 download_trajectory_files.py -o 'path_to_the_output_folder' -i 'path_to_the_amrvac_parameter_file'
#
# -o is not needed, in which case the files are downloaded to the folder './orbit/YYYY-MM-DD/'
#
# the script also produces a file 'download-TVJAFE-...' which contains information about the result of downloading


# Extends script by Jonathan Dan
__email__ = "rdm at kuleuven.be"
__license__ = "CC BY-SA 4.0 - https://creativecommons.org/licenses/by-sa/4.0/"


import hashlib
import json
import logging
import os
from pathlib import Path
import requests
from requests.adapters import HTTPAdapter, Retry
import argparse
import logging
import datetime
import os
import urllib3
import pandas as pd
from scipy.io import FortranFile

DATAVERSE_URL = ""
API_KEY = ""
DOI = ""
VERSION = ""
OUTPUT_DIR = "./datasets"
LOG_PATH = ''
DATE = ''


def config_logger():
    global LOG_PATH
    log_format = logging.Formatter(fmt='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S')

    filename = DOI.split('/')[-1] 
    LOG_PATH = os.path.join(
        '.', 
        f'download-{filename}-{datetime.datetime.today().strftime("%Y%m%d%H%M%S")}.log'
        )
    file_handler = logging.FileHandler(filename=LOG_PATH, mode='a')
    file_handler.setFormatter(log_format)
    file_handler.setLevel(logging.INFO)
    logger = logging.getLogger()
    logger.addHandler(file_handler)
    logger.setLevel(logging.INFO)


def fix_doi(doi):
    if 'doi.org' in doi:
        doi = 'doi:' + '/'.join(doi.split('/')[-2:])
    return doi

def fix_version(version):
    try:
        return str(float(version))
    except:
        logging.error(f'Invalid version :{version}')
        
def fix_date(date, input_file):
    if date is None:
        with open(input_file, 'r') as parameter_file:
            pd_date = parameter_file.read().split('magnetogram_time')[-1]
            pd_date = pd_date.strip('= ')
            pd_date = pd_date[:11] + pd_date[0]
            #print(pd_date)
        return pd_date
    else:
        return date      
        
def fix_output(output_path, date_string):
    if output_path is None:
        return os.path.join('.', 'orbit/', date_string.strip('\'\"') , '')
    else:
        return output_path
    

def get_args():
    global DATAVERSE_URL
    global API_KEY
    global DOI 
    global VERSION
    global OUTPUT_DIR
    global DATE
    parser = argparse.ArgumentParser()
    parser._optionals.title = 'arguments'
    parser.add_argument('-o', '--output', required=False, help='the output directory.')
    parser.add_argument('--host', required=False, default='https://rdr.kuleuven.be/', help='Dataverse URL [REQUIRED]')
    parser.add_argument('--key', required=False,  help='Dataverse API token')
    parser.add_argument('--dataset', required=False, default='https://doi.org/10.48804/TVJAFE' , help='dataset DOI in the format of  doi:......../...... or https://doi.org/......../......, [REQUIRED]')
    parser.add_argument('--version', required=False, default='1.0', help='dataset version number [REQUIRED]')
    
    date_group = parser.add_mutually_exclusive_group(required=True)
    date_group.add_argument('--date', help='date of the (first) magnetogram used in the simulation in the format YYYY-MM-DD')
    date_group.add_argument('-i', '--input', help='parameter file where the magnetogram time is stored')
    
    args = parser.parse_args()
    DATAVERSE_URL = args.host
    API_KEY = args.key
    DOI = fix_doi(args.dataset) 
    VERSION = fix_version(args.version)
    DATE = fix_date(args.date, args.input)
    OUTPUT_DIR = fix_output(args.output, DATE) #or os.path.join(*OUTPUT_DIR.split('/'), DOI.replace('doi:', 'doi-').replace('/', '-'), f'v{VERSION}')
    DATE = pd.Timestamp(DATE)

urllib3.disable_warnings()
if __name__=='__main__':
    get_args()
    config_logger()
    logging.info(f'Dataverse: {DATAVERSE_URL}')
    logging.info(f'Dataset: {DOI}')
    logging.info(f'Version: {VERSION}')
    logging.info(f'Output folder: {OUTPUT_DIR}')
    logging.info(f'Date: {DATE}')
    print(f'Dataverse: {DATAVERSE_URL}')
    print(f'Dataset: {DOI}')
    print(f'Version: {VERSION}')
    print(f'Output folder: {OUTPUT_DIR}')
    print(f'Date: {DATE}')
    print(f'Downloading dataset...')

    s = requests.Session()
    retries = Retry(total=5,
                backoff_factor=1)
    s.mount('http://', HTTPAdapter(max_retries=retries))
    s.mount('https://', HTTPAdapter(max_retries=retries))

    # Get list of FILES
    url = DATAVERSE_URL + '/api/datasets/:persistentId/versions/' + VERSION + '/files'
    respFileList = s.get(url, params={'persistentId': DOI}, verify=False)
    failed_files = []
    restricted_files = []
    hash_difference = []
    downloaded = []
    present = []
    satellite_list = ['earth', 'mars', 'mercury', 'venus', 'sta', 'stb', 'psp', 'SolO', 'mpo'] # 'jupiter', 'juno']
    
    if respFileList.ok:
        fileList = json.loads(respFileList.content)
        
        # Download files
        files = fileList['data']
        total=len(files)
        for i, file in enumerate(files):
            
            # Check if the file should be downloaded based on DATE and the name of the satellite
            filename = file['dataFile']['filename']
            if ('readme' in filename.lower()): 
                continue
            satellite, start_date, end_date = filename.split('__')
            end_date = end_date.split('.')[0]
            
            start_date = start_date.replace('_', '-')
            end_date = end_date.replace('_', '-')
            
            start_date=pd.Timestamp(start_date)
            end_date=pd.Timestamp(end_date)
            # print(start_date, DATE)
            
            if(start_date+pd.Timedelta(15, unit='d') < DATE < end_date-pd.Timedelta(60, unit='d') and (satellite in satellite_list)):
            
                # Create directory structure
                file_path = OUTPUT_DIR.split('/') + file.get('directoryLabel', '').split('/') 
                description = f"{'/'.join([file.get('directoryLabel'),filename]) if file.get('directoryLabel') else filename}"
                filename = os.path.join(*(file_path + [filename]))
                logging.info('Downloading file: ' + description + f"({i+1}/{total})")
                print('Downloading file: ' + description + f"({i+1}/{total})")
                Path(os.path.join(*file_path)).mkdir(parents=True, exist_ok=True)

                # if file exists skip file:
                # if os.path.exists(filename):
                    # md5Hash = hashlib.md5(open(filename, 'rb').read()).hexdigest()
                    # if md5Hash == file['dataFile']['md5']:
                        # logging.info('Skipping already downloaded file {}'.format(filename))
                        # present.append(file['dataFile']['md5'])
                        # downloaded.append(file['dataFile']['md5'])
                        # continue

                # Download file
                try:
                    url = DATAVERSE_URL + '/api/access/datafile/' + str(file['dataFile']['id'])
                    if API_KEY:
                        headers = {'X-Dataverse-key': API_KEY}
                        respFile = s.get(url, headers=headers)
                    else:
                        respFile = s.get(url)
                except Exception as e:
                    logging.error('{} - {}'.format(file['dataFile']['id'], e))
                if respFile.ok:
                    # Chec md5 hash
                    md5Hash = hashlib.md5(respFile.content).hexdigest()
                    if md5Hash == file['dataFile']['md5']:
                        # Write file to disk
                        with open(filename, 'wb') as edfFile:
                            edfFile.write(respFile.content)
                        downloaded.append(file['dataFile']['md5'])
                        
                        # Add the name of the file to the list of downloaded files
                        with open(os.path.join(OUTPUT_DIR, 'satellite_list.txt'), 'a') as list_downloaded_satellites:
                            list_downloaded_satellites.write(description + '\n')
                        
                        
                    else:
                        logging.error('{} - MD5 hash difference'.format(file['dataFile']['id']))
                elif respFile.status_code == 403:
                    failed_files.append(file['dataFile']['filename'])
                    restricted_files.append(file['dataFile']['filename'])
                    logging.error(f"{file['dataFile']['id']} - Failed to download file. {respFile.status_code}: access to file is restricted.")
                else:
                    failed_files.append(file['dataFile']['filename'])
                    try:
                        message =  json.loads(respFile.text).get('message')
                    except:
                        message = respFile.text
                    logging.error(f"{file['dataFile']['id']} - Failed to download file. {respFile.status_code}: {message}")
        logging.warning(f"{len(present)} of {len(files)} files were already present.")
        logging.warning(f"{len(failed_files)} of {len(files)} files failed.")
        logging.warning(f"{len(restricted_files)} skipped due to restricted access.")
        logging.info(f"{len(downloaded)} of {len(files)} files downloaded correctly")
        if failed_files:
            print(f"{len(failed_files)} downloads failed.")
        if restricted_files:
            print(f"{len(restricted_files)} were skipped due to restricted access.")
        print(f"{len(downloaded)} of {len(files)} files downloaded correctly.")
        print(f'See log file {LOG_PATH} for details.')
    else:
        try:
            message =  json.loads(respFileList.text).get('message')
        except:
            message = respFileList.text
            
        logging.error(f'Failed to retrieve dataset file list. {respFileList.status_code}: {message}')
