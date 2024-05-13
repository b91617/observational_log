# %%
import argparse
from tqdm import tqdm
import ccdproc as ccdp
from astropy.io import fits
from astropy.time import Time
from pathlib import Path
import os
import csv

parser = argparse.ArgumentParser(description='Log observation parameters.')

parser.add_argument('-i', '--file-input',
                    default=os.getcwd(),
                    help='Input working folder.')
parser.add_argument('-n', '--header-no',
                    type=int,
                    default=0,
                    help='Input the header number of where IMAGETYP are in.')

args = parser.parse_args()

loc = args.file_input
no  = args.header_no

def log(root_dir, header_no):
    
    print(f'Current working directory: {root_dir}')

    folder_ls = [dirs for dirs in os.listdir(root_dir) if os.path.isdir(os.path.join(root_dir, dirs))]
    true_folder_ls = [dirs for dirs in folder_ls if not dirs.startswith('.')]

    with open('log.csv', 'w', newline='') as log:
        csv.writer(log).writerow(['File Name', 'Date', 'Time', 'MJD', 'Image Type', 'Filter', 'Exposure Time', 'X Bin', 'Y Bin', 'MAGZERO', 'MAGZPT', 'ZMAG'])
        
        for folder in true_folder_ls:
            working_dir = os.path.join(root_dir + '/' + folder) + '/'

            fits_files = ccdp.ImageFileCollection(location=Path(working_dir),
                                                  find_fits_by_reading=True,
                                                  glob_exclude='.*',
                                                  ext=header_no)
            imgtyp_set = set(fits_files.summary['imagetyp'])

            for image_type in imgtyp_set:
                fits_list = fits_files.files_filtered(imagetyp=image_type)
                
                for fits_name in tqdm(fits_list):
                    data = fits.open(working_dir + fits_name)

                    try:
                        obs_date = data[header_no].header['DATE-OBS']
                    except KeyError:
                        obs_date = '-'

                    try:
                        obs_time = data[header_no].header['TIME-OBS']
                    except KeyError:
                        obs_time = '-'
                    
                    try:
                        obs_jd = data[header_no].header['JD']
                        obs_mjd = Time(obs_jd, format='jd').mjd
                    except KeyError:
                        obs_mjd = '-'

                    try:
                        imgtyp = data[header_no].header['IMAGETYP']
                    except KeyError:
                        imgtyp = '-'

                    try:
                        obs_filter = data[header_no].header['FILTER']
                    except KeyError:
                        obs_filter = '-'
                    
                    try:
                        exp_time = data[header_no].header['EXPTIME']
                    except KeyError:
                        exp_time = '-'

                    try:
                        xbin = data[header_no].header['XBINNING']
                    except KeyError:
                        xbin = '-'

                    try:
                        ybin = data[header_no].header['YBINNING']
                    except KeyError:
                        ybin = '-'

                    try:
                        magzero = data[header_no].header['MAGZERO']
                    except KeyError:
                        magzero = '-'

                    try:
                        magzpt = data[header_no].header['MAGZPT']
                    except KeyError:
                        magzpt = '-'
                    
                    try:
                        zmag = data[header_no].header['ZMAG']
                    except KeyError:
                        zmag = '-'

                    csv.writer(log).writerow([fits_name, obs_date, obs_time, obs_mjd, imgtyp, obs_filter, exp_time, xbin, ybin, magzero, magzpt, zmag])
                    
                    data.close()

    print('DONE')

log(loc, no)