#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
    Author: Max Mahlke
    Date: 03 December 2019

    Extracts SSO observations from the complete observational
    sets of the MPC

    Call as:	python extract_observations.py AST_IDENTIFIER
'''

import os
import string
import subprocess
import sys
from time import time
from urllib.request import urlretrieve

from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
import click
import numpy as np
import pandas as pd
from tqdm import tqdm


@click.group()
def cli():
    pass


@cli.command()
@click.option('--unnumbered-only', is_flag=True,
              help='Only retrieve observations of unnumbered minor planets')
@click.option('--numbered-only', is_flag=True,
              help='Only retrieve observations of numbered minor planets')
def retrieve(unnumbered_only, numbered_only):
    ''' Retrieve observations to file
    '''

    for filename in ['NumObs.txt.gz', 'UnnObs.txt.gz']:

        if filename == 'NumObs.txt.gz':
            if unnumbered_only:
                continue
            click.echo(f'\nRetrieving observations of '
                       f'numbered minor planets..\n')

        if filename == 'UnnObs.txt.gz':
            if numbered_only:
                continue
            click.echo(f'\nRetrieving observations of '
                       f'unnumbered minor planets..\n')

        # Download
        output_path = _retrieve(filename)

        # Gunzip
        click.echo('Gunzipping..')
        subprocess.call(['gunzip', output_path])
        output_path = 'data/' + os.path.splitext(filename)[0]

        # Split files into chunks for faster grepping later
        if filename == 'NumObs.txt.gz':
            # Split the numbered observations by 100k
            click.echo('Splitting file into chunks..')

            # Firs 100k start with digit
            subprocess.call(['grep', f'^[0-9]', output_path],
                            stdout=open(f'{os.path.dirname(output_path)}/'
                                        f'0_99999.txt', 'w'))
            # Up to 620k, the two leading digits are endoced in ascii letters
            subprocess.call(['grep', f'^[A-J]', output_path],
                            stdout=open(f'{os.path.dirname(output_path)}/'
                                        f'100000_199999.txt', 'w'))
            subprocess.call(['grep', f'^[K-T]', output_path],
                            stdout=open(f'{os.path.dirname(output_path)}/'
                                        f'200000_299999.txt', 'w'))
            subprocess.call(['grep', f'^[T-Za-f]', output_path],
                            stdout=open(f'{os.path.dirname(output_path)}/'
                                        f'300000_399999.txt', 'w'))
            subprocess.call(['grep', f'^[g-p]', output_path],
                            stdout=open(f'{os.path.dirname(output_path)}/'
                                        f'400000_499999.txt', 'w'))
            subprocess.call(['grep', f'^[q-z]', output_path],
                            stdout=open(f'{os.path.dirname(output_path)}/'
                                        f'500000_620000.txt', 'w'))
            # Now, we're swtiching to ~=620k, and base64 encoding
            # second character is x*64**3
            subprocess.call(['grep', f'^~0', output_path],
                            stdout=open(f'{os.path.dirname(output_path)}/'
                                        f'620001_882144.txt', 'w'))
            subprocess.call(['grep', f'^~1', output_path],
                            stdout=open(f'{os.path.dirname(output_path)}/'
                                        f'620001_882144.txt', 'w'))
            subprocess.call(['grep', f'^~2', output_path],
                            stdout=open(f'{os.path.dirname(output_path)}/'
                                        f'882145_144288.txt', 'w'))

        if filename == 'UnnObs.txt.gz':
            # Split the unnumbered observations by century
            click.echo('Splitting file into chunks..')
            for year in ['18XX', '19XX', '20XX']:
                letter = string.ascii_uppercase[int(year[:2]) - 10]
                subprocess.call(['grep', f'^\ *{letter}', output_path],
                                stdout=open(f'{os.path.dirname(output_path)}/'
                                            f'{year}.txt', 'w'))

        # Remove cached observations file
        os.remove(output_path)


def _retrieve(filename):
    ''' Downloads the observations from MPC to the
    path/to/module/data directory

    :filename: str - filename at MPC
    :returns: str - absolute path to file output
    '''

    url = 'https://minorplanetcenter.net/iau/ECS/MPCAT-OBS/'
    output_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                              'data/')
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, filename)

    # Download observations
    with ProgressBar(unit='B', unit_scale=True,
                     miniters=1, desc=filename.split('.')[0]) as t:
        urlretrieve(url + filename, filename=output_path,
                    reporthook=t.update_to)
    return output_path


class ProgressBar(tqdm):
    # Simple progress bar for data download
    # https://github.com/tqdm/tqdm#hooks-and-callbacks
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)


@click.option('--number', type=int,
              help='Select observations by asteroid number')
@click.option('--designation', type=str,
              help='Select observations by asteroid number')
@click.option('--band', type=str, default='*',
              help='Select observations by asteroid number')
@click.option('--observatory', type=str, default='*',
              help=f'Select observations by observatory code.')
@click.option('--csv', is_flag=True,
              help=f'Save output to csv file')
@click.option('--raw', is_flag=True,
              help=f'Print original 80-column format')
@cli.command()
def obs(number, designation, band, observatory, csv, raw):
    ''' Echo MPC observations
    '''

    path_data = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                             'data/')

    # Locate the correct file chunk to grep
    if number:
        if number <= 99999:
            filename = '0_99999.txt'
        elif 100000 <= number <= 199999:
            filename = '100000_199999.txt'
        elif 200000 <= number <= 299999:
            filename = '200000_299999.txt'
        elif 300000 <= number <= 399999:
            filename = '300000_399999.txt'
        elif 400000 <= number <= 499999:
            filename = '400000_499999.txt'
        elif 500000 <= number <= 620000:
            filename = '500000_620000.txt'
        elif 620001 <= number <= 882144:
            filename = '620001_882144.txt'
        elif 882145 <= number <= 144288:
            filename = '882145_144288.txt'

        # Build the identifier. It should always equal 12 characters
        if number <= 99999:
            # The MPC identifier is just the zero-padded number
            ident = str(number).zfill(5)
            ident = ident.ljust(12, '.')
        else:
            # We need to convert the two first digits into a letter
            ident = str(number)
            ident = _number_to_letter(ident[:2]) + ident[2:]
            ident = ident.ljust(12, '.')

    elif designation:
        sys.exit()
        filename = ''
    else:
        click.echo(f'Need to provide either SSO number or designation of '
                   f'unnumbered minor planet')

    path_mpc = os.path.join(path_data, filename)

    if not os.path.isfile(path_mpc):
        click.echo(f'Could not find observations in path. '
                   f'Consider running "mpc retrieve" and getting a coffee.')
        sys.exit()

    if (time() - os.path.getmtime(path_mpc)) / (3600*24) >= 30:
        click.echo('\nObservation files are older than one month. '
                   'Consider running "mpc retrieve" and getting a coffee.\n')

    # Build grep command: The identifier is at the beginning of the line,
    # the observatory (if provided) at the end. At position 70, we should have
    # the band
    grep = subprocess.Popen(['grep', f'^{ident}.*{band}......{observatory}$',
                             path_mpc], stdout=subprocess.PIPE)

    # If raw output is requested, echo to console or write to file
    if raw:
        if not csv:
            for obs in grep.stdout:
                print(obs.decode().strip())
            sys.exit()
        else:
            with open(ident.strip('.') + '.csv', 'w') as out:
                for obs in grep.stdout:
                    out.write(obs.decode())
            sys.exit()

    # Else, parse output and repeat
    parsed = pd.DataFrame(columns=['number', 'desig', 'discovery', 'note1',
                                   'note2', 'epoch', 'ra', 'dec', 'mag',
                                   'band', 'observatory'])

    for obs in grep.stdout:
        obs = obs.decode().strip()

        # Satellite data is skipped because the coordinate treatment
        # is not trivial
        if obs[-3:] in ['C51']:
            continue

        parsed = parsed.append(_parse_observation(obs),
                               ignore_index=True)

    if not csv:
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', None)
        click.echo(parsed)
    else:
        parsed.to_csv(ident.strip('.') + '.csv', index=False)
    return obs


def _number_to_letter(number):
    ''' Converts a two-digit number into a letter,
    following the MPC packed format '''

    number = int(number)

    assert number >= 10, f'Number {number} should be >= 10'

    letters = [*string.ascii_uppercase, *string.ascii_lowercase]
    return letters[number - 10]


def _letter_to_number(letter):
    ''' Converts a letter into a two-digit number,
    following the MPC packed format '''

    letters = [*string.ascii_uppercase, *string.ascii_lowercase]
    assert letter in letters, f'Letter {letter} not recognised'

    letters = np.array(letters)
    return np.where(letters == letter)[0][0] + 10


def _parse_observation(line):
    ''' Converts the 80-column MPC format to CSV data column

    :line: str - MPC formatted observation line
    :returns: pd.Series - observation in more practical format
    '''

    # ------
    # Unpack either the number or the designation

    # Convert the number
    number = line[:5]

    if not number.isspace():
        # Check if first character needs to be converted
        if number[0].isalpha():
            number = _letter_to_number(number[0]) + number[1:]

        number = int(number)
        desig = ''

    else:
        # Unpack the designation
        desig = ''
        number = ''

    # ------
    discovery = line[12] if line[12] != ' ' else ''
    note1 = line[13] if line[13] != ' ' else ''
    note2 = line[14] if line[14] != ' ' else ''

    # ------
    # Convert the epoch to MJD
    epoch = Time(line[15:25].replace(' ', '-')).mjd  # day to MJD
    epoch = f'{int(epoch)}.{line[26:32]}' # add decimal day

    # ------
    # Convert ra and dec to degree
    coord = SkyCoord(ra=line[32:44], dec=line[44:57],
                     unit=(u.hourangle, u.deg),
                     frame='icrs')
    ra = coord.ra.deg
    dec = coord.dec.deg

    # ------
    if not line[65:69].isspace():
        mag = float(line[65:69])
    else:
        mag = ''
    band = line[70]
    observatory = line[-3:]

    obs = pd.Series({'number': number, 'desig': desig, 'discovery': discovery,
                     'note1': note1, 'note2': note2, 'epoch': epoch, 'ra': ra,
                     'dec': dec, 'mag': mag, 'band': band,
                     'observatory': observatory})
    return obs

    # # convert asteroid designations
    # # old designation style, e.g.: 1989AB
    # ident = data['pdesig'][0]
    # if isinstance(ident, np.ma.masked_array) and ident.mask:
        # ident = ''
    # elif (len(ident) < 7 and ident[:4].isdigit() and
            # ident[4:6].isalpha()):
        # ident = ident[:4]+' '+ident[4:6]
    # # Palomar Survey
    # elif 'PLS' in ident:
        # ident = ident[3:] + " P-L"
    # # Trojan Surveys
    # elif 'T1S' in ident:
        # ident = ident[3:] + " T-1"
    # elif 'T2S' in ident:
        # ident = ident[3:] + " T-2"
    # elif 'T3S' in ident:
        # ident = ident[3:] + " T-3"
    # # standard MPC packed 7-digit designation
    # elif (ident[0].isalpha() and ident[1:3].isdigit() and
          # ident[-1].isalpha() and ident[-2].isdigit()):
        # yr = str(conf.pkd.find(ident[0]))+ident[1:3]
        # let = ident[3]+ident[-1]
        # num = str(conf.pkd.find(ident[4]))+ident[5]
        # num = num.lstrip("0")
        # ident = yr+' '+let+num
    # data.add_column(Column([ident]*len(data), name='desig'),
                    # index=1)
    # data.remove_column('pdesig')

    # elif all([o['object_type'] != 'M' for o in src]):
        # # comets
        # data = ascii.read("\n".join([o['original_record']
                                     # for o in src]),
                          # format='fixed_width_no_header',
                          # names=('number', 'comettype', 'desig',
                                 # 'note1', 'note2', 'epoch',
                                 # 'RA', 'DEC', 'mag', 'phottype',
                                 # 'observatory'),
                          # col_starts=(0, 4, 5, 13, 14, 15,
                                      # 32, 44, 65, 70, 77),
                          # col_ends=(3, 4, 12, 13, 14, 31,
                                    # 43, 55, 69, 70, 79))

        # # convert comet designations
        # ident = data['desig'][0]

        # if (not isinstance(ident, (np.ma.masked_array,
                                   # np.ma.core.MaskedConstant))
                # or not ident.mask):
            # yr = str(conf.pkd.find(ident[0]))+ident[1:3]
            # let = ident[3]
            # # patch to parse asteroid designations
            # if len(ident) == 7 and str.isalpha(ident[6]):
                # let += ident[6]
                # ident = ident[:6] + ident[7:]
            # num = str(conf.pkd.find(ident[4]))+ident[5]
            # num = num.lstrip("0")
            # if len(ident) >= 7:
                # frag = ident[6] if ident[6] != '0' else ''
            # else:
                # frag = ''
            # ident = yr+' '+let+num+frag
            # # remove and add desig column to overcome length limit
            # data.remove_column('desig')
            # data.add_column(Column([ident]*len(data),
                                   # name='desig'), index=3)
    # else:
        # raise ValueError(('Object type is ambiguous. "{}" '
                          # 'are present.').format(
                              # set([o['object_type'] for o in src])))

    # # convert dates to Julian Dates
    # dates = [d[:10].replace(' ', '-') for d in data['epoch']]
    # times = np.array([float(d[10:]) for d in data['epoch']])
    # jds = Time(dates, format='iso').jd+times
    # data['epoch'] = jds

    # # convert ra and dec to degrees
    # coo = SkyCoord(ra=data['RA'], dec=data['DEC'],
                   # unit=(u.hourangle, u.deg),
                   # frame='icrs')
    # data['RA'] = coo.ra.deg
    # data['DEC'] = coo.dec.deg

    # # convert Table to QTable
    # data = QTable(data)
    # data['epoch'].unit = u.d
    # data['RA'].unit = u.deg
    # data['DEC'].unit = u.deg
    # data['mag'].unit = u.mag

    # return line
