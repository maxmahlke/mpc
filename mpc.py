#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
    Author: Max Mahlke
    Date: 03 December 2019

    Extracts SSO observations from the complete observational
    sets of the MPC

    Call as:    mpc --help
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
import requests
from sbpy.data import Names
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
        output_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   output_path)

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
              help='Select observations by band. Can be comma-separated list.')
@click.option('--observatory', type=str, default='*',
              help=f'Select observations by observatory code.'
                   f'Can be comma-separated list.')
@click.option('--csv', is_flag=True,
              help=f'Save output to csv file')
@click.option('--raw', is_flag=True,
              help=f'Print original 80-column format')
@click.option('--ephemerides', is_flag=True,
              help=f'Query asteorid ephemerides with Miriade')
@cli.command()
def obs(number, designation, band, observatory, csv, raw, ephemerides):
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
        for cent in ['18', '19', '20']:
            if designation.startswith(cent):
                filename = f'{cent}XX.txt'
                break
            ident = ' ' * 5 + Names.to_packed(designation)

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

    if ',' in band:
        band = f'[{band.replace(",", "")}]'

    grep_string = r'\|'.join([f'^{ident}.*{band}......{obscode}$'
                              for obscode in observatory.split(',')])

    # Build grep command: The identifier is at the beginning of the line,
    # the observatory (if provided) at the end. At position 70, we should have
    # the band
    grep = subprocess.Popen(['grep', grep_string, path_mpc],
                            stdout=subprocess.PIPE)

    # If raw output is requested, echo to console or write to file
    if raw:
        if not csv:
            for obs in grep.stdout:
                print(obs.decode().strip('\n'))
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
        obs = obs.decode().strip('\n')
        observation = _parse_observation(obs)

        if observation is False:
            continue

        parsed = parsed.append(observation,
                               ignore_index=True)

    if parsed.empty:
        click.echo('No observations found!')
        parsed.to_csv(ident.strip().strip('.') + '.csv', index=False)
        sys.exit()

    if ephemerides:
        parsed = _miriade_ephems_for_obs(parsed)
        if parsed is False:
            print(f'Miriade query failed for {ident}')
            sys.exit()

    if not csv:
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', None)
        click.echo(parsed)
    else:
        parsed.to_csv(ident.strip().strip('.') + '.csv', index=False)
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
            number = str(_letter_to_number(number[0])) + str(number[1:])

        number = int(number)
        desig = ''

    else:
        # Unpack the designation
        desig = Names.from_packed(line[5:12])
        number = ''

    # ------
    discovery = line[12] if line[12] != ' ' else ''
    note1 = line[13] if line[13] != ' ' else ''
    note2 = line[14] if line[14] != ' ' else ''

    # ------
    # Convert the epoch to MJD
    epoch = Time(line[15:25].replace(' ', '-')).mjd  # day to MJD
    epoch = f'{int(epoch)}.{line[26:32]}'  # add decimal day

    # ------
    # Convert ra and dec to degree
    try:
        coord = SkyCoord(ra=line[32:44], dec=line[44:57],
                         unit=(u.hourangle, u.deg),
                         frame='icrs')
    except ValueError:  # arises when satellite position is given. skip those
        return False
    ra = coord.ra.deg
    dec = coord.dec.deg

    # ------
    if not line[65:69].isspace():
        mag = float(line[65:70])
    else:
        mag = ''
    band = line[70]
    observatory = line[-3:]

    obs = pd.Series({'number': number, 'desig': desig, 'discovery': discovery,
                     'note1': note1, 'note2': note2, 'epoch': epoch, 'ra': ra,
                     'dec': dec, 'mag': mag, 'band': band,
                     'observatory': observatory})
    return obs


def _miriade_ephems_for_obs(obs):
    '''Gets asteroid ephemerides from IMCCE Miriade for MPC formatted
    observations of a single SSO

    :obs: pd.DataFrame - MPC observations
    :returns: pd.DataFrame - Input dataframe with ephemerides columns appended
                     False - If query failed somehow
    '''
    try:
        ident = obs['number'].values[0]  # What SSO are we talking about 
    except IndexError:
        ident = obs['desig'].values[0]

    obs = obs.sort_values('epoch')  # increases performance in Miriade
    obs = obs.reset_index()

    # ------
    # Query Miriade for phase angles
    url = 'http://vo.imcce.fr/webservices/miriade/ephemcc_query.php'

    params = {

        '-name': f'a:{ident}',
        '-mime': 'json',
        '-tcoor': '5',
        '-output': '--jul',
        '-tscale': 'UTC'
    }

    # Pass sorted list of epochs to speed up query
    # Have to convert them to JD
    epochs = [Time(e, format='mjd').jd for e in obs.epoch.astype(float)]
    files = {'epochs': ('epochs',
                        '\n'.join(['%.6f' % epoch
                                   for epoch in epochs]))}
    # Execute query
    try:
        r = requests.post(url, params=params, files=files, timeout=50)
    except requests.exceptions.ReadTimeout:
        return False
    j = r.json()

    # Read JSON response
    try:
        ephem = pd.DataFrame.from_dict(j['data'])
    except KeyError:
        return False

    obs = pd.merge(obs, ephem, how='outer', left_on=obs.index,
                   right_on=ephem.index)
    return obs
