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
import re
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

    output_path = os.path.realpath(os.path.join(os.path.dirname(
                                                    os.path.abspath(__file__)),
                                               'data/'))

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

        click.echo('Splitting file into chunks..')

        # We divide the observations files into 100 (20) roughly
        # equally sized chunks and use an index later for lookups
        # This greatly increases grep speed

        if filename == 'NumObs.txt.gz':
            N = 100  # number of chunks
            prefix = 'num'
            position = (0, 5)
        else:
            N = 20
            prefix = 'unnum'
            position = (5, 12)

        subprocess.call(['split', '-n', f'l/{N}', output_path,
                         f'{os.path.dirname(output_path)}/{prefix}'])

        # Now create lookup index
        click.echo('Creating lookup index..')
        index_path = f'data/index{prefix}.txt'
        index_path = os.path.realpath(os.path.join(os.path.dirname(__file__),
                                                   index_path))

        # To create the index, we check the first and last
        # designations in each file
        chunks_path = os.path.realpath(os.path.dirname(output_path))
        chunks = [os.path.join(chunks_path, c) for c in os.listdir(chunks_path)
                  if c.startswith(prefix)]

        with open(index_path, 'w') as index:
            for chunk in chunks:
                with open(chunk, 'r') as f:
                    first_desi = f.readline()[position[0]:position[1]]
                    for line in f:
                        last_desi = line[position[0]:position[1]]
                index.write(f'{chunk} {first_desi} {last_desi}\n')

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
        index_path = f'data/indexnum.txt'
        index_path = os.path.realpath(os.path.join(os.path.dirname(__file__),
                                                   index_path))

        to_grep = []  # most likely just one file, but could be two

        with open(index_path, 'r') as index:

            for line in index:

                filename, fnumber, lnumber = line.split()

                # Check if first character needs to be converted
                if fnumber[0].isalpha():
                    fnumber = str(_letter_to_number(fnumber[0])) \
                                  + str(fnumber[1:])
                fnumber = int(fnumber)

                if lnumber[0].isalpha():
                    lnumber = str(_letter_to_number(lnumber[0])) \
                                  + str(lnumber[1:])
                lnumber = int(lnumber)

                if fnumber <= number <= lnumber:
                    to_grep.append(filename)

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
        index_path = f'data/indexunnum.txt'
        index_path = os.path.realpath(os.path.join(os.path.dirname(__file__),
                                                   index_path))

        to_grep = []  # most likely just one file, but could be two

        # Convert passed designation
        check_for = Names.to_packed(designation)

        with open(index_path, 'r') as index:
            for line in index:

                filename, fdesi, ldesi = line.split()

                if fdesi <= check_for <= ldesi:
                    to_grep.append(line.split()[0])

        ident = ' ' * 5 + check_for

    else:
        click.echo(f'Need to provide either SSO number or designation of '
                   f'unnumbered minor planet')

    # Check if data exists and is not outdated
    if not os.path.isfile(to_grep[0]):
        click.echo(f'Could not find observations in path. '
                   f'Consider running "mpc retrieve" and getting a coffee.')
        sys.exit()

    if (time() - os.path.getmtime(to_grep[0])) / (3600*24) >= 30:
        click.echo('\nObservation files are older than one month. '
                   'Consider running "mpc retrieve" and getting a coffee.\n')

    if ',' in band:
        band = f'[{band.replace(",", "")}]'

    grep_string = r'\|'.join([f'^{ident}.*{band}......{obscode}$'
                              for obscode in observatory.split(',')])

    # Build grep command: The identifier is at the beginning of the line,
    # the observatory (if provided) at the end. At position 70, we should have
    # the band
    output = []
    for file_ in to_grep:
        grep = subprocess.Popen(['grep', grep_string, file_],
                                stdout=subprocess.PIPE)
        for line in grep.stdout:
            output.append(line)

    # If raw output is requested, echo to console or write to file
    if number:
        path_output = f'{number}.csv'
    elif designation:
        clean_designation = re.sub(r'[^\w]', '', designation)
        path_output = f'{clean_designation}.csv'

    if raw:
        if not csv:
            for obs in output:
                print(obs.decode().strip('\n'))
            sys.exit()
        else:
            with open(path_output, 'w') as out:
                for obs in output:
                    out.write(obs.decode())
            sys.exit()

    # Else, parse output and repeat
    parsed = pd.DataFrame(columns=['number', 'desig', 'discovery', 'note1',
                                   'note2', 'epoch', 'ra', 'dec', 'mag',
                                   'band', 'observatory'])

    for obs in output:
        obs = obs.decode().strip('\n')
        observation = _parse_observation(obs)

        if observation is False:
            continue

        parsed = parsed.append(observation,
                               ignore_index=True)

    if parsed.empty:
        parsed.to_csv(path_output, index=False)
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
        parsed.to_csv(path_output, index=False)
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
    ident = obs['number'].values[0]  # What SSO are we talking about

    if not isinstance(ident, (int, float, np.int64)):
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
