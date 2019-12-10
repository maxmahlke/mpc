# mpc

Local queries of MPC observations. When disk-space is no issue and remote-queries fail too often. Else, use `astroquery.mpc`

### Install

Clone the git repository or download the zip. Then run 

` $ [sudo] pip install --editable .`

to install it. An executable `mpc` is added to your path and can be called system-wide.

### Usage

Running `mpc retrieve` will download the latest observations to the package directory. Unzipped, this is more than 20GB of text files.
You can specify `--numbered-only` or `--unnumbered-only` to only update the respective observations.

Once the download has finished, run `mpc obs` to grep through the observation files. Target asteroids can be specified by `--number` or  `--designation`. By default, the observations are echoed to the terminal. Using `--csv` will save them to file in the current working directory. Using the `--raw` flag, the original MPC format is returned, else, the observations are converted to a more readable format.

