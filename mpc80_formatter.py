import pandas as pd
from astropy.coordinates import Angle
from astropy.time import Time
from os.path import isfile

import argparse

__all__ = ["format_obs80"]

def format_obs80(obj_id=None, time=None, ra=None, dec=None, mag=None, filter=None, obs_code="I11",
                 obs_file=None, columns={"time": "FieldMJD_TAI", "ra": "AstRA(deg)", "dec": "AstDec(deg)", 
                                         "mag": "observedTrailedSourceMag", "filter": "optFilter",
                                         "id": "hex_id", "obs_code": "obs_code"},
                 out_path=None):
    """Format a set of observations *of minor planets* into the Minor Planet Center 80 column format

    Note that this function is specifically designed to format observations of minor planets, I'm not sure if
    it applies to comets or natural satellites.

    Parameters
    ----------
    obj_id : :class:`numpy.ndarray`, optional
        The unique IDs (7 character strings) of the objects being observed
    time : :class:`numpy.ndarray`, optional
        The observation times in MJD
    ra : :class:`numpy.ndarray`, optional
        The right ascensions of the objects in degrees
    dec : :class:`numpy.ndarray`, optional
        The declinations of the objects in degrees
    mag : :class:`numpy.ndarray`, optional
        The magnitudes of the objects
    filter : :class:`numpy.ndarray`, optional
        The filters used to observe the objects (single character each, e.g. "g", "r", "i")
    obs_code : :class:`str` or :class:`numpy.ndarray`, optional
        The observatory code(s) to use for the observations (if a single string is given, it will be used for
        all observations)
    obs_file : :class:`str` or :class:`pandas.DataFrame`, optional
        Either the path to a .h5 file containing a Pandas table of the observations or the table itself as a 
        Pandas DataFrame. If given, the `obj_id`, `time`, `ra`, `dec`, `mag`, and `filter` parameters will be
        ignored and the data will be read from the file or DataFrame.
    columns : :class:`dict`, optional
        A dictionary containing the column names used to read the data from the `obs_file`. The keys should be
        the same as the parameter names and the values should be the column names in the file
    out_path : :class:`str`, optional
        The path to the file to write the formatted observations to

    Notes
    -----
    The 80 column format is described in detail on the Minor Planet Center's website:
    https://www.minorplanetcenter.net/iau/info/OpticalObs.html

    Raises
    ------
    FileNotFoundError
        If the `obs_file` parameter is given and the file does not exist
    """
    # check if obs_file is given and read in the data if it is
    if obs_file is not None:
        # check if obs_file is a string or a DataFrame
        if not isinstance(obs_file, (pd.DataFrame, str)):
            raise TypeError("obs_file must be a string or a Pandas DataFrame")
        
        # if obs_file is a DataFrame, use it directly
        if isinstance(obs_file, pd.DataFrame):
            df = obs_file

        # if obs_file is a string, read in the data from the file
        else:
            if not isfile(obs_file):
                raise FileNotFoundError(f"File {obs_file} does not exist")
            df = pd.read_hdf(obs_file, key="df")

        # get the data from the DataFrame
        obj_id, time, ra, dec, mag, filter, obs_code = df[columns["id"]].values, df[columns["time"]].values,\
            df[columns["ra"]].values, df[columns["dec"]].values, df[columns["mag"]].values,\
            df[columns["filter"]].values, df[columns["obs_code"]].values

    # check if all the parameters have the same length    
    n_rows = len(obj_id)
    if not all(len(x) == n_rows for x in [time, ra, dec, mag, filter]):
        raise ValueError("All parameters must have the same length")
    
    # if obs_code is a single string, repeat it over a list of the same length as the other parameters
    if isinstance(obs_code, str):
        obs_code = [obs_code] * n_rows

    # convert RA and Dec to hourangles and MJD to regular dates
    ra_degrees = Angle(ra, unit="deg").hms
    dec_degrees = Angle(dec, unit="deg").dms
    t = Time(time, format="mjd").datetime

    # each line stars with 5 spaces
    lines = [" " * 5 for _ in range(n_rows)]
    for i in range(n_rows):
        # save the ID
        lines[i] += obj_id[i]

        # add two spaces and a C (for CCD, I think?)
        lines[i] += " " * 2 + "C"

        # convert time to HH MM DD.ddddd format
        lines[i] += f"{t[i].year:4.0f} {t[i].month:02.0f} {t[i].day + time[i] % 1.0:08.5f} "

        # convert RA to HH MM SS.ddd
        lines[i] += f"{ra_degrees.h[i]:02.0f} {ra_degrees.m[i]:02.0f} {ra_degrees.s[i]:06.3f}"

        # convert Dec to sDD MM SS.dd
        lines[i] += f"{dec_degrees.d[i]:+03.0f} {abs(dec_degrees.m[i]):02.0f} {abs(dec_degrees.s[i]):05.2f}"

        # leave 9 blank columns
        lines[i] += " " * 9

        # add the magnitude and filter (right aligned)
        lines[i] += f"{mag[i]:04.1f}  {filter[i]}"

        # add 5 some more spaces and an observatory code
        lines[i] += " " * 5 + obs_code[i] + "\n"

    if out_path is not None:
        # write that to a file
        with open(out_path, "w") as obs_file:
            obs_file.writelines(lines)
    else:
        return "".join(lines)


def main():
    # setup argparse to run format_obs80 from the command line
    parser = argparse.ArgumentParser(description="Format a set of observations into the 80 column format used by the Minor Planet Center")
    parser.add_argument("--obj_id", type=str, nargs="+", help="The unique IDs (7 character strings) of the objects being observed")
    parser.add_argument("--time", type=float, nargs="+", help="The observation times in MJD")
    parser.add_argument("--ra", type=float, nargs="+", help="The right ascensions of the objects in degrees")
    parser.add_argument("--dec", type=float, nargs="+", help="The declinations of the objects in degrees")
    parser.add_argument("--mag", type=float, nargs="+", help="The magnitudes of the objects")
    parser.add_argument("--filter", type=str, nargs="+", help="The filters used to observe the objects")
    parser.add_argument("--obs_code", type=str, help="The observatory code(s) to use for the observations")
    parser.add_argument("--obs_file", type=str, help="The path to a .h5 file containing a Pandas table of the observations")
    parser.add_argument("--out_path", type=str, help="The path to the file to write the formatted observations to")
    args = parser.parse_args()

    output = format_obs80(obj_id=args.obj_id, time=args.time, ra=args.ra, dec=args.dec, mag=args.mag,
                          filter=args.filter, obs_code=args.obs_code, obs_file=args.obs_file,
                          out_path=args.out_path)
    if args.out_path is None:
        print(output)

if __name__ == "__main__":
    main()