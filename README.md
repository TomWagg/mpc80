# MPC 80-column formatter
This is a simple script for converting a set of observations of minor planets into the Minor Planet Center 80 column format.

## Usage in Python scripts
You can import this function to use in your python scripts directly as follows.
```python
# optionally, if your script is in a different location
import sys
sys.path.append("path/to/mpc80")

import pandas as pd
import mpc80_formatter
formatted_output = mpc80_formatter.format_obs80(obj_id=["ABCDEFG"], time=[60372],
                                                ra=[1.0], dec=[2.0], mag=[21],
                                                filter=["r"], obs_code="I11")
print(formatted_output)

# creating a small dummy file for demo
# you can customise these column names with `columns` parameter
data = {
    "hex_id": ["ABCDEFG"],
    "FieldMJD_TAI": [60372],
    "AstRA(deg)": [1.0],
    "AstDec(deg)": [2.0],
    "observedTrailedSourceMag": [21],
    "optFilter": ["r"],
    "obs_code": ["I11"]
}
df = pd.DataFrame(data)
df.to_hdf("neo.h5", key="df")
del df

# or format a whole file
mpc80_formatter.format_obs80(obs_file="neo.h5", out_path="neo.obs")

# or use an open pandas dataframe
df = pd.read_hdf("neo.h5", key="df")
mpc80_formatter.format_obs80(obs_file=df, out_path="neo.obs")
```

## Command-line usage
You can also use it on the command line. Run ``python mpc80_formatter.py -h`` to get this help message
```
usage: mpc80_formatter.py [-h] [--obj_id OBJ_ID [OBJ_ID ...]] [--time TIME [TIME ...]] [--ra RA [RA ...]] [--dec DEC [DEC ...]] [--mag MAG [MAG ...]] [--filter FILTER [FILTER ...]] [--obs_code OBS_CODE] [--obs_file OBS_FILE] [--out_path OUT_PATH]

Format a set of observations into the 80 column format used by the Minor Planet Center

options:
  -h, --help            show this help message and exit
  --obj_id OBJ_ID [OBJ_ID ...]
                        The unique IDs (7 character strings) of the objects being observed
  --time TIME [TIME ...]
                        The observation times in MJD
  --ra RA [RA ...]      The right ascensions of the objects in degrees
  --dec DEC [DEC ...]   The declinations of the objects in degrees
  --mag MAG [MAG ...]   The magnitudes of the objects
  --filter FILTER [FILTER ...]
                        The filters used to observe the objects
  --obs_code OBS_CODE   The observatory code(s) to use for the observations
  --obs_file OBS_FILE   The path to a .h5 file containing a Pandas table of the observations
  --out_path OUT_PATH   The path to the file to write the formatted observations to
```

Here's an example to convert a single line and print it out on the command line
```bash
python obs80_formatter.py --obj_id ABCDEFG --time 60372 --ra 1.0 --dec 2.0 --mag 21 --filter r --obs_code I11
```
You could imagine piping this directly to `digest2`.

Alternatively, here's one to convert a file full of observations
```bash
python mpc80_formatter.py --obs_file neo.h5 --out_path neo.obs
```
