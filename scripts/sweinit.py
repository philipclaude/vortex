"""
Program to generate an initial condition for vortex
from swe-python solution.
"""
import argparse
import json
import netCDF4 as nc
import numpy as np

def main(ref_file, out_file):
  """
  Generates an initial condition from swe-python solution in ref_file,
  writing a .json to out_file with x,y,z,h.
  """
  data = nc.Dataset(ref_file, "r", format="NETCDF4")
  a = data.sphere_radius

  hc = data.variables["hh_cell"][0, :, :].squeeze()
  xc = np.array(data["xCell"]) / a
  yc = np.array(data["yCell"]) / a
  zc = np.array(data["zCell"]) / a
  n_cells = len(xc)
  assert len(yc) == n_cells and len(zc) == n_cells

  ic = {
      'h': hc.tolist(),
      'x': xc.tolist(),
      'y': yc.tolist(),
      'z': zc.tolist()
  }
  with open(out_file, 'w', encoding='utf-8') as f:
    f.write(json.dumps(ic))


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('--ref', help='reference NetCDF file from swe-python')
  parser.add_argument('--out', help='output file with comparison')
  args = parser.parse_args()
  main(args.ref, args.out)
