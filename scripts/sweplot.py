'''
Program to create images from vortex .json particle output.
'''
import json
import math
import os
import argparse
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Delaunay # pylint: disable=no-name-in-module
from scipy.spatial import KDTree
# pylint: disable=no-name-in-module
from netCDF4 import Dataset

UNITS = {
  'h': 'm',
  'd': 'm',
  'rv': '1/s'
}

# for Delaunay implementation
def pad_corners(lon: np.ndarray, lat: np.ndarray, quantity: np.ndarray):
  '''
  Extends the longitude and latitude arrays to include the corners
  '''
  corners = np.array([[-math.pi, -math.pi/2], [math.pi, -math.pi/2],
                      [-math.pi, math.pi/2], [math.pi, math.pi/2]], dtype=float)
  corner_values = []
  for lonc, latc in corners:
    d2 = (lon - lonc)**2 + (lat - latc)**2
    corner_values.append(quantity[np.argmin(d2)])
  corner_values = np.asarray(corner_values, float)
  lon_ext = np.concatenate([lon, corners[:, 0]])
  lat_ext = np.concatenate([lat, corners[:, 1]])
  q_ext = np.concatenate([quantity, corner_values])
  return lon_ext, lat_ext, q_ext

def get_quantity(data, field, qref, tree, a):
  """
  Retrieves the quantity to plot from data.
  """
  q = np.array(data['h'])
  assert q is not None

  if 'hs' in data and field == 'h':
    q += np.array(data['hs'])
  elif 'rv' in data:
    q = np.array(data['rv'])
  if tree:
    x = data['x']
    y = data['y']
    z = data['z']
    if 'hs' in data:
      q -= np.array(data['hs']) # compare depth saved by swe-python
    # pylint: disable=consider-using-enumerate
    for i in range(len(q)):
      info = tree.query([a * x[i], a * y[i], a * z[i]])
      q[i] = (q[i] - qref[info[1]]) * 100 / qref[info[1]]
  return q

def main(days, plot_type, src, out, field, diff):
  '''
  Runs the main plotting program.
  '''
  # import reference solution at the last day if provided
  ref = None
  tree = None
  a = None
  if diff:
    assert field == 'h'
    ref = Dataset(diff, "r", format="NETCDF4")
    a = ref.sphere_radius
    xc = np.array(ref["xCell"])
    yc = np.array(ref["yCell"])
    zc = np.array(ref["zCell"])
    n_cells = len(xc)
    assert len(yc) == n_cells and len(zc) == n_cells

    # build a kdtree from the stationary points (cell centers)
    points = np.zeros([n_cells, 3])
    points[:, 0] = xc
    points[:, 1] = yc
    points[:, 2] = zc
    tree = KDTree(points)

  minval = float('inf')
  maxval = -minval
  if len(days) == 1:
    days = range(days[0] + 1)
  for day in days:
    print(f"Processing day {day}")
    href = None
    if ref:
      assert field == 'h'
      href = ref.variables["hh_cell"][day, :, :].squeeze()
    assert os.path.exists(f"{src}/particles{24 * day}.json")
    with open(f"{src}/particles{24 * day}.json", encoding='utf-8') as f:
      plt.figure()

      # import the particle locations and height
      data = json.loads(f.read())
      q = get_quantity(data, field, href, tree, a)
      qmin = min(q)
      qmax = max(q)
      minval = min(qmin, minval)
      maxval = max(qmax, maxval)
      print(f"Day {day}: qmin = {qmin}, qmax = {qmax}")
  print(f"min = {minval}, max = {maxval}")

  for day in days:
    print(f"Processing day {day}")
    href = None
    if ref:
      assert field == 'h'
      href = ref.variables["hh_cell"][day, :, :].squeeze()
    with open(f"{src}/particles{24 * day}.json", encoding='utf-8') as f:
      plt.figure()

      # import the particle locations and height
      data = json.loads(f.read())
      x = np.array(data['x'])
      y = np.array(data['y'])
      z = np.array(data['z'])
      q = get_quantity(data, field, href, tree, a)
      t = np.asin(z) # latitude (theta)
      l = np.atan2(y, x) # longitude (lambda)

      # plot
      if plot_type == 'point':
        s = plt.scatter(l, t, c=q, cmap='coolwarm', s=5, edgecolors='none')
        cbar = plt.colorbar(s, orientation='vertical', location='right',
                            fraction=0.05, shrink=0.675)
        if diff:
          cbar.ax.set_title("%", pad=10)
        else:
          cbar.ax.set_title(f"[{UNITS[field]}]", pad=10)
      elif plot_type == 'tri':
        lon_ext, lat_ext, height_ext = pad_corners(l, t, q)
        tri = Delaunay(np.column_stack([lon_ext, lat_ext]))
        s = plt.tripcolor(
          lon_ext, lat_ext, tri.simplices, height_ext, vmin=minval, vmax=maxval,
          cmap='coolwarm', edgecolors='none'
        )
        cbar = plt.colorbar(s, orientation='vertical', location='right',
                            fraction=0.05, shrink=0.675)
        if diff:
          cbar.ax.set_title("%", pad=10)
        else:
          cbar.ax.set_title(f"[{UNITS[field]}]", pad=10)
      else:
        raise TypeError(f"unknown plot type {plot_type}")

      plt.xlabel('$\\lambda$')
      plt.ylabel('$\\theta$', rotation=0)
      plt.xticks(ticks=[-math.pi, -math.pi/2, 0, math.pi/2, math.pi],
                 labels=['$-180^\\circ$', '$-90^\\circ$', '0', '$90^\\circ$', '$180^\\circ$'])
      plt.yticks(ticks=[-math.pi/2, -math.pi/4, 0, math.pi/4, math.pi/2],
                 labels=['$-90^\\circ$', '$-45^\\circ$', '0', '$45^\\circ$', '$90^\\circ$'])

      ax = plt.gca()
      ax.tick_params(axis='x', which='both', length=0)
      ax.tick_params(axis='y', which='both', length=0)

      ax.spines[:].set_visible(False)
      ax.spines['left'].set_position(('outward', -10))

      plt.tight_layout()
      ax.set_aspect('equal')
      plt.savefig(f"{out}-{field}-day{day}.png", dpi=200, bbox_inches='tight', pad_inches=0.0)
      plt.close()


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('--src', help='directory containing particles*.json files')
  parser.add_argument('--out', help='output prefix for images')
  parser.add_argument('--type', help='either point or tri', default='tri')
  parser.add_argument('--field', help='which quantity to plot (h, rv)', default='h')
  parser.add_argument('--days', type=int, nargs='+', help='how many days to plot', default=1)
  parser.add_argument('--diff', default='',
                      help="plots the difference in the field, given a reference solution")
  args = parser.parse_args()
  assert args.src and args.out
  main(args.days, args.type, args.src, args.out, args.field, args.diff)
