'''
Program to compare vortex solutions with swe-python.
'''
import argparse
import json
import netCDF4 as nc
from scipy.spatial import KDTree
import numpy as np
from matplotlib import pyplot as plt

def main(src, ref, n_days, out):
  '''
  Compares all vortex solutions (in src) with the reference solution.
  '''
  plt.style.use('tableau-colorblind10')
  days = list(range(n_days + 1))
  for s in src:

    h_error = []
    for day in range(0, n_days + 1):

      hour = day * 24
      vtx_file = f"{s}/particles{hour}.json"
      with open(vtx_file, encoding='utf-8') as f:
        vtx_data = json.loads(f.read())
        xp = vtx_data["x"]
        yp = vtx_data["y"]
        zp = vtx_data["z"]
        hp = vtx_data["h"]

      data = nc.Dataset(ref, "r", format="NETCDF4")
      a = data.sphere_radius

      h = data.variables["hh_cell"][day, :, :]
      xc = np.array(data["xCell"])
      yc = np.array(data["yCell"])
      zc = np.array(data["zCell"])
      n_cells = len(xc)
      assert len(yc) == n_cells and len(zc) == n_cells

      # build a kdtree from the stationary points (cell centers)
      points = np.zeros([n_cells, 3])
      points[:, 0] = xc
      points[:, 1] = yc
      points[:, 2] = zc
      tree = KDTree(points)

      error = 0
      n_particles = len(xp)
      sum_h = 0
      for i in range(n_particles):
        # get closest point and height value
        info = tree.query([a * xp[i], a * yp[i], a * zp[i]])
        ha = h[info[1]][0]
        error += (ha - hp[i]) ** 2
        sum_h += ha ** 2
      eh = (error / sum_h) ** 0.5
      print(f"Day {day}: error = {eh}")
      h_error.append(eh)

    plt.semilogy(days, h_error, '-o', label=f'{n_particles} particles')

  plt.xlabel('day', size=14)
  plt.ylabel(r'$\Delta I_h$', rotation=0, loc='top', labelpad=-20, size=14)
  plt.xlim([0, n_days])
  days_label = list(range(0, n_days + 1, n_days // 3))
  plt.xticks(ticks=days_label, labels=days_label, size='large')
  plt.yticks(size='large')
  plt.grid('major', axis='y')
  plt.legend(frameon=False, prop={'size': 14}, loc='lower right')
  plt.savefig(out)
  plt.show()

if __name__ == "__main__":
  # Example run using vortex solutions for Williamson test case 5 (w5) with a time step of 30 seconds (t30)
  # using 5, 6 and 7 icosahedron subdivisions to initialize the particles (i5, i6, i7):
  # python3.11 swediff.py --ref out_wtc5-cvt8.nc --src w5-i5-t30/ w5-i6-t30/ w5-i7-t30/ --out w5-difference.pdf --days 15
  parser = argparse.ArgumentParser()
  parser.add_argument('--src', type=str, nargs='+',
                      help='list of directories containing particles*.json files')
  parser.add_argument('--ref', help='reference NetCDF file from swe-python')
  parser.add_argument('--days', type=int, help='# days to plot')
  parser.add_argument('--out', help='output file with comparison')
  args = parser.parse_args()
  assert args.src and args.ref and args.out
  main(args.src, args.ref, args.days, args.out)
