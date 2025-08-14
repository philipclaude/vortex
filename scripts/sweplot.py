'''
Program to create images from vortex .json particle output.
'''
import json
import math
import argparse
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Delaunay # pylint: disable=no-name-in-module

# times are in hours, height_labels are in metres
SETUP = {
  'w2': { 'times': [0, 120, 288], 'height_labels': [1500, 2000, 2500]},
  'w5': { 'times': [0, 120, 240, 360], 'height_labels': [5250, 5500, 5750]},
  'w6': { 'times': [0, 120, 240, 360], 'height_labels': [8500, 9000, 9500, 10000]},
}

def hs5(l, t):
  '''
  Calculates the surface height for Williamson Case 5 because the vortex
   .json files only save the depth.
  '''
  hs0 = 2000
  r = math.pi / 9.0
  lambda_c = -math.pi / 2
  theta_c = math.pi / 6
  d_lambda = l - lambda_c
  d_theta = t - theta_c
  d = min(r * r, d_lambda * d_lambda + d_theta * d_theta)
  hs = hs0 * (1 - d ** 0.5 / r)
  return hs

lon_min, lon_max = -math.pi, math.pi
lat_min, lat_max = -math.pi/2, math.pi/2
# for Delaunay implementation
def pad_corners(lon, lat, h):
  lon_arr = np.asarray(lon, float)
  lat_arr = np.asarray(lat, float)
  height_arr = np.asarray(h, float)
  corners = np.array([[lon_min, lat_min], [lon_max, lat_min],
                      [lon_min, lat_max], [lon_max, lat_max]], dtype=float)
  corner_heights = []
  for lonc, latc in corners:
    d2 = (lon_arr - lonc)**2 + (lat_arr - latc)**2
    corner_heights.append(height_arr[np.argmin(d2)])
  corner_heights = np.asarray(corner_heights, float)
  lon_ext = np.concatenate([lon_arr, corners[:, 0]])
  lat_ext = np.concatenate([lat_arr, corners[:, 1]])
  h_ext = np.concatenate([height_arr, corner_heights])
  return lon_ext, lat_ext, h_ext

def main(name, plot_type, src, out):
  '''
  Runs the main plotting program.
  '''
  for step in SETUP[name]['times']:
    with open(f"{src}/particles{step}.json", encoding='utf-8') as f:
      plt.figure()

      # import the particle locations and height
      data = json.loads(f.read())
      x = data['x']
      y = data['y']
      z = data['z']
      d = data['h']

      # compute latitude and longitude
      t = [math.asin(value) for value in z]
      l = [math.atan2(y[i], x[i]) for i in range(len(x))]

      # determine height to plot
      if name == 'w5':
        h = [d[i] + hs5(l[i], t[i]) for i in range(len(x))]
      else:
        h = d

      # plot
      if plot_type == 'point':
        s = plt.scatter(l, t, c=h, cmap='coolwarm', s=5, edgecolors='none')
        cbar = plt.colorbar(s, orientation='vertical', location='right',
                            fraction=0.05, shrink=0.675)
        cbar.ax.set_title('[m]', pad=10)
        color_values = [round(min(h))] + SETUP[name]['height_labels'] + [round(max(h))]
        cbar.set_ticks(color_values)
      elif plot_type == 'tri':
        l_arr = np.asarray(l, float)
        t_arr = np.asarray(t, float)
        h_arr = np.asarray(h, float)
        lon_ext, lat_ext, height_ext = pad_corners(l_arr, t_arr, h_arr)
        tri = Delaunay(np.column_stack([lon_ext, lat_ext]))
        vmin = float(h_arr.min())
        vmax = float(h_arr.max())
        s = plt.tripcolor(
          lon_ext, lat_ext, tri.simplices, height_ext,
          cmap='coolwarm', vmin=vmin, vmax=vmax, edgecolors='none'
        )
        cbar = plt.colorbar(s, orientation='vertical', location='right',
                            fraction=0.05, shrink=0.675)
        cbar.ax.set_title('[m]', pad=10)
        cbar.set_ticks([round(vmin)] + SETUP[name]['height_labels'] + [round(vmax)])

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
      plt.savefig(f"{out}/{name}-t{step}.png", dpi=200, bbox_inches='tight', pad_inches=0.0)
      plt.close()


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('--name', help='case to plot: either w2, w5 or w6')
  parser.add_argument('--src', help='directory containing particles*.json files')
  parser.add_argument('--out', help='output directory for images')
  parser.add_argument('--type', help='either point or tri', default='point')
  args = parser.parse_args()
  assert args.name and args.src and args.out
  main(args.name, args.type, args.src, args.out)
