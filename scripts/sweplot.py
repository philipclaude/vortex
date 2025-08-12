'''
Program to create images from vortex .json particle output.
'''
import json
import math
import argparse
import matplotlib.pyplot as plt

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

def main(name, src, out):
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
      s = plt.scatter(l, t, c=h, cmap='coolwarm', s=5, edgecolors='none')
      cbar = plt.colorbar(s, orientation='vertical', location='right', fraction=0.05, shrink=0.675)
      cbar.ax.set_title('[m]', pad=10)
      color_values = [round(min(h))] + SETUP[name]['height_labels'] + [round(max(h))]
      cbar.set_ticks(color_values)

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
  args = parser.parse_args()
  assert args.name and args.src and args.out
  main(args.name, args.src, args.out)
