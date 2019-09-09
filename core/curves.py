import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os
from .correlation import compute_correlation


def make_curve(data_curves, directory):
  fig = plt.figure()
  ax = fig.gca()
  for i in data_curves:
    s_ = min(len(data_curves[i][0]), len(data_curves[i][1]))
    ax.plot(data_curves[i][0][0:s_], data_curves[i][1][0:s_], label=i)
  #plt.legend()
  fig_name = directory + "/fig.pdf"
  if os.path.isfile(fig_name):
   os.remove(fig_name)
  plt.savefig(fig_name)
  plt.close()

  data_curves_keys = list(data_curves.keys())
  if len(data_curves_keys) < 2: return None
  param_keys = list(data_curves[data_curves_keys[0]][2].keys())
  units = data_curves[data_curves_keys[0]][3]
  s_ = len(param_keys)
  for k1 in range(s_):
    for k2 in range(s_):
      v1 = []
      v2 = []
      vr1 = []
      vr2 = []
      vb1 = []
      vb2 = []
      for i in data_curves:
        v1.append(data_curves[i][2][param_keys[k1]])
        v2.append(data_curves[i][2][param_keys[k2]])
        if "%" in i:
          vr1.append(data_curves[i][2][param_keys[k1]])
          vr2.append(data_curves[i][2][param_keys[k2]])
        else:
          vb1.append(data_curves[i][2][param_keys[k1]])
          vb2.append(data_curves[i][2][param_keys[k2]])
      corr = compute_correlation(v1, v2)
      print(param_keys[k1], " - ", param_keys[k2], corr)
      fig = plt.figure()
      ax = fig.gca()
      title = param_keys[k1] + " - " + param_keys[k2]
      ax.scatter(vr1, vr2, color="red")
      ax.scatter(vb1, vb2, color="blue")
      plt.xlabel(param_keys[k1] + " (" + units[param_keys[k1]] + ")")
      plt.ylabel(param_keys[k2] + " (" + units[param_keys[k2]] + ")")
      plt.title(title + "(correlation = " + str(corr) + ")")
      fig_name = directory + "/" + title + ".pdf"
      if os.path.isfile(fig_name):
        os.remove(fig_name)
      plt.savefig(fig_name)
      plt.close()
