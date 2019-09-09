import math

def compute_correlation(x, y):
  if len(x) != len(y):
    raise Error
  s_ = len(x)
  mean_x = 0
  mean_y = 0
  for i in range(s_):
    mean_x += x[i]
    mean_y += y[i]
  mean_x /= s_
  mean_y /= s_

  sigma_x = 0
  sigma_y = 0
  for i in range(s_):
    sigma_x += (x[i] - mean_x) * (x[i] - mean_x)
    sigma_y += (y[i] - mean_y) * (y[i] - mean_y)
  sigma_x = math.sqrt(sigma_x / s_)
  sigma_y = math.sqrt(sigma_y / s_)

  sigma_xy = 0
  for i in range(s_):
    sigma_xy += (x[i] - mean_x) * (y[i] - mean_y)
  sigma_xy /= s_
  r = sigma_xy / sigma_x / sigma_y
  return r * r

