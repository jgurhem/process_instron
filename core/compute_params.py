
def find_point(start, end, lissage, depla, my_lambda):
  y = 0
  x = 0
  for i in range(start, end):
    if my_lambda(i):
      y = lissage[i]
      x = depla[i]
      break
  results = dict()
  results["y"] = y
  results["x"] = x
  results["index"] = i
  return results
