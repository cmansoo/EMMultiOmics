import pandas as pd
import sys
from mpmath import pcfd, log


mypath = "C:/Users/mcho1/Desktop/umn/R/EMMultiOmics/"
print("mypath: ", mypath)
a0 = 0.1
print("a0: ", a0)
gstr = 4.16233090530697e-05
print("gstr: ", gstr)

zv_path= "".join(["zv", str(a0), "_", str(gstr), ".txt"])
print("zv_path: ", zv_path)

D = pd.read_csv(zv_path)
print("D: ", D)

print(D.columns.tolist())
i = int(D.columns.tolist()[0])
print("i: ", i)

result = []
# read
with open(file=zv_path, mode="r") as f:
    for line in f:
        result.append(list(map(float, line.split(","))))
        print("line: ", line)
        print("result: ", result)



v = result[0][1:(i+1)]
print("v: ", v)

z = result[0][(i+1):(2*i+1)]
print("z: ", z)

A = i * [0]
print("A: ", A)

for j in range(i):
    qty = float(log(pcfd(v[j], z[j])))
    A[j] = qty
    print("qty: ", qty)

# write
out1_path = "".join(["out1", str(a0), "_", str(gstr), ".txt"])
print("out1_path: ", out1_path)
out1 = open(out1_path, mode="w")
export_values = " ".join(repr(e) for e in A)
print("export_values: ", export_values)

# write export values to the file 
print(out1_path, file=out1)
out1.close()

# import mpmath
# print(mpmath.__file__)