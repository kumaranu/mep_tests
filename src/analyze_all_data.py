import numpy as np

f = open('all_data.txt', 'r')
INF = f.readlines()
for line in INF:
    # print(line.split())
    tmp_data_list = line.split()
    index = tmp_data_list[0].strip(':')
    # print(index)
    ts_e_err = float(tmp_data_list[4])
    if abs(ts_e_err) > 0.15:
        print(index, ts_e_err)
f.close()

