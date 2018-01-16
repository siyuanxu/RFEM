import numpy as np
import os
from rfem import model
from multiprocessing import Pool


def post_process(work_path):

    os.chdir(work_path)
    with open('DASAC1.DAT') as disp:
        data = disp.readlines()
        data = [i for i in data if len(i) > len(data[0])]
    with open('displacement.dat', 'w') as disp:
        [disp.write(i[1:]) for i in data]

    disp = np.loadtxt('displacement.dat')

    ele_num = int(max([i[0]for i in disp]))
    load_num = int(len(disp) / ele_num)

    disp_raw = disp.T[1:].T  # no sn number

    disp_1 = disp_raw[ele_num:ele_num * 2]
    disp_2 = disp_raw[ele_num * (load_num - 1):ele_num * (load_num - 0)]

    final_disp = disp_2 - disp_1  # dam axies |down river| settlement

    ele_nod = model('elenod.dat').mesh_data
    ele_sn = ele_nod.T[0]

    # settlement
    ele_settle = []
    for i in ele_nod:
        jackie = np.mean([final_disp[ii - 1][2] for ii in i[1:9]])
        ele_settle.append(jackie)
    max_settlement = -1 * min(ele_settle)
    id = ele_settle.index(-1 * max_settlement)
    ms_ele_id = ele_nod[id][0]

    os.chdir('..')
    return max_settlement, ms_ele_id


if __name__ == '__main__':
    os.chdir('rands')
    rands_list = os.listdir()

    work_path_list = [i for i in rands_list]
    res = []
    for i in work_path_list:
        res.append(post_process(i))
    os.chdir('..')
    res = np.array(res)
    print(res)
    np.savetxt('settlement.txt', res)

    np.savetxt('no_rand.txt', np.array(post_process('no_rand')))
