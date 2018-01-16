import os
import shutil
import numpy as np
import copy
from multiprocessing import Pool


class model(object):
    def __init__(self, ori_elenod='ori/elenod.dat'):
        super(model, self).__init__()
        self.ori_elenod = ori_elenod
        self.model_info()

    def model_info(self):
        with open(self.ori_elenod, 'r') as f:
            data = f.readlines()

        ele_num = int(data[0])
        self.ele_num = ele_num
        self.mesh_data = np.loadtxt(data[1:ele_num + 1], dtype=int)
        self.other_data = data[ele_num + 1:]

        mats = list(self.mesh_data.T[9])
        self.mat_num = max(mats)
        print('there are {} kind of materials in the model'.format(
            self.mat_num))

        self.mats_sn = list(set(mats))
        self.mat_eles = []  # how many elements does one mat have
        for mat in self.mats_sn:
            m_e = mats.count(mat)
            self.mat_eles.append(m_e)
        mats_sn = np.array(self.mats_sn, dtype=int)
        mat_eles = np.array(self.mat_eles, dtype=int)
        self.mat_matrix = np.vstack([mats_sn, mat_eles]).T
        print('material seriel number and elements they have')
        print(self.mat_matrix)

    def random_model(self, rand_elenod='rand_elenod.dat'):
        mat_matrix = self.mat_matrix
        jack = []  # stop number
        jackie = 1
        for i in mat_matrix:
            jackie += i[1]
            jack.append(jackie)
        mat_matrix = np.vstack([mat_matrix.T, jack]).T
        joe = mat_matrix.T[2] - mat_matrix.T[1]  # start number
        mat_matrix = np.vstack([mat_matrix.T, joe]).T.tolist()
        self.mat_matrix_all = copy.deepcopy(mat_matrix)
        new_model = []
        for elements in self.mesh_data.tolist():
            old_mat = elements[9]
            new_mat = mat_matrix[old_mat - 1][-1]
            mat_matrix[old_mat - 1][-1] += 1
            new_ele = elements
            new_ele[9] = new_mat
            new_model.append(new_ele)
        # print(new_model)
        self.new_model = new_model
        rand_model_str = ''.join([''.join(
            [str(ii) + ' 'for ii in i]) + '\n' for i in self.new_model])
        with open(rand_elenod, 'w') as f:
            f.write('{}\n'.format(self.ele_num))
            f.write(rand_model_str)
            f.write(''.join(self.other_data))


class mat(object):
    '''make a sdas control file first
    then split into 3 files
    load_control.txt|mat_control.txt|mat.txt
    and additionally a text file to configure the SD
    for every real material is needed'''

    def __init__(self, random_seed_number, mat_matrix,
                 random_config='random.config',
                 ori_load_control='ori/load_control.txt',
                 ori_mat_control='ori/mat_control.txt',
                 ori_mat='ori/mat.txt'):
        super(mat, self).__init__()
        self.mat_matrix = mat_matrix
        self.random_config = np.loadtxt(random_config)
        self.random_seed_number = random_seed_number
        self.ori_load_control = open(ori_load_control).read()
        self.ori_mat_control = np.loadtxt(ori_mat_control, dtype=int)
        self.ori_mat = np.loadtxt(ori_mat)

    def random_mat(self):
        # give every mat a mat control
        mat_matrix = self.mat_matrix

        new_mat_control = []
        for i in self.ori_mat_control.tolist():
            start = mat_matrix[i[0] - 1][-1]
            stop = mat_matrix[i[0] - 1][-2]
            jack = [[ii] + i[1:] for ii in range(start, stop)]
            new_mat_control.extend(jack)
        new_mat_control_str = ''.join([''.join(
            [str(ii) + ' 'for ii in i]) + '\n' for i in new_mat_control])

        # store the material serial numbers
        sn = np.array(new_mat_control, dtype=int).T[0]

        # mess up the materials :-D
        random_seeds = [i for i in range(self.random_seed_number)]
        for seed in random_seeds:
            np.random.seed(seed)
            rand_mats = []
            for mat in self.ori_mat:
                mat_num = int(mat[0])
                # print(mat_num)
                ori_mat = mat[1:]
                rand_para = np.random.normal(
                    1, self.random_config[mat_num - 1],
                    mat_matrix[mat_num - 1][1])
                rand_mat = [list(i * ori_mat) for i in rand_para]
                rand_mats.extend(rand_mat)
            # print(rand_mats)
            rand_mats = np.vstack([sn, np.array(rand_mats).T]).T.tolist()
            for i in range(len(rand_mats)):
                rand_mats[i][0] = int(rand_mats[i][0])
            rand_mats_str = ''.join([''.join(
                [str(round(ii, 2)) + ' 'for ii in i]) +
                '\n' for i in rand_mats])
            with open('rand_{}.txt'.format(seed), 'w') as f:
                f.write(self.ori_load_control)
                f.write(new_mat_control_str)
                f.write(rand_mats_str)


def run_rand(seeds):
    seeds_num = seeds
    test_model = model()
    test_model.random_model()
    test_mat_matrix = test_model.mat_matrix_all
    test_mat = mat(seeds_num, test_mat_matrix)
    test_mat.random_mat()

    if os.path.exists('rands'):
        shutil.rmtree('rands')
        pass
    else:
        os.mkdir('rands')

    for i in range(0, seeds_num):
        path = 'rands/{}'.format(i)
        os.makedirs(path)
        shutil.copy('rand_elenod.dat', '{}/elenod.dat'.format(path))
        shutil.copy('rand_{}.txt'.format(i), '{}/mat.txt'.format(path))

    [os.remove('rand_{}.txt'.format(i)) for i in range(0, seeds_num)]
    os.remove('rand_elenod.dat')


def run_sdas(work_path):
    os.chdir(work_path)
    os.system('sdas-no-read.exe')
    os.chdir('..')


if __name__ == '__main__':
    run_rand(100)
    os.chdir('rands')
    rands_list = os.listdir()
    work_path_list = [i for i in rands_list]
    print(work_path_list)
    pool = Pool(8)
    pool.map(run_sdas, work_path_list)
