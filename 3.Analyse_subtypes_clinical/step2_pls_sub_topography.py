# import panda as pd
import os

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import hdf5storage
import scipy.io as sio
from multiprocessing import Pool, cpu_count
import os, time, random
import os.path
import numpy as np
import os
import sys

sys.path.append('/home/cxpang/matlab/code/4.HAMD_correlation/Subscale/Functions');
import cca
import cca_sub2
import numpy
from sklearn.linear_model import Lasso
import pandas as pd

if __name__ == '__main__':
    ResultsFolder = '/HeLabData2/cxpang/DIDA/NGSR/result_xy/hamd';
    # matFolder='/home/cxpang/matlab/code/4.HAMD_correlation/using_mat'
    matFolder = '/home/cxpang/matlab/code/8.brain_age/using_mat'

    AtlasLoading_Mat = hdf5storage.loadmat(ResultsFolder + '/loading_sub_sub1_remove65.mat');
    SubjectsData = AtlasLoading_Mat['loading_sub1'];

    Psychopathology_Mat = hdf5storage.loadmat(matFolder + '/hamd_sub_sub1_remove65.mat');
    ExternalizingCorrtraits = Psychopathology_Mat['hamd'];

    covariate = hdf5storage.loadmat(matFolder + '/covariate_hamd_sub_sub1_remove65.mat');
    Covariates = np.zeros((len(ExternalizingCorrtraits), 3));  # original is 3

    Covariates[:, 0] = covariate['age'].reshape(len(ExternalizingCorrtraits));
    Covariates[:, 1] = covariate['sex'].reshape(len(ExternalizingCorrtraits));
    Covariates[:, 2] = covariate['mfd'].reshape(len(ExternalizingCorrtraits));
    select_p = 0.05

    # #Range of parameters
    # ComponentNumber_Range = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17];
    ComponentNumber_Range = 3;  # here represent the component number, in permutation only 1 component is enough
    FoldQuantity = 10;
    Parallel_Quantity = 10;
    AtlasLoading_Folder = ResultsFolder + '/cca/AtlasLoading';
    CVRepeatTimes = 1
    ResultantFolder = AtlasLoading_Folder + '/loading_sub_sub1_0.05';
    cca.PLSr1_KFold_RandomCV_MultiTimes(select_p, SubjectsData, ExternalizingCorrtraits, Covariates, FoldQuantity,
                                        ComponentNumber_Range, CVRepeatTimes, ResultantFolder, Parallel_Quantity, 0,
                                        'he_queue.q')