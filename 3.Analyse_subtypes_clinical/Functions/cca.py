# -*- coding: utf-8 -*-
import hdf5storage
import os
import scipy.io as sio
import numpy as np
import time
from sklearn import linear_model
from sklearn import preprocessing
from sklearn import cross_decomposition
from joblib import Parallel, delayed
import statsmodels.formula.api as sm
from sklearn.cross_decomposition import CCA
from sklearn.cross_decomposition import PLSRegression
from sklearn.cross_decomposition import PLSCanonical
from scipy.stats import spearmanr
import pandas as pd
import numpy as np
from sklearn.linear_model import Lasso
from sklearn.metrics import r2_score
import math
# need to combine with caculate_corr.m (/home/cxpang/matlab/code/4.HAMD_correlation/Subscale)
CodesPath = '/HeLabData2/cxpang/DIDA/NGSR/hamd/pncsinglefuncparcel_psychopathology-main/Functions';


def PLSr1_KFold_RandomCV_MultiTimes(select_p,Subjects_Data, Subjects_Score, Covariates, Fold_Quantity, ComponentNumber_Range,
                                    CVRepeatTimes, ResultantFolder, Parallel_Quantity, Permutation_Flag, Queue,
                                    RandIndex_File_List=''):
    if not os.path.exists(ResultantFolder):
        os.makedirs(ResultantFolder);
    # only one time
    for i in np.arange(CVRepeatTimes):
        ResultantFolder_TimeI = ResultantFolder + '/Time_' + str(i)
        if not os.path.exists(ResultantFolder_TimeI):
            os.makedirs(ResultantFolder_TimeI);
        if RandIndex_File_List != '':
            RandIndex_File = RandIndex_File_List[i]
        else:
            RandIndex_File = '';
        PLSr1_KFold_RandomCV(select_p,Subjects_Data, Subjects_Score, Covariates, Fold_Quantity,
                                            ComponentNumber_Range, 1, ResultantFolder_TimeI, Parallel_Quantity,
                                            Permutation_Flag, RandIndex_File)






def PLSr1_KFold_RandomCV(select_p,Subjects_Data, Subjects_Score, Covariates, Fold_Quantity, ComponentNumber_Range,
                         CVRepeatTimes_ForInner, ResultantFolder, Parallel_Quantity, Permutation_Flag,
                         RandIndex_File=''):
    
   
    Subjects_Quantity = len(Subjects_Score)
    # not random, using sorting to divide training and test
    Features_Quantity = np.shape(Subjects_Data)[1];
    for j in np.arange(Fold_Quantity):# ten-cv
        print('fold' + str(j))
        tmp=np.sum(Subjects_Score,axis=1)# according the hamd subscale summed score
        SortedID = sorted(range(len(tmp)), key=lambda k: tmp[k])
        Fold_J_Index = SortedID[j: Subjects_Quantity : Fold_Quantity] ;
        Subjects_Data_test = Subjects_Data[Fold_J_Index, :]
        Subjects_Score_test = Subjects_Score[Fold_J_Index,:]


        Covariates_test = Covariates[Fold_J_Index, :]
        Subjects_Data_train = np.delete(Subjects_Data, Fold_J_Index, axis=0)
        Subjects_Score_train = np.delete(Subjects_Score, Fold_J_Index, axis=0)
        Covariates_train = np.delete(Covariates, Fold_J_Index, axis=0)
        Covariates_Quantity = np.shape(Covariates)[1]
        Fold_J_result = {'Subjects_Data_test': Subjects_Data_test, 'Subjects_Data_train': Subjects_Data_train}
        Fold_J_FileName = 'Fold_' + str(j) + '_remove65.mat'
        maskFile = os.path.join('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/sub1', Fold_J_FileName)
        hdf5storage.savemat(maskFile, Fold_J_result)

        Fold_J_result = {'Subjects_Score_test': Subjects_Score_test,'Subjects_Score_train':Subjects_Score_train}
        Fold_J_FileName = 'Fold_' + str(j) + '_score_remove65.mat'
        maskFile = os.path.join('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/sub1', Fold_J_FileName)
        hdf5storage.savemat(maskFile, Fold_J_result)
    
        # Controlling covariates from brain data
        df = {};
        for k in np.arange(Covariates_Quantity):
            df['Covariate_' + str(k)] = Covariates_train[:, k];
        # Construct formula
        Formula = 'Data ~ Covariate_0';
        for k in np.arange(Covariates_Quantity - 1) + 1:
            Formula = Formula + ' + Covariate_' + str(k)
        # Regress covariates from each brain features
        a=0
        if (j <= 9 and a==0):#step2
            print('loading regresss')
            Fold_J_FileName = 'Fold_' + str(j) + '_regress.mat'
            maskFile = os.path.join('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/sub1', Fold_J_FileName)
            tmp_mat=hdf5storage.loadmat(maskFile)
            Subjects_Data_train= tmp_mat['Subjects_Data_train'];
            Subjects_Data_test= tmp_mat['Subjects_Data_test'];
            # for score
            Covariates_score=Covariates[:, :2]
            Covariates_Quantity_s = np.shape(Covariates_score)[1]
            Covariates_train_score=Covariates_train[:, :2]
            Covariates_test_score=Covariates_test[:, :2]
            df = {};
            for k in np.arange(Covariates_Quantity_s):
                df['Covariate_' + str(k)] = Covariates_train_score[:, k];
        # Construct formula
            Formula = 'Data ~ Covariate_0';
            for k in np.arange(Covariates_Quantity_s - 1) + 1:
                Formula = Formula + ' + Covariate_' + str(k)
            Features_Quantity = np.shape(Subjects_Score)[1];
            for k in np.arange(Features_Quantity):
           # print(k)
                df['Data'] = Subjects_Score_train[:, k];
                # Regressing covariates using training data
                LinModel_Res = sm.ols(formula=Formula, data=df).fit()
                # Using residuals replace the training data
                Subjects_Score_train[:, k] = LinModel_Res.resid;

                # Calculating the residuals of testing data by applying the coeffcients of training data
                Coefficients = LinModel_Res.params;
                Subjects_Score_test[:, k] = Subjects_Score_test[:, k] - Coefficients[0];
                for m in np.arange(Covariates_Quantity_s):
                    Subjects_Score_test[:, k] = Subjects_Score_test[:, k] - Coefficients[m + 1] * Covariates_test_score[:, m]
            Fold_J_result = {'Subjects_Score_train': Subjects_Score_train}
            Fold_J_FileName = 'Fold_' + str(j) + '_regress_score.mat'
            maskFile = os.path.join('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/sub2', Fold_J_FileName)
            hdf5storage.savemat(maskFile, Fold_J_result)
        else:#sstep1 regress, saving time
            for k in np.arange(Features_Quantity):
                df['Data'] = Subjects_Data_train[:, k];
                # Regressing covariates using training data
                LinModel_Res = sm.ols(formula=Formula, data=df).fit()
                # Using residuals replace the training data
                Subjects_Data_train[:, k] = LinModel_Res.resid;
                # Calculating the residuals of testing data by applying the coeffcients of training data
                Coefficients = LinModel_Res.params;
                Subjects_Data_test[:, k] = Subjects_Data_test[:, k] - Coefficients[0];
                for m in np.arange(Covariates_Quantity):
                    Subjects_Data_test[:, k] = Subjects_Data_test[:, k] - Coefficients[m + 1] * Covariates_test[:, m]

            Fold_J_result = {'Subjects_Data_test': Subjects_Data_test,'Subjects_Data_train':Subjects_Data_train}
            Fold_J_FileName = 'Fold_' + str(j) + '_regress.mat'
            maskFile = os.path.join('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/sub1', Fold_J_FileName)
            hdf5storage.savemat(maskFile, Fold_J_result)

        if j <= 9 and select_p<1:#step4
            print('loading mask')
            Fold_J_FileName = 'fold_' + str(j) + '_mask.mat'
            maskFile = os.path.join('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/sub2/p'+str(select_p), Fold_J_FileName)
            tmp=hdf5storage.loadmat(maskFile)
            mask=tmp['pos']
            mask=mask.tolist()
            mask = [i for j in mask for i in j] ;
            b = []
            for k in np.arange(Features_Quantity):
                    b.append(k)
            my_mask = list(set(b) - set(mask))
        else:
            my_mask=[]

        if Permutation_Flag:
            # If do permutation, the training scores should be permuted, while the testing scores remain
            Subjects_Index_Random = np.arange(len(Subjects_Score_train))
            np.random.shuffle(Subjects_Index_Random)
            Subjects_Score_train = Subjects_Score_train[Subjects_Index_Random,:]
            if j == 0:
                PermutationIndex = {'Fold_0': Subjects_Index_Random}
            else:
                PermutationIndex['Fold_' + str(j)] = Subjects_Index_Random

        Subjects_Data_train_select=np.delete(Subjects_Data_train,my_mask,axis=1)
        Subjects_Data_test_select=np.delete(Subjects_Data_test,my_mask,axis=1)


        normalize = preprocessing.MinMaxScaler()
        Subjects_Data_train_select = normalize.fit_transform(Subjects_Data_train_select)
        Subjects_Data_test_select = normalize.transform(Subjects_Data_test_select)
        
        Optimal_ComponentNumber= ComponentNumber_Range#PLSr1_OptimalComponentNumber_KFold(Subjects_Data_train, Subjects_Score_train, Fold_Quantity, ComponentNumber_Range, 1, ResultantFolder, 1)
        #clf = CCA(n_components=Optimal_ComponentNumber)
        clf = PLSCanonical(n_components=Optimal_ComponentNumber)
        #clf = PLSRegression(n_components=Optimal_ComponentNumber)
        
        clf.fit(Subjects_Data_train_select, Subjects_Score_train)
        X_test_r, Y_test_r = clf.transform(Subjects_Data_test_select, Subjects_Score_test)
    
        #for num in range (Optimal_ComponentNumber):
        Fold_J_result = {'u': X_test_r, 'v': Y_test_r}
        Corr = np.corrcoef(X_test_r[:,0], Y_test_r[:,0])[0, 1]
        print(Corr)

        Fold_J_FileName = 'Fold_' + str(j) +'comp_uv.mat'
        ResultantFile = os.path.join(ResultantFolder, Fold_J_FileName)
        
        hdf5storage.savemat(ResultantFile, Fold_J_result)


