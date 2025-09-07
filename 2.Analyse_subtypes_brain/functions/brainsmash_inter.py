import time
import os
import multiprocessing
from multiprocessing import Pool
from brainsmash.workbench.geo import volume
from brainsmash.mapgen.sampled import Sampled
import scipy.io as sio
from multiprocessing import Pool, cpu_count
import os, time, random
import os.path
def brainsmash(output_dir,data_path):
    start = time.time()
    coord_file = "/HeLabData2/cxpang/DIDA/NGSR/brainsmash/coor.txt"  # 要mask一下"
    filenames = volume(coord_file, output_dir)
    for i in  range(1, 19):
        print(i)
        brain_map =data_path +"outputt_pos_sub1_"+str(i)+".txt"  
        gen = Sampled(x=brain_map, D=filenames['D'], index=filenames['index'])
        surrogate_maps = gen(n=100)
        sio.savemat(output_dir+"outputt_pos_sub1_"+str(i)+".mat" , {'surrogate_maps': surrogate_maps})
        
        brain_map =data_path +"outputt_neg_sub1_"+str(i)+".txt"  
        gen = Sampled(x=brain_map, D=filenames['D'], index=filenames['index'])
        surrogate_maps = gen(n=100)
        sio.savemat(output_dir+"outputt_neg_sub1_"+str(i)+".mat" ,  {'surrogate_maps': surrogate_maps})
        
        brain_map =data_path +"outputt_pos_sub2_"+str(i)+".txt"  
        gen = Sampled(x=brain_map, D=filenames['D'], index=filenames['index'])
        surrogate_maps = gen(n=100)
        sio.savemat(output_dir+"outputt_pos_sub2_"+str(i)+".mat", {'surrogate_maps': surrogate_maps})
        
        brain_map =data_path +"outputt_neg_sub2_"+str(i)+".txt"  
        gen = Sampled(x=brain_map, D=filenames['D'], index=filenames['index'])
        surrogate_maps = gen(n=100)
        sio.savemat(output_dir+"outputt_neg_sub2_"+str(i)+".mat", {'surrogate_maps': surrogate_maps})
if __name__ == '__main__':
    output_path = '/home/cxpang/t_corr/'
    data_path = '/home/cxpang/t_corr/'
    brainsmash( output_path, data_path)
