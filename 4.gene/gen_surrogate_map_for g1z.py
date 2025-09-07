from brainsmash.workbench.geo import volume
from brainsmash.mapgen.sampled import Sampled
import scipy.io as sio

coord_file = "D:/study/sub2/code/code/gene/coor.txt"
output_dir = "D:/study/sub2/gene/sub1"
filenames = volume(coord_file, output_dir)
brain_map = "D:/study/sub2/gene/sub1/output_mask.txt"
gen = Sampled(x=brain_map, D=filenames['D'], index=filenames['index'])
surrogate_g1 = gen(n=100)
sio.savemat('surrogate_g1.mat',{'surrogate_maps':surrogate_g1})