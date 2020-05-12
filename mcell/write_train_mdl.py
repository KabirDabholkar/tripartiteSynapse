from itertools import product
import os



#vdcc_num=[ 70.,  78.,  84.,  89.,  94.,  98., 103., 107., 111., 115., 120., \
#        124., 128., 133., 138., 143., 149., 157.]
vdcc_num=[ 69,  81,  89,  94,  97, 100, 103, 106, 110, 115, 119, 124, 129, \
       134, 140, 145, 151, 158, 169]

sims=["R150control","R150ER2x","R300ER2x","R150ER3x","R300ER3x"]

for sim,v,hz in product(sims,vdcc_num,[5,10,20]):
    path=sim+"_clamped"
    #if not os.path.exists(path):
    #    os.makedirs(path)
    with open(sim+"_clamped/RS20p%dhz.mdl" %hz,'r') as f:
        lines=f.readlines()
        
    with open(path+"/RS20p%d"%hz +"hz%d.mdl" %v,'w') as f:
        new_lines=lines
        for i,line in enumerate(new_lines):
            if "fname =" in line: 
                new_lines[i]='\tfname = "freq/'+sim+'/RS20p" & f & "hz"\n'
            if "/* Modifications in Parameters */" in line:
                new_lines.insert(i+1,"\tVDCC_number_presynaptic = %d\n" %v)
            new_lines[i]=new_lines[i].replace("/s_","/v%ds_" %v)
        f.writelines(new_lines)

