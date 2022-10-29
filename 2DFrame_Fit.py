import numpy as np
import scipy
import pylab as plt
import scipy.optimize
import scipy.io
import scipy.interpolate
import os

#Minimum number of frames for each interpolation iteration
k_low=[22,50,100,170,250]
#Intermedaite number of frames to be obtained in each interpolation iteration
k_up=[50,100,170,250,500]

#Iterating over the minimum number of frames list
for n in np.arange(0,len(k_low),1):
    #For the first iteration we create new dictionaries from the scratch
    if n==0:
        lista=np.arange(1,k_low[n]+1,1)
        #Creating Python dictonaries to store the intermediate frames
        data={}
        image_fin={}

        #Loop to read the available number of frames in each iteration to be interpolated
        for i in lista:
            #Change the corresponding directory where the images are stored with            #the name in your local directory. The intermediate frames will be 
            #stored in the directories created by the code
            data[i]=scipy.io.loadmat('cell_2_con_mask/con_'+str(k_low[n])+'/MaskTemp2_'+str(i)+'.mat')
            for k in data[i].keys():
                if type(data[i][k])==np.ndarray:
                    key_arr=k
            image_fin[i]=data[i][key_arr]
    #For subsequent iterations, we set image_fin as the previous im_fin
    #i.e. the available frames will be the intermediate frames from the 
    #previous iteration.
    else:
        image_fin=im_fin
        lista=np.arange(1,k_low[n]+1,1)
    vec_im=np.linspace(1,k_low[n],k_up[n])
    shape=np.shape(image_fin[1])
    im_fin={}
    for k in np.arange(1,len(vec_im)+1,1):
        im_fin[k]=np.zeros(shape)
    for i in np.arange(0,shape[0],1):
        for j in np.arange(0,shape[1],1):
            inter=np.array([])
            for k in lista:
                #Creating the 1D interpolators along the i,j position of each frame
                inter=np.append(inter,image_fin[k][i,j])
                inter[np.isnan(inter)]=0
            #Interpolating along the intermediate frames
            int_im=scipy.interpolate.interp1d(lista,inter)
            for l in np.arange(0,len(vec_im),1):
                #Storing the interpolated frames in the im_fin dictionary
                im_fin[l+1][i,j]=int_im(vec_im[l])

    print(k_low[n],'step')
    #Creating the directory, where the intermediate frames in each
    #iteration will be stored.
    os.system('mkdir cell_2_con_mask/con_'+str(k_up[n]))
    for i in np.arange(1,len(im_fin)+1,1):
        mdic={"mat":im_fin[i],"label":'mat'}
        scipy.io.savemat('cell_2_con_mask/con_'+str(k_up[n])+'/MaskTemp2_'+str(i)+'.mat',mdic)
