import glob
from glob import glob
import PIL


dir_new = '/Users/rsimons/Desktop/clear/figures/linepng'
dir_old = '/Users/rsimons/Desktop/linepng'


out_dir = '/Users/rsimons/Desktop/clear/figures/linepng/out'
#imgs_comb = np.vstack((np.asarray( i.resize(min_shape) ) for i in imgs ))





fls = glob(dir_new + '/GN2_*_*.sed.png')


ids = unique(array([fl.strip('.sed.png')[-5:] for fl in fls]))



typ = 'sed'
for i, di in enumerate(ids):
    imgs1 = []
    for ii in arange(-1, 4):
        try: 
            im1 = PIL.Image.open(dir_new + '/GN2_%i_%s.%s.png'%(ii,di, typ))
        except:
            im1 = PIL.Image.fromarray(zeros((300,1100)), mode = 'RGBA')
        imgs1.append(im1)

    if len(imgs1) > 0:
        min_shape = sorted([(np.sum(i.size), i.size ) for i in imgs1])[0][1]
        imgs_comb1 = np.vstack((np.asarray( i.resize(min_shape) ) for i in imgs1 ))
        imgs_comb_new = PIL.Image.fromarray( imgs_comb1)

        imgs2 = []
        for ii in arange(-1, 4):
            try: im2 = PIL.Image.open(dir_old + '/GN2_%i_%s.%s.png'%(ii,di, typ))
            except: im2 = PIL.Image.fromarray(zeros((300,1100)), mode = 'RGBA')
            imgs2.append(im2)           

        if len(imgs2) > 0:
            min_shape = sorted([(np.sum(i.size), i.size ) for i in imgs2])[0][1]
            imgs_comb2 = np.vstack((np.asarray( i.resize(min_shape) ) for i in imgs2 ))
            imgs_comb_old = PIL.Image.fromarray( imgs_comb2)

            imgs = [imgs_comb_old, imgs_comb_new]
            min_shape = sorted([(np.sum(i.size), i.size ) for i in imgs])[0][1]
            imgs_comb = np.hstack((np.asarray( i.resize(min_shape) ) for i in imgs ))
            imgs_comb_all = PIL.Image.fromarray( imgs_comb)



            imgs_comb_all.save(out_dir + '/GN2_%s.%ss.png'%(di, typ))








