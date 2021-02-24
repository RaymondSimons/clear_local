import pidly
import os
import numpy as np
from numpy import *
from hri import hri
from joblib import Parallel, delayed

def izi(fluxes, errors, lines, logzprior = None, idl=None, dosave=False, savfile='res.sav', 
            grid=os.path.join(os.environ['IZI_DIR'],'grids','d13_kappa20.fits')) :

            #idl = pidly.IDL()
            print_output = False

            idl('fluxes = {0}'.format(np.array2string(fluxes, separator=',',max_line_width=1000)), print_output = print_output)
            idl('errors = {0}'.format(np.array2string(errors, separator=',',max_line_width=1000)), print_output = print_output)
            idl('lines = {0}'.format(np.array2string(lines, separator=',',max_line_width=1000)), print_output = print_output)
            idl('logzprior = {0}'.format(np.array2string(logzprior, separator=',',max_line_width=1000)).replace('\n', ''), print_output = print_output)
            idl('res=izi(fluxes, errors, lines, LOGZPRIOR=logzprior, NZ=100, gridfile="{0}")'.format(grid), print_output = print_output)
            if dosave :
                idl('save, file="{0}", res'.format(savfile))
            res = idl.ev('res', use_cache=True)
            return(res)




def run_parallel(test, run = 'R23'):

      z_arr = np.arange(8.5, 9.5, 0.01)
      prior = zeros(len(z_arr)) + 1.
      logzprior = vstack((z_arr, prior))
      idl_path = '/Applications/harris/idl87/bin/idl'#'/Applications/harris/idl87'#'/grp/software/Linux/itt/idl/idl84/idl/bin/idl'
      idl_path = '/grp/software/Linux/itt/idl/idl84/idl/bin/idl'
      idl = pidly.IDL(idl_path)
      np.random.seed(1)

      samples = 100

      SNs = [ 0.1, 1., 1.5,  2., 2.5,  3.,  3.5,  4. ,  4.5,  5. ,  5.5,  6. ,6.5,  7. ,  7.5,  8. ,  8.5,  9. ,  9.5, 10., 100.]

      lines_for_izi = np.array(['oii3726;oii3729', 'oiii4959;oiii5007', 'hbeta'])

      if test == 92:
            fluxes_intrinsic = np.array([1.0, 0.07, 1.0])      
      elif test == 90:
            fluxes_intrinsic = np.array([3.50, 0.5, 1.0])
      elif test == 88:
            fluxes_intrinsic = np.array([2.5, 4.0, 1.0])


      if run =='R3':
            lines_for_izi = lines_for_izi[1:]
            fluxes_intrinsic = fluxes_intrinsic[1:]

      Z = {}
      Z['case1'] = {}
      Z['case2'] = {}

      for key in Z.keys():
            for SN in SNs:
                  Z[key][SN] = {}
                  Z[key][SN] = {}
                  Z[key][SN]['Zmd'] = []
                  Z[key][SN]['Zle'] = []
                  Z[key][SN]['Zue'] = []
                  Z[key][SN]['Zp'] = []
                  Z[key][SN]['qmd'] = []
                  Z[key][SN]['qle'] = []
                  Z[key][SN]['que'] = []
                  Z[key][SN]['qp'] = []



      for SN in SNs:
            print (test, SN)
            for s in np.arange(samples):
                  errors_for_izi = fluxes_intrinsic/SN
                  fluxes_for_izi = np.array([fluxes_intrinsic[i] + np.random.normal(0, errors_for_izi[i]) for i in np.arange(len(fluxes_intrinsic))])

                  #fluxes_for_izi[fluxes_for_izi < 0.] = -666 
                  gd = np.where(fluxes_for_izi > 0.)[0]
                  if len(gd) > 1: 
                        res = izi(fluxes_for_izi[gd], errors_for_izi[gd], lines_for_izi[gd], logzprior=logzprior, idl=idl, dosave=False, savfile=None,
                                    grid=os.environ['IZI_DIR']+'/grids/d13_kappa20.fits')
                        (Zmd, Zlo, Zhi, Zp) = hri( res['zarr'][0], res['zpdfmar'][0])
                        (qmd, qlo, qhi, qp) = hri( res['qarr'][0], res['qpdfmar'][0])
                        #print (Zmd, Zlo, Zhi)
                        Z['case1'][SN]['Zmd'].append(Zmd)
                        Z['case1'][SN]['Zle'].append(Zmd - Zlo)
                        Z['case1'][SN]['Zue'].append(Zhi - Zmd)
                        Z['case1'][SN]['Zp'].append(Zp)
                        Z['case1'][SN]['qmd'].append(qmd)
                        Z['case1'][SN]['qle'].append(qmd - qlo)
                        Z['case1'][SN]['que'].append(qhi - qmd)
                        Z['case1'][SN]['qp'].append(qp)


                  errors_for_izi = fluxes_intrinsic/3.
                  errors_for_izi[-1] = fluxes_intrinsic[-1]/SN

                  fluxes_for_izi = np.array([fluxes_intrinsic[i] + np.random.normal(0, errors_for_izi[i]) for i in np.arange(len(fluxes_intrinsic))])
                  #fluxes_for_izi[fluxes_for_izi < 0.] = -666 
                  gd = np.where(fluxes_for_izi > 0.)[0]
                  if len(gd) > 1: 
                        res = izi(fluxes_for_izi[gd], errors_for_izi[gd], lines_for_izi[gd], logzprior=logzprior, idl=idl, dosave=False, savfile=None,
                                    grid=os.environ['IZI_DIR']+'/grids/d13_kappa20.fits')
                        (Zmd, Zlo, Zhi, Zp) = hri( res['zarr'][0], res['zpdfmar'][0])
                        (qmd, qlo, qhi, qp) = hri( res['qarr'][0], res['qpdfmar'][0])
                        #print (Zmd, Zlo, Zhi)
                        Z['case2'][SN]['Zmd'].append(Zmd)
                        Z['case2'][SN]['Zle'].append(Zmd - Zlo)
                        Z['case2'][SN]['Zue'].append(Zhi - Zmd)
                        Z['case2'][SN]['Zp'].append(Zp)
                        Z['case2'][SN]['qmd'].append(qmd)
                        Z['case2'][SN]['qle'].append(qmd - qlo)
                        Z['case2'][SN]['que'].append(qhi - qmd)
                        Z['case2'][SN]['qp'].append(qp)


      np.save('/user/rsimons/git/clear_local/metallicity/ref_response/izi_SN_Z%s_%s.npy'%(test, run), Z)



if __name__ == '__main__':
      Parallel(n_jobs = -1)(delayed(run_parallel)(test, run = 'R3') for test in [88, 90, 92])







