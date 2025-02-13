from astropy.io import ascii
import numpy as np

outdir   = '/user/rsimons/submit_scripts/jobs'

cat = ascii.read('/user/rsimons/grizli_extractions/Catalogs/sample_cats/any_sample.cat')





fields = ['GS1','GS2', 'GS3', 'GS4', 'GS5', 'GN1', 'GN2', 'GN3', 'GN4', 'GN5', 'GN7', 'ERSPRIME']


sub_all = open(outdir + '/submit_all.sh', 'w+')
for field in fields:
    sub_all.write('sh ' + outdir + '/%s_submit_all.sh\n'%field)
    sf = open(outdir + '/%s_submit_all.sh'%field, 'w+')

    good = np.where(cat['col1'] == field)[0]
    for (field, di) in cat[good]:
        f = open(outdir + '/%s_%.5i.job'%(field, di), 'w+')

        f.write('Name = %s_%.5i_metalmcmc\n'%(field, di))
        f.write('Universe = Vanilla\n')
        f.write('Priority = 19\n')
        f.write('getenv = true\n')
        f.write('Executable = /home/rsimons/git/clear_local/metallicity/run_izi.py\n')
        f.write('Arguments = %s %.5i\n'%(field, di))
        f.write('Log = /user/rsimons/submit_scripts/logs/$(Name)_$(Cluster).log\n')
        f.write('Error = /user/rsimons/submit_scripts/logs/$(Name)_$(Cluster)_error.log\n')
        f.write('Output = /user/rsimons/submit_scripts/logs/$(Name)_$(Process).out\n')
        f.write('Queue\n')

        f.close()

        sf.write('condor_submit '+ outdir + '/%s_%.5i.job\n'%(field, di))
    sf.close()

sub_all.close()