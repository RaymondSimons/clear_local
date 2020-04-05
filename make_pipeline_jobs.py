from astropy.io import ascii
import numpy as np

outdir   = '/user/rsimons/submit_scripts/jobs/grizli_pipeline'



fields = ['GS1','GS2', 'GS3', 'GS4', 'GS5', 'GN1', 'GN2', 'GN3', 'GN4', 'GN5', 'GN7', 'ERSPRIME']


sub_all = open(outdir + '/submit_all.sh', 'w+')
for field in fields:
    sub_all.write('qsub ' + outdir + '/%s.job'%(field))
    for (field, di) in cat[good]:
        f = open(outdir + '/%s.job'%(field), 'w+')

        f.write('Name = %s_pipeline\n'%(field))
        f.write('Universe = Vanilla\n')
        f.write('Priority = 19\n')
        f.write('getenv = true\n')
        f.write('Executable = /home/rsimons/git/clear_local/grizli_reduction.py\n')
        f.write('Arguments = --field %s --do_fit\n'%(field))
        f.write('Log = /user/rsimons/submit_scripts/logs/$(Name)_$(Cluster).log\n')
        f.write('Error = /user/rsimons/submit_scripts/logs/$(Name)_$(Cluster)_error.log\n')
        f.write('Output = /user/rsimons/submit_scripts/logs/$(Name)_$(Process).out\n')
        f.write('Queue\n')

        f.close()

sub_all.close()