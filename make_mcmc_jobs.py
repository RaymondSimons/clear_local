from astropy.io import ascii

outdir   = '/user/rsimons/submit_scripts/jobs'

cat = ascii.read('/user/rsimons/grizli_extractions/Catalogs/sample_cats/any_sample.cat')





for (field, di) in cat:

    print field, di

'''
f = open(outdir + '%s_%.5i.job'%(field, di))


Name = test_condor
Universe = Vanilla
Priority = 19
getenv = true
Executable = /home/rsimons/git/clear_local/test_condor.py
Arguments = 1 2 3 4
Log = logs/$(Name)_$(Cluster).log
Error = logs/$(Name)_$(Cluster)_error.log
Output = logs/$(Name)_$(Process).out
Queue

'''