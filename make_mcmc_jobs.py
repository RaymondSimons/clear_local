from astropy.io import ascii

outdir   = '/user/rsimons/submit_scripts/jobs'

cat = ascii.read('/user/rsimons/grizli_extractions/Catalogs/sample_cats/any_sample.cat')





fields = ['GS1','GS2', 'GS3', 'GS4', 'GS5', 'GN1', 'GN2', 'GN3', 'GN4', 'GN5', 'GN7', 'ERSPRIME']


for field in fields:
    sf = open(outdir + '/%s_submit_all.sh'%field, 'w+')

    good = np.where(cat['col1'] == field)[0]
    for (field, di) in cat[good]:
        f = open(outdir + '/%s_%.5i.job'%(field, di), 'w+')

        f.write('Name = %s_%.5i_metalmcmc\n'%(field, di))
        f.write('Universe = Vanilla\n')
        f.write('Priority = 19\n')
        f.write('getenv = true\n')
        f.write('Executable = /home/rsimons/git/clear_local/metal_maps_mcmc.py\n')
        f.write('Arguments = %s %.5i\n'%(field, di))
        f.write('Log = /user/rsimons/submit_scripts/logs/$(Name)_$(Cluster).log\n')
        f.write('Error = /user/rsimons/submit_scripts/logs/$(Name)_$(Cluster)_error.log\n')
        f.write('Output = /user/rsimons/submit_scripts/logs/$(Name)_$(Process).out\n')
        f.write('Queue\n')

        f.close()

        sf.write('qsub '+ outdir + '/%s_%.5i.job\n')
    sf.close()