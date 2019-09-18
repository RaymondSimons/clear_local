import matplotlib.pyplot as plt
import numpy as np
from numpy import *

plt.ioff()
plt.close('all')
plt.rcParams['text.usetex'] = True
fig, ax = plt.subplots(1,1,figsize = (7,7))




trimmed = True
do_min = True
if trimmed: figname = 'G102_G141_lines_trimmed.png'
else: figname = 'G102_G141_lines.png'
if trimmed:ax.set_xlabel('z$_{\t{line}}$', fontsize = 35)

else:ax.set_xlabel('z$_{\t{line}}$', fontsize = 25)


lines = [(r'Ly$\alpha$-1215' 		, 1215.4),
		('NV-1240' 		, 1240.81),
		('NIV-1487' 		, 1487.0),
		('CIV-1549' 		, 1549.48),
		('HeII-1640' 		, 1640.4),
		('OIII-1663' 		, 1665.85),
		('NIII-1750' 		, 1750.0),
		('CIII-1908' 		, 1908.734),
		('MgII-2799' 		, 2799.117),
		('NeV-3346' 		, 3346.8),
		('NeVI-3426' 		, 3426.85),
		('OII-3727,3729' 	, 3727.092),
		('NeIII-3867' 		, 3867.5),
		('H9-3890' 			, 3890.166),
		('H8-3971' 			, 3971.198),
		('Hd-4102' 			, 4102.892),
		('Hg-4341' 			, 4341.692),
		('OIII-4363' 		, 4364.436),
		(r'H$\beta$-4862' 			, 4862.71),
		('OIII-4960,5008' 	, 5008.24),
		('HeI-5877' 		, 5877.2),
		('OI-6302, 6363' 	, 6302.046),
		(r'H$\alpha$+NII-6563' 			, 6564.61),
		('SII-6718,6732' 	, 6718.29),
		('ArIII-7138' 		, 7138.0),
		('OII-7325' 		, 7322.0),
		('SIII-9531,9069' 	, 9530.6),
		('HeI-1083' 		, 10830.0),
		('PaB-12822' 		, 12821.6)]

if trimmed:
	lines = [('[OII]-3727,3729' 	, 3727.092),
			(r'H$\beta$-4862' 			, 4862.71),
			('[OIII]-4960,5008' 	, 5008.24),
			(r'H$\alpha$ + [NII]-6563' 			, 6564.61),
			('[SII]-6718,6732' 	, 6718.29)]



'''
lines = [('CIV .1550', 0.1550),
		 ('Fe II .2380', 0.2380),
		 ('Fe II .2586/600', 0.2586),
		 ('Mg II .2796/803', 0.2796),
		 ('[OII] .3728', 0.3728),
		 ('Ca II H .3934', 0.3934),
		 ('Ca II K .3969', 0.3969),
		 (r'H$\delta$ .4102', 0.4102),
		 (r'H$\gamma$ .4341', 0.4341),
		 ('[OIII] 0.4363', 0.4363),
		 (r'H$\beta$ .4862', 0.4862),
		 ('[OIII] .4960', 0.4960),
		 ('[OIII] .5007', 0.5007),
		 ('Mgb .5168/75/83', 0.5168),
		 (r'Fe$\lambda$ .5270/335', 0.5270),
		 ('Na D .5890/6', 0.5890),
		 ('[NII] .6548', 0.6548),
		 (r'H$\alpha$ .5663', 0.6563),
		 ('[NII] .6585', 0.6585),
		 ('[SII] .6717/31', 0.6717),
		 ('Ca II .8542/62/98', 0.8542),
		 ('[SII] .9096', 0.9096),
		 ('[SII] .9532', 0.9532),
		 ('[FeII] 1.257', 1.257),
		 (r'Pa$\beta$ 1.282', 1.282),
		 ('[FeII] 1.644', 1.644),
		 (r'Pa$\alpha$ 1.875', 1.875),
		 (r'H$_2$ 1-0 S(1) 2.212', 2.212)
		 ]
'''
lines = lines[::-1]

filters = [('G102', 8000, 11500, 'darkblue'),
		   ('G141', 10750, 17000, 'darkred')]

#For the regions where the filters overlap

filters_plot = [('G102', 8000, 11500, 'darkblue'),
		   		('G141', 10750, 17000, 'darkred')]



filters_overlap = [(10750, 11500, 'darkred')]

ytick_labels = np.array(['{:>12}'.format(line[0].split('-')[0])  for line in lines])


if do_min:
	ytick_labels = []
	yticks = []
	for l, line in enumerate(lines):
		#if ('OII-3727,3729' in line[0]) | (r'H$\alpha$' in line[0]) | ('OIII-4960,5008' in line[0]) | (r'H$\beta$-4862' in line[0]) | (r'SII-6718' in line[0]):
		ytick_labels.append('{:>12}'.format(line[0].split('-')[0]))
		yticks.append(l)
	ytick_labels = np.array(ytick_labels)
	yticks = np.array(yticks)

ax.set_yticks(yticks)
ax.set_ylim(-1, len(lines))

if trimmed: ax.set_yticklabels(ytick_labels, fontsize = 20)
else: ax.set_yticklabels(ytick_labels)




for f, filt in enumerate(filters_plot[::-1]):
	filt_name = filt[0]
	filt_lam0 = filt[1]
	filt_lam1 = filt[2]
	txt_color = filt[3]

	ax.annotate(filt_name, (0.75, 0.1 + 0.1*f), xycoords = 'axes fraction', fontsize = 40, fontweight = 'bold', color = txt_color)

	for l, line in enumerate(lines):
		#if ('OII-3727,3729' in line[0]) | (r'H$\alpha$' in line[0]) | ('OIII-4960,5008' in line[0]) | (r'H$\beta$-4862' in line[0]) | (r'SII-6718' in line[0]):
		line_lam_rest = line[1]
		z_min_filt = filt_lam0/line_lam_rest - 1.
		z_max_filt = filt_lam1/line_lam_rest - 1.
		xmn = log10((1+z_min_filt))
		xmx = log10((1+z_max_filt))
		ymn = (l + 0.5)/(len(lines)+1)
		ymx = (l + 1.5)/(len(lines)+1)


		if trimmed: fs = 12
		else: fs = 8
		if trimmed:ax.annotate('{:>12} \AA'.format(line[0].split('-')[1]), (-0.01, ymn+0.03), ha = 'right', va = 'center', xycoords = 'axes fraction', fontsize = fs, color = 'grey')
		else: ax.annotate('{:>12}'.format(line[0].split('-')[1]), (-0.12, np.mean([ymn, ymx])), ha = 'right', va = 'center', xycoords = 'axes fraction', fontsize = fs, color = 'grey')
		ax.axvspan(xmin = xmn, xmax = xmx, ymin = ymn, ymax = ymx, facecolor = txt_color, edgecolor = 'black', linewidth = 3)


for fo, filt_o in enumerate(filters_overlap):
	filt_lam0 = filt_o[0]
	filt_lam1 = filt_o[1]
	txt_color = filt_o[2]

	for l, line in enumerate(lines):
		#if ('OII-3727,3729' in line[0]) | (r'H$\alpha$' in line[0]) | ('OIII-4960,5008' in line[0]) | (r'H$\beta$-4862' in line[0]) | (r'SII-6718,6732' in line[0]):
		line_lam_rest = line[1]
		print (line)
		z_min_filt = filt_lam0/line_lam_rest - 1.
		z_max_filt = filt_lam1/line_lam_rest - 1.
		xmn = log10(z_min_filt + 1)
		xmx = log10(z_max_filt + 1)
		ymn = (l + 0.50)/(len(lines)+1)
		ymx = (l + 1.0)/(len(lines)+1)
		ax.axvspan(xmin = xmn, xmax = xmx, ymin = ymn, ymax = ymx, facecolor = txt_color, edgecolor = 'black', linewidth = 3)
		#This line is for cosmetics and is the most likely to break if the numb er of lines used changes
		ax.axvspan(xmin = xmx-0.01, xmax = xmx + 0.1, ymin = ymn+0.003, ymax = ymx-0.004, facecolor = txt_color, edgecolor = None, linewidth = 0)








if trimmed != True:
	ax.annotate('$\lambda$ (\AA)', (-0.12, 1.00), ha = 'right', xycoords = 'axes fraction', color = 'grey')
	ax.annotate('line', (-0.02, 1.00), ha = 'right', xycoords = 'axes fraction', fontweight = 'bold')


if trimmed: 
	ztck_min = 0.
	ztck_max = 4.
else:
	ztck_min = 0.
	ztck_max = 11.

ax.set_xlim(log10(1+ztck_min), log10(1.+ztck_max))
zticks = np.arange(ztck_min, ztck_max)

zticks_str = ['%i'%z for z in zticks]

zticks_log10_p1 = np.log10(zticks + 1)

ax.set_xticks(zticks_log10_p1)

if trimmed:ax.set_xticklabels(zticks_str, fontsize = 20)
else:ax.set_xticklabels(zticks_str, fontsize = 12)


if trimmed:fig.subplots_adjust(left = 0.20, bottom = 0.15, top = 0.95, right = 0.98)

else:fig.subplots_adjust(left = 0.20, bottom = 0.12, top = 0.95, right = 0.98)



fig.savefig('/Users/rsimons/Desktop/clear/figures/%s'%figname, dpi = 300)




