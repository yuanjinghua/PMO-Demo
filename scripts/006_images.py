'''
This script is for producing moment 0-2, NH2 and Tex 
maps from CO FITS cubes.

By Jinghua Yuan

Originally written in Dec 11th 2016
'''
import numpy as np
import aplpy as ap
import matplotlib as mpl
from matplotlib import pyplot as plt
from astropy import units as u
from astropy.io import fits, ascii
from spectral_cube import SpectralCube as sc

mpl.rc("font", family="serif", size=12)
mpl.rc("axes", linewidth = 1.0)
mpl.rc("lines", linewidth = 1.0)
mpl.rc("xtick.major", pad = 5, width = 1)
mpl.rc("ytick.major", pad = 5, width = 1)
mpl.rc("xtick.minor", width = 1)
mpl.rc("ytick.minor", width = 1)

fitsDir = '../fitsDir/'
epsDir  = '../epsDir/'
sourList = ascii.read('../sourList.txt')

def twoclick(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata

    print('x = %.2f, y = %.2f'%(
        ix, iy))

    # assign global variable to access outside of function
    global coord
    coord.append((ix, iy))

    # Disconnect after 2 clicks
    if len(coord) == 2:
        fig1.canvas.mpl_disconnect(cid0)
        plt.close(1)
    return

for isour in range(len(sourList)):
	sourName =  sourList['source'][isour]
	cubeCO = sc.read(fitsDir+sourName+'_co_cube.fits')
	cube13CO = sc.read(fitsDir+sourName+'_13co_cube.fits')
	cubeC18O = sc.read(fitsDir+sourName+'_c18o_cube.fits')

	cubeKMS_co = cubeCO.with_spectral_unit(u.km / u.s)
	cubeKMS_13co = cube13CO.with_spectral_unit(u.km / u.s)
	cubeKMS_c18o = cubeC18O.with_spectral_unit(u.km / u.s)

	rmsCO = np.std(cubeKMS_co.hdu.data[0:10,:,:])
	rms13CO = np.std(cubeKMS_13co.hdu.data[0:10,:,:])
	rmsC18O = np.std(cubeKMS_c18o.hdu.data[0:10,:,:])

	hdr = cubeKMS_co.hdu.header.copy()
	centPix = [int(hdr['CRPIX1']-1), int(hdr['CRPIX2']-1)]

	vels_co  = cubeKMS_co.spectral_axis.data
	slabDataCO = cubeKMS_co.hdu.data[:,centPix[1]-10:centPix[1]+11, 
	                                 centPix[0]-10:centPix[0]+11]
	inten_co = np.mean(slabDataCO, axis=1)
	inten_co = np.mean(inten_co, axis=1)
	
	vels_13co  = cubeKMS_13co.spectral_axis.data
	slabData13CO = cubeKMS_13co.hdu.data[:,centPix[1]-10:centPix[1]+11, 
	                                     centPix[0]-10:centPix[0]+11]
	inten_13co = np.mean(slabData13CO, axis=1)
	inten_13co = np.mean(inten_13co, axis=1)
	
	vels_c18o  = cubeKMS_c18o.spectral_axis.data
	slabDataC18O = cubeKMS_c18o.hdu.data[:,centPix[1]-10:centPix[1]+11, 
	                                     centPix[0]-10:centPix[0]+11]
	inten_c18o = np.mean(slabDataC18O, axis=1)
	inten_c18o = np.mean(inten_c18o, axis=1)

	fig0 = plt.figure(0, figsize = (10,7))
	ax0 = fig0.add_subplot(111)
	ax0.plot(vels_co, inten_co, drawstyle = 'steps-mid', 
			color = 'black')
	ax0.plot(vels_13co, inten_13co, drawstyle = 'steps-mid', 
				color = 'red')
	ax0.plot(vels_c18o, inten_c18o, drawstyle = 'steps-mid', 
			color = 'blue')
	ax0.set_xlim(min(vels_co),max(vels_co))
	plt.show(0)
	print('Please specify how many velocity components '+
		'you are interested in')
	nCom = int(input())
	plt.close(0)

	for i in range(nCom):
		print("Please click on the image to define the lower "+
			"and upper limits for component "+'%.2d' %(i+1))
		fig1 = plt.figure(1, figsize = (10,7))
		ax1 = fig1.add_subplot(111)
		ax1.plot(vels_co, inten_co, drawstyle = 'steps-mid', 
				color = 'black')
		ax1.plot(vels_13co, inten_13co, drawstyle = 'steps-mid', 
				color = 'red')
		ax1.plot(vels_c18o, inten_c18o, drawstyle = 'steps-mid', 
				color = 'blue')
		ax1.set_xlim(min(vels_co),max(vels_co))
		coord = []
		cid0 = fig1.canvas.mpl_connect('button_press_event', twoclick)
		plt.show(1)

		slabCO = cubeKMS_co.spectral_slab(coord[0][0]*u.km/u.s, coord[1][0]*u.km/u.s)
		slab13CO = cubeKMS_13co.spectral_slab(coord[0][0]*u.km/u.s, coord[1][0]*u.km/u.s)
		slabC18O = cubeKMS_c18o.spectral_slab(coord[0][0]*u.km/u.s, coord[1][0]*u.km/u.s)

		print("Producing moment 0, 1, and 2 maps")

		m0_co = slabCO.moment0()
		m0_13co = slab13CO.moment0()
		m0_c18o = slabC18O.moment0()

		masked13CO = slab13CO.with_mask(slab13CO>rms13CO*3 * u.K)
		m1_13co = masked13CO.moment1()
		m2_13co = masked13CO.moment2()

		print("producing column density and excitation temperature maps")

		Tex = m0_co.hdu.copy()
		Tex.header['BUNIT'] = 'K'
		Peak_co = np.max(slabCO.hdu.data, axis=0)
		Tex.data = 5.53/np.log(5.532/(Peak_co+0.819)+1)

		Nh2 = m0_13co.hdu.copy()
		Nh2.header['BUNIT'] = 'cm^-2'
		Nh2.data = (89*1e4*2.42e14*(Tex.data+0.92)/(1-np.exp(-5.29/Tex.data))*
					m0_13co.hdu.data/5.29/(1/(np.exp(5.29/Tex.data)-1)-
					1/(np.exp(5.29/2.73)-1)))

		print("saving fits files.")

		m0_co.write(fitsDir+sourName+'_'+'%.2d' %(i+1)+'_co_m0.fits',
			overwrite=True)
		m0_13co.write(fitsDir+sourName+'_'+'%.2d' %(i+1)+'_13co_m0.fits',
			overwrite=True)
		m0_c18o.write(fitsDir+sourName+'_'+'%.2d' %(i+1)+'_c18o_m0.fits',
			overwrite=True)
		m1_13co.write(fitsDir+sourName+'_'+'%.2d' %(i+1)+'_13co_m1.fits',
			overwrite=True)
		m2_13co.write(fitsDir+sourName+'_'+'%.2d' %(i+1)+'_13co_m2.fits',
			overwrite=True)
		Tex.writeto(fitsDir+sourName+'_'+'%.2d' %(i+1)+'_Tex.fits',
			clobber = True)
		Nh2.writeto(fitsDir+sourName+'_'+'%.2d' %(i+1)+'_Nh2.fits',
			clobber = True)

		fig2 = plt.figure(2, figsize = (6,6))
		fb1  = ap.FITSFigure(m0_co.hdu, figure = fig2)
		fb1.show_colorscale(vmin=0, cmap = 'viridis')
		fb1.show_colorbar(axis_label_text = 'K km s$^{-1}$')
		fb1.tick_labels.set_xformat("hhmmss")
		fb1.tick_labels.set_yformat("ddmm")
		fb1.show_beam()
		levs = np.linspace(0.3,0.9,7)*np.max(m0_13co.hdu.data)
		fb1.show_contour(m0_13co.hdu, levels = levs, colors = 'red')
		fb1.add_label(0, 1.03, sourName+'_'+'%.2d' %(i+1), relative = True, 
			horizontalalignment='left', color = "black")
		fb1.add_label(1, 1.03, '$^{13}$CO on $^{12}$CO', relative = True, 
			horizontalalignment='right', color = "black")

		fig2.savefig(epsDir+sourName+'_'+'%.2d' %(i+1)+'_13co_on_co.eps', 
			bbox_inches='tight', papertype='a2')
		fig2.savefig(epsDir+sourName+'_'+'%.2d' %(i+1)+'_13co_on_co.pdf', 
			bbox_inches='tight', papertype='a2')

		plt.close(2)

		fig2 = plt.figure(2, figsize = (6,6))
		fb1  = ap.FITSFigure(Tex, figure = fig2)
		fb1.show_colorscale(cmap = 'viridis')
		fb1.show_colorbar(axis_label_text = 'K')
		fb1.tick_labels.set_xformat("hhmmss")
		fb1.tick_labels.set_yformat("ddmm")
		fb1.show_beam()
		levs = np.linspace(0.3,0.9,7)*np.max(Nh2.data)
		fb1.show_contour(Nh2, levels = levs, colors = 'red')
		fb1.add_label(0, 1.03, sourName+'_'+'%.2d' %(i+1), relative = True, 
			horizontalalignment='left', color = "black")
		fb1.add_label(1, 1.03, 'N$_{H_2}$ on T$_{ex}$', relative = True, 
			horizontalalignment='right', color = "black")

		fig2.savefig(epsDir+sourName+'_'+'%.2d' %(i+1)+'_Nh2_on_Tex.eps', 
			bbox_inches='tight', papertype='a2')
		fig2.savefig(epsDir+sourName+'_'+'%.2d' %(i+1)+'_Nh2_on_Tex.pdf', 
			bbox_inches='tight', papertype='a2')

		plt.close(2)

		fig2 = plt.figure(2, figsize = (6,6))
		fb1  = ap.FITSFigure(Nh2, figure = fig2)
		fb1.show_colorscale(cmap = 'viridis')
		fb1.show_colorbar(axis_label_text = 'cm$^{-2}$')
		fb1.tick_labels.set_xformat("hhmmss")
		fb1.tick_labels.set_yformat("ddmm")
		fb1.show_beam()
		levs = np.linspace(0.3,0.9,7)*np.max(Nh2.data)
		fb1.show_contour(Nh2, levels = levs, colors = 'red')
		fb1.add_label(0, 1.03, sourName+'_'+'%.2d' %(i+1), relative = True, 
			horizontalalignment='left', color = "black")
		fb1.add_label(0.9, 1.03, 'N$_{H_2}$', relative = True, 
			horizontalalignment='right', color = "black")

		fig2.savefig(epsDir+sourName+'_'+'%.2d' %(i+1)+'_Nh2.eps', 
			bbox_inches='tight', papertype='a2')
		fig2.savefig(epsDir+sourName+'_'+'%.2d' %(i+1)+'_Nh2.pdf', 
			bbox_inches='tight', papertype='a2')

		plt.close(2)



		





    
