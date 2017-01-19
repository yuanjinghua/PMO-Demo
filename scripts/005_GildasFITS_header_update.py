'''
This script can modify header information 
of cube fits converted from gildas format.

By Jinghua Yuan

version   1: Nov 27th 2016
version 1.1: Dec 10th 2016

'''

from astropy.io import fits, ascii

fitsDir = '../fitsDir/'
sourList = ascii.read('../sourList.txt')

for isour in range(len(sourList)):
	sourName =  sourList['source'][isour]

	for line in ['co', '13co', 'c18o']:
		hduRaw = fits.open(fitsDir+sourName+'_'+line+'_cube.fits')
		hdrRaw = hduRaw[0].header
		
		xclip = int((hduRaw[0].data.shape[3] - 30)/2-0.5)
		yclip = int((hduRaw[0].data.shape[2] - 30)/2-0.5)
		data = hduRaw[0].data[0,:,xclip:-xclip,yclip:-yclip]
		hdu = fits.PrimaryHDU(data)
		keys = ['CTYPE1', 'CRVAL1', 'CDELT1', 'CRPIX1', 
			'CTYPE2', 'CRVAL2', 'CDELT2', 'CRPIX2', 'CTYPE3', 
			'CRVAL3', 'CDELT3', 'CRPIX3', 'OBJECT', 'RA', 'DEC', 
			'EQUINOX', 'LINE', 'RESTFREQ', 'BMAJ', 'BMIN', 'BPA']
		for ikey in keys:
			hdu.header[ikey] = hdrRaw[ikey]
	
		hdu.header['CTYPE1'] = 'RA---SFL'
		hdu.header['CTYPE2'] = 'DEC--SFL'
		hdu.header['CRPIX1'] = hdu.header['CRPIX1']-xclip
		hdu.header['CRPIX2'] = hdu.header['CRPIX2']-yclip
		hdu.header['CTYPE3'] = 'VELOCITY'
		hdu.header['CUNIT3'] = 'm/s'
		hdu.header['BUNIT']  = 'K'
	
		hdu.writeto(fitsDir+sourName+'_'+line+'_cube.fits', 
			clobber = True)

