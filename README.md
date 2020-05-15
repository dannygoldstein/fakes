# fakes

`$ pip install fakes`

inject realistic fake point sources onto ztf science images using psfex to model the PSF

requires `psfex >= 3.21.1`, `sextractor >= 2.18.0`, I recommend installing these via conda-forge:

`$ conda install -c conda-forge astromatic-source-extractor astromatic-psfex`

if you get an error building galsim when installing via pip, install galsim via conda-forge:

`$ conda install -c conda-forge galsim`

To see the tests pass:

`$ py.test`

To use:

```python
import fakes
from astropy.coordinates import SkyCoord

coord = SkyCoord(213.56123324, 38.02614853, unit='deg')
fakes.inject_psf(image='ztf_20181120510683_000718_zg_c01_o_q1_sciimg.fits', mag=15, coord=coord)
```

To inject multiple fakes into a single image, `mag` and `coord` can be floating-point / coordinate arrays respectively. The program  reads in the pixel information from the first HDU of the passed image. It does not modify the original HDU. Instead, it writes the fake-injected image to an extension HDU. A fits binary table providing the truth table for the injected fakes is written to another extension HDU, following the fake-injected image.

![fake image](https://user-images.githubusercontent.com/2769632/81816987-0e3f8480-94fa-11ea-81a5-ccd81cee8bf0.png)

 
