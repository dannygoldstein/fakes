# fakes

`$ pip install fakes`

inject realistic fakes onto ztf science images

requires `psfex >= 3.21.1`, `sextractor >= 2.18.0`, I recommend installing these via conda-forge:

`$ conda install -c conda-forge astromatic-source-extractor astromatic-psfex`

To see the tests pass:

`$ py.test`

To use:

```python
import fakes
from astropy.coordinates import SkyCoord

coord = SkyCoord(213.56123324, 38.02614853, unit='deg')
fakes.inject_psf(image='ztf_20181120510683_000718_zg_c01_o_q1_sciimg.fits', mag=15, coord=coord)
```

![fake image](https://user-images.githubusercontent.com/2769632/81816987-0e3f8480-94fa-11ea-81a5-ccd81cee8bf0.png)
