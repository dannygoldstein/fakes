# fakes
inject realistic fakes onto ztf science images

```python
import fakes
from astropy.coordinates import SkyCoord

coord = SkyCoord(213.56123324, 38.02614853, unit='deg')
fakes.inject_psf(image='ztf_20181120510683_000718_zg_c01_o_q1_sciimg.fits', mag=15, coord=coord)
```

