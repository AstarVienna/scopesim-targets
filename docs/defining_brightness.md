# Defining target brightness
All targets with separate spectral information (basically anything except a full datacube) usually also need some information on how to scale the spectrum.
The only exception to this are user-supplied spectra, which are already flux-calibrated to the desired level.

## Nomenclature
The name _brightness_ was chosen as the least ambigous of several rejected alternatives and because it matched with the ESO standard.
It is not usually associated with any particular unit or physical type, but rather more of a general "how bright is my source".

### Rejected alternatives
Some of these were used in `ScopeSim_Templates`, with no particular standard.

#### Magnitude
Physical flux units can also be used to specify brightness.
Calling it _magnitude_ would suggest it only accepts magnitudes.

#### Flux
The opposite as with _magnitude_, this would suggest using a magnitude value is _not_ supported.

#### Amplitude
This name seems more fitting for the scaling of normalized profiles, can be confusing for point sources.

#### Luminosity
The second-best alternative, but still more associated with specific units, like $L_\odot$.

## Possible way to specify brightness
### Band, Magnitude
Currently the only implemented way to specify target brightness.
Consists of a 2-tuple with a string identifying the bandpass ("filter") in which the magnitude is given, and a numerical value or `astropy.Quantity` in `u.mag`.
In this case, `mag` refers to Vega magnitudes, `ABmag` and `STmag` will be implemented in a future version.

### Band, Flux
_To be implemented._

### Wavelength, Flux
_To be implemented._

### Frequency, Flux
_To be implemented._

## To be discussed
Brightness from spectral type?

Absolute magnitudes?

Surface brightness vs integrated brightness?
