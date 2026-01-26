# Defining spectral information
Any target that is not a full datacube needs some kind of spectral information.
This can range in complexity from a simple stellar spectral type, for which we have basic template spectra, all the way to a user-supplied spectrum.
This page gives a quick overview on how spectral information is handled in ScopeSim-Targets.

## Nomenclature
The name _spectrum_ was chosen over _SED_ for the following reasons, in no particular order:
* _SED_ is an acronym, whereas _spectrum_ is more explicit and easier understood.
* _SED_ (being an acronym) is usually capitalized, which goes against the common Python conventions of naming function arguments or class attributes (PEP 8).
  However, calling it _sed_ instead seems even less clear, especially because it can be pronounced in English like a word.
* In most contexts here, we're talking about the _spectrum_ of a star, planet or galaxy.
  _SED_ is more often used for the "larger scale shape" of the spectral energy distribution, which is often only coarsely sampled.

Side note: We consider _spectrum_ to be the singular case, and _spectra_ the plural form.
Thus, there is no such thing as "a spectra", because _spectra_ always refers to multiple.
If there is a case where more than one spectrum is meant (such as multiple point sources in a cluster), the word _spectra_ shall be used.

## Possible way to specify spectral information
In short:
* `synphot.SourceSpectrum` instance
* `astar_utils.SpectralType` instance, in which case the closest available template will be used.
* A string able to be parsed into a `astar_utils.SpectralType` instance, see above.
* "spex:spextra_name", where "spextra_name" resolves to a `spextra.Spextrum` database entry.
* "file:file_name", where "file_name" points to a local file containing spectral information, see below.

### Spectra from file
Currently, a "file:file_name" identifier will be forwarded as-is to `synphot.SourceSpectrum.from_file()`, meaning any format supported by that constructor is supported here.
Support for more formats is planned down the road.

## Wavelengths
Any spectra should be given in rest-frame and vacuum wavelengths.
All shifts are done by the tools (e.g. `ScopeSim`) using the target definitions.
Line-of-sight velocity or cosmological redshift ("z") can be supplied as part of the positional information, see [defining positions](defining_positions.md).
