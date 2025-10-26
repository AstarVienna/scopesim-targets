# Defining spectral information

## Nomenclature
The name _spectrum_ was chosen over _SED_ for the following reasons, in no particular order:
* _SED_ is an acronym, wheras _spectrum_ is more explicit and easier understood.
* _SED_ (being an acronym) is usually capitalized, which goes against the common Python conventions of naming function arguments or class attributes (PEP 8).
  However, calling it _sed_ instead seems even less clear, especially because it can be pronounced in English like a word.
* In most contexts here, we're talking about the _spectrum_ of a star, planet or galaxy.
  _SED_ is more often used for the "larger scale shape" of the spectral energy distribution, which is often only coarsly sampled.

Side note: We consider _spectrum_ to be the singular case, and _spectra_ the plural form.
Thus, there is no such thing as "a spectra", because _spectra_ always means multiple.
If there is a case where more than one spectrum is meant (such as multiple point sources in a cluster), the word _spectra_ shall be used.

## Possible way to specify spectral information
In short:
* `synphot.SourceSpectrum` instance
* `astar_utils.SpectralType` instance, in which case the closest available template will be used.
* A string able to be parsed into a `astar_utils.SpectralType` instance, see above.
* "spex:spextra_name", where "spextra_name" resolves to a `spextra.Spextrum` database entry.

More TBA.
