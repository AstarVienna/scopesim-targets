# YAML syntax for Targets

## Stellar
Some examples for stellar objects or groups of such.
For detailed explanations of the parameters, see [defining positions](defining_positions.md),
[defining spectra](defining_spectra.md) and  [defining brightness](defining_brightness.md) respectively.

### Single stars
```{literalinclude} example_yamls/stellar/star0.yaml
:name: yaml_stellar_star0
:caption: Simple single star offset from field center
```
```{literalinclude} example_yamls/stellar/star1.yaml
:name: yaml_stellar_star1
:caption: Simple single star with alternative position syntax
```
```{literalinclude} example_yamls/stellar/star2.yaml
:name: yaml_stellar_star2
:caption: Simple single star with full coordinates
```

### Binaries
```{literalinclude} example_yamls/stellar/binary0.yaml
:name: yaml_stellar_binary0
:caption: Binary with physical separation and distance
```
```{literalinclude} example_yamls/stellar/binary1.yaml
:name: yaml_stellar_binary1
:caption: Binary with angular separation
```
```{literalinclude} example_yamls/stellar/binary2.yaml
:name: yaml_stellar_binary2
:caption: Binaries with explicit positions can be defined simply as a star field with just two stars
```

### Star field
```{literalinclude} example_yamls/stellar/star_field0.yaml
:name: yaml_stellar_star_field0
:caption: Star field
```
```{literalinclude} example_yamls/stellar/star_field1.yaml
:name: yaml_stellar_star_field1
:caption: Star field
```
```{literalinclude} example_yamls/stellar/star_ield2.yaml
:name: yaml_stellar_star_field2
:caption: Star field including binaries
```
