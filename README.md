# acceleratum <img src="man/figures/logo.png" align="right" height="138" alt="" />

An R package to aid in bench calibration of triaxial accelerometery devices.  

Install with
```
remotes::install_github("loregalle/acceleratum",
                        build_vignettes = TRUE,
                        force = TRUE)
```

Get started with `vignette("acceleratum", package = "acceleratum")`  

Note that building the vignette requires a LaTeX installation. If
you don't have one you can try:
```
tinytex::install_tinytex()
```

Or you can avoid building the vignette and download the pdf from the
`rendered_vignettes` folder.

Currently in early development stages, be aware of bugs and possible
future changes.

# Scope of the package
`acceleratum` has been developed to aid with bench calibration of (triaxial)
accelerometer devices prior to field deployment. While originally developed
and used in a biologging context, its application are not restricted to it.  

A triaxial accelerometer needs calibration to correct systematic errors that
prevent raw measurements from matching true acceleration. In `acceleratum` we
consider three sources of errors that are structural to the device:

*   Scale (S): Each axis may respond differently to the same acceleration
(e.g., +1g on the x-axis might read as 0.98, while +1g on y reads as 1.02).
Calibration equalizes these sensitivities.
*   Misalignment (M): The sensor's internal axes may not be perfectly orthogonal
(90° apart). This causes acceleration along one axis to leak
into another.
*   Bias (b): Each axis can have a constant offset (e.g., reading 0.1g even
at rest), which must be subtracted to center measurements at zero.

`acceleratum` also allows to correct for mispositioning
of a deployed device due to user error (e.g. a device mounted backwards), and
unavoidable small differences in positioning of a device deployed on
living animals.

While bench calibration of a device is preferred, `acceleratum` also provides
methods to attempt the estimation of calibration parameters from
deployed devices.
