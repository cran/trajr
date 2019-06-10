# Changes to the `trajr` package

## trajr 1.3.0

* Detect and report some invalid parameter values.
* Altered handling of parameter `readcsvFn` to `TrajsBuild` to make it 
  possible to use `readr::read_csv` without a wrapper function.
* Added function `TrajResampleTime` to resample a trajectory to fixed step times.
* Added parameters `start.pt.pch` and `start.pt.col` to plotting functions.
* Added parameter `dt` to `TrajTranslate`
* Fix vertical extents of rectangles in `plot.TrajSpeedIntervals` to handle non-default ylim values.
* Added optional progressbar to `TrajsMergeStats`.
* TrajsMergeStats now passes the arguments stringsAsFactors = FALSE to rbind. This prevents incorrect 
  behaviour and the warning "invalid factor level, NA generated" if one or more of your statistics are characters.
* Enhanced `TrajRotate` to allow absolute rotation and arbitrary origin of rotation.

## trajr 1.2.0

* Added start.pt.cex parameter to function `lines.Trajectory`.
* Added function `TrajConvertTime`.

## trajr 1.1.0

* Added correct citation.
* Fixed: `plot.TrajSpeedIntervals` was not passing additional arguments (`...`) to `plot`.
* Added: functions `TrajDuration`, `TrajMeanVelocity`, `TrajTranslate`.
* Added `translateToOrigin` parameter to function `TrajsBuild`.