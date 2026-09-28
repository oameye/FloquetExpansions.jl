# Changelog

All notable changes to this project are documented in this file. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Fixed

- Automatic dissipative-frame discovery seeds the frame with the Fourier harmonics of each
  channel operator, time average first, rather than with the time-dependent operator itself. A
  `jump` or `collapse` channel whose operator carries drive phases now completes with
  `positive_completion(vv, Gram())` without an explicit `DissipativeFrame`.
