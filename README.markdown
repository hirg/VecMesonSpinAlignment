# phi-meson Spin Alignment Code

### TreeProduction:
> - reconstrut event plane => ZDC-SMD (1st) and TPC (2nd)
> - reconstrut phi-meson => Same Event and Mixed Event
> - save event plane and reconstruted phi-meson into TTree

### FillSpinAlignment:
> - read in phi-meson TTree
> - boost K+ back into phi-meson rest frame
> - correlate phi-meson and event plane with K+ (cos(theta\*))
> - save histogram

### CalSpinAlignment:
> - subtract background from Mixed Event
> - extract raw spin alignment signal
> - apply TPC efficiency correction
> - apply MC event plane resolution correction
> - calculate systematic error

### Embedding
> - read in STAR embedding TTree
> - find MC tracks and corresponding RC tracks
> - Save all relavent information to TNtuple

### RcPhiEffCorr
> - read in TNtuple from Embedding process
> - apply same cut as picoDst data analysis
> - extract TPC efficiency for K+/K-
> - generate phi-meson through MC and let it decay with PHYTIA
> - apply K+/K- TPC effciency to decay daugthers and reconstruct phi-meson
> - extract TPC efficiency for phi-meson

### McPhiResCorr
> - generate phi-meson through MC with STAR published spectra and v2
> - let phi-meson decay with PHYTIA and boost K+ back to phi-meson rest frame
> - correlate phi-meson with fixed event plane (Psi = 0)
> - smear event plane with measured event plane resolution
> - correlate rho_00^phy with rho_00^obs to extract event plane resolution correction factor

### Utility
> - constant used in `VecMesonSpinAlignment`
> - functions used in `VecMesonSpinAlignment`
> - custom defined fype used in `VecMesonSpinAlignment`

### PlotMacro
> - macros to plot QA and final figures for rho_00 vs energy with sysErrors

### figures
> - place to save QA and final plots
