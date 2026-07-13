project(
  'karpos-surfex',
  'c', 'fortran',
  version: '0.0.1.dev0',
  meson_version: '>=1.2',
  default_options: ['warning_level=1', 'buildtype=release'],
)

# Cython comme langage (compilation des .pyx)
add_languages('cython', native: false, required: true)

# --- Python + numpy ----------------------------------------------------------
py = import('python').find_installation(pure: false)
numpy_inc = run_command(
  py, '-c', 'import numpy; print(numpy.get_include())',
  check: true,
).stdout().strip()
inc_np = include_directories(numpy_inc)

# --- NetCDF (Fortran) via nf-config ------------------------------------------
nfconfig = find_program('nf-config', required: true)
nc_flibs = run_command(nfconfig, '--flibs', check: true).stdout().split()

# --- Flags Fortran alignés sur les Rules SURFEX ------------------------------
# Promotion real-8 + compat gfortran >=10 (cf. src/Rules.MCgfortran.mk).
fc = meson.get_compiler('fortran')
fflags = fc.get_supported_arguments(
  '-fdefault-real-8', '-fdefault-double-8', '-fno-second-underscore',
  '-fconvert=swap', '-fallow-argument-mismatch', '-fallow-invalid-boz',
)

# --- Couche C-ABI : shim + stub grib -----------------------------------------
capi = static_library(
  'sfx_capi',
  ['capi/grib_api_stub.F90', 'capi/sfx_capi.F90'],
  fortran_args: fflags,
  install: false,
)

# --- Lib SURFEX pré-construite (optionnelle tant que le shim n'est pas câblé) -
# Produite par la chaîne SURFEX (../build.sh) avec le stub grib_api.
surfex_lib_path = get_option('surfex_lib')
link_args = nc_flibs
if surfex_lib_path != ''
  link_args += [surfex_lib_path]
endif

# --- Extension Cython --------------------------------------------------------
ext = py.extension_module(
  '_surfex',
  'karpos_surfex/_surfex.pyx',
  include_directories: inc_np,
  link_with: capi,
  link_args: link_args,
  install: true,
  subdir: 'karpos_surfex',
)

py.install_sources(
  ['karpos_surfex/__init__.py'],
  subdir: 'karpos_surfex',
)
