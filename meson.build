project('HeRoam', 'c')

# Source Files
main_srcs = files(
     'src/main.c',
     'src/vector.c',
     'src/h5file.c',
     'src/simulation.c'
)

reader_srcs = files(
     'src/reader.c',
     'src/h5file.c',
)


# Build Type options
debug_flags = ['-g', '-pg', '-Wall', '-Wextra']
optimized_flags = ['-O3', '-march=native', '-Wall', '-Wextra',
                    '-fopt-info-vec-optimized', '-fmax-errors=5',
                    '-ffast-math', '-funroll-loops']

# Flags based on buildtype
if get_option('buildtype') == 'debug'
  cflags = debug_flags
elif get_option('buildtype') == 'release'
  cflags = optimized_flags
endif

# Linker flags
debug_ldflags = ['-pg', '-lm']
optimized_ldflags = ['-lm']

# Flags based on buildtype
if get_option('buildtype') == 'debug'
   ldflags = debug_ldflags
elif get_option('buildtype') == 'release'
   ldflags = optimized_ldflags
endif

# Dependencies
omp_dep = dependency('openmp', required: false)
h5_dep = dependency('hdf5', language: 'c')

# Build HeRoam binary
if omp_dep.found()
   executable('HeRoam',
            main_srcs,
            dependencies: [omp_dep, h5_dep],
            c_args: cflags,
            link_args: ldflags)
else
    executable('HeRoam',
            main_srcs,
            dependencies: h5_dep,
            c_args: cflags,
            link_args: ldflags)
endif

# Build h5read binary
executable('h5read',
          reader_srcs,
          dependencies: [omp_dep, h5_dep],
          c_args: cflags,
          link_args: ldflags)
