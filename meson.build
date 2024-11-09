project('HeRoam', 'c')

# Source Files
srcs = files(
     'src/main.c',
     'src/vector.c',
     'src/simulation.c'
)

# Build Type options
debug_flags = ['-g', '-pg', '-Wall', '-Wextra']
optimized_flags = ['-O3', '-march=native', '-Wall', '-Wextra',
                    '-fopt-info-vec-optimized', '-fmax-errors=5']

# Flags based on buildtype
if get_option('buildtype') == 'debug'
  cflags = debug_flags
elif get_option('buildtype') == 'release'
  cflags = optimized_flags
endif

# Linker flags
debug_ldflags = ['-pg', '-lm', '-lgsl', '-lgslcblas']
optimized_ldflags = ['-lm', '-lgsl', '-lgslcblas']

# Flags based on buildtype
if get_option('buildtype') == 'debug'
   ldflags = debug_ldflags
elif get_option('buildtype') == 'release'
   ldflags = optimized_ldflags
endif

# Dependencies
omp_dep = dependency('openmp', required: false)

# Build
if omp_dep.found()
   exe = executable('HeRoam',
            srcs,
            dependencies: omp_dep,
            c_args: cflags,
            link_args: ldflags)
else
    exe = executable('HeRoam',
            srcs,
            c_args: cflags,
            link_args: ldflags)
endif
