project('matrix', 'c',
  version: '0.1.0',
  license: [ 'UNLICENSE' ],
  default_options: [ 'warning_level=3', 'build=debug' ]
)

inc = include_directories('dist')
sources = files('test/main.c', 'dist/matrix.c')
executable('test_lib', sources, include_directories: inc)
