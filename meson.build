project('matrix', 'c',
  version: '0.1.0',
  license: [ 'UNLICENSE' ],
  default_options: [ 'warning_level=3', 'build=debug' ]
)

inc = include_directories('dist')

test_executables = {
  'test-solve' : ['test/solve.c', 'dist/matrix.c'],
  'test-solve-numerical' : ['test/solve_numerical.c', 'dist/matrix.c'],
  'test-reduce-inverse' : ['test/reduce_inverse.c', 'dist/matrix.c'],
}

foreach name, sources : test_executables
  executable(name, sources, include_directories: inc)
endforeach
