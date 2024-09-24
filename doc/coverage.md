# Coverage

Sometimes it can be useful to get an insight into which parts of the code are run for a certain test case. There is a compilation target specifically for this: `coverage`.

## Prerequisites

- GFortran compiler
- LCov coverage post processing tool

## How-to use

```shell
setup.pl -d=2 -arch=coverage
make
```

The `coverage` target compiles with `gfortran` and uses `gcov` to measure the coverage while running. After running, you can run `make cov` to post-process the collected data. For example,

```shell
mpirun -np=4 ./amrvac
make cov
```

This will generate a report in the `./html` directory. Run,

```shell
xdg-open html/index.html
```

Or point your browser to the correct location manually.

## Editor integration

Getting any of these to work properly may be a bit of a hassle:

- VSCode has several plugins that can handle LCov information for interactive use.
- Emacs has a `coverage-mode` minor mode available on MELPA.
- Vim has a plugin at [github.com/google/vim-coverage](https://github.com/google/vim-coverage).
- NeoVim has [andythigpen/nvim-coverage](https://github.com/andythigpen/nvim-coverage)
