# Getting Started

## Dependencies

junopy uses several external open-source Python packages that are listed in the
`pyproject.toml` file under `dependencies`. If installing using pip, these
dependencies will be installed automatically, or you can choose to install them
individually.

## Installation

I use [`uv`](https://docs.astral.sh/uv/) for installing `junopy` and it is very 
easy to install this way. `uv` will create a virtual environment for this 
project, so you needn't worry about corrupting your global environment. (You
can also explicitly create a virtual environment using `uv venv`)


If you do not forsee making changes to the `junopy` source code, clone the git repository, `cd` into
the `junopy` directory and run - 
```
$ uv pip install .
```

If you do forsee making changes to `junopy`, use the `-e` flag. (for editable install)
```
$ uv pip install -e .
```

`uv` will install all necessary dependencies.

## Usage

I suggest you go through the examples directory that shows how you could use `junopy`. I store 
all `user` code in a different directory (which could be anywhere on your system), and use the 
`config.ConfigVars.set_root("<insert-path-to-juno-data-directory-here>")`. That way, I have all 
Juno data stored in one place for different projects.  An example of this is the examples directory
itself.
