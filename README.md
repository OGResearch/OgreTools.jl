# OgreTools

Linux: [![Build Status](https://travis-ci.org/nul0m/OgreTools.jl.svg?branch=master)](https://travis-ci.org/nul0m/OgreTools.jl)

Windows: [![Build status](https://ci.appveyor.com/api/projects/status/nvnlxyjhco48cjfs/branch/master?svg=true)](https://ci.appveyor.com/project/nul0m/ogretools-jl/branch/master)


Code coverage: [![Coverage Status](https://coveralls.io/repos/github/nul0m/OgreTools.jl/badge.svg?branch=master)](https://coveralls.io/github/nul0m/OgreTools.jl?branch=master)

## Installation

Install within Julia using

```jl
    Pkg.clone("https://github.com/nul0m/OgreTools.jl.git")
    Pkg.build("OgreTools")
```

## Testing

Test the package by running
```jl
    Pkg.test("OgreTools")
```

## Examples

### Parsing

Parse a model file and create an object containing parsed model:

```jl
using OgreTools
modpath = joinpath(Pkg.dir(),"OgreTools","test","test_model.mod")
# by default dynamic part of equations is used
mParsed = parseFile(modpath)
print(mParsed)
# to use steady-state part of the equation (when present) use isSstate option
mParsedSstate = parseFile(modpath,isSstate=true)
print(mParsedSstate)
```
