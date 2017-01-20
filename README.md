# OgreTools

Linux: [![Build Status](https://travis-ci.org/nul0m/OgreTools.jl.svg?branch=master)](https://travis-ci.org/nul0m/OgreTools.jl)

Windows: [![Build status](https://ci.appveyor.com/api/projects/status/nvnlxyjhco48cjfs/branch/master?svg=true)](https://ci.appveyor.com/project/nul0m/ogretools-jl/branch/master)

Code coverage: [![Coverage Status](https://coveralls.io/repos/github/nul0m/OgreTools.jl/badge.svg?branch=master)](https://coveralls.io/github/nul0m/OgreTools.jl?branch=master)

OgreTools is a set of tools for macroeconomic modeling in Julia. Currently, it is not an official Julia package and supports only the latest stable version of Julia.

## Installation

If you have Internet connection on your machine it is enough to execute the following Julia commands to install the package:

```jl
    Pkg.clone("https://github.com/nul0m/OgreTools.jl.git")
    Pkg.build("OgreTools")
```

It will automatically download and install all the dependencies.

If there were no errors during installation the package is likely to be properly installed. You may check if it passes all the tests by running `Pkg.test(“OgreTools”)` or just try some examples from the next section right away.

## Examples

### Parsing

Parse a model file and create an object containing parsed model:

```jl
    # load the package
    using OgreTools
    
    # generate the path to the existing model file in the `tests` folder
    modpath = joinpath(Pkg.dir(),"OgreTools","test","test_model.mod")
    
    # by default parser ignores steady-state parts of the model equations
    mParsed = parseFile(modpath)

    # result is the parsed dynamic model
    print(mParsed)

    # to use steady-state part of the equation (when present) use isSstate option
    mParsedSstate = parseFile(modpath,isSstate=true)

    # in that case result contains steady-state parts of the equations 
    print(mParsedSstate)
```
### Solving parsed model

Parsed model object is one of the inputs for the solver problem. The other three are initial and terminal conditions and the date range.

Initial and terminal conditions are specified in a form of databases with timeseries for all the model symbols (variables, shocks and parameters). `DataBase`, `TimeSeries` and `Date` types are defined in this package. Here is how one can use them to define date range as well as initial and terminal conditions:

```jl
    # Define the solution range and create the database
    fd  = yy(2000)
    ld  = yy(3000)
    db  = DataBase(fd:ld)

    # Put variables/parameters into the database
    db["C"] = 1 # Constants are automatically expanded to the database range
    db["I"] = TimeSeries(fd:ld,1) # But it can be done manually as well
    db["K"] = 10
    db["A"] = 1

    db["alpha"] = 0.4
    db["beta"]  = 0.99
    db["gamma"] = 0.5
    db["delta"] = 0.9
    db["rho"]   = 0.5
    db["xi"]    = 0

    # We need two different databases:
    # 1. Initial conditions:
    db_init = deepcopy(db)
    # 2. Terminal conditions:
    db_targ = deepcopy(db)
```

Now we can define the solver problem and solve it. In this example we will use steady-state model `mParsedSstate` created in the parsing exercise above:

```jl
    # Adjust solution range for the lags and leads
    rng = (fd + mParsedSstate.maxlag) : (ld - mParsedSstate.maxlead)
    
    # Create the solver object with the steady-state model
    m_ss = SolveModel(mParsedSstate,db_init,db_targ,rng)

    # Tell solver to display iterations, by default nothing is displayed
    m_ss.display = "iter" 
    
    # Solve the steady-state model
    solve!(m_ss)
```

### Retrieving solution results

After the `solve!()` command solver object `m_ss` contains solution. You may retrieve it and print out the following way:

```jl
    # get the database with results and database with equations' discrepancies
    db_ss, db_resid_ss = get_result(m_ss)
    
    # Print the results
    print_rng   = yy(2490):yy(2510)
    print_vars  = ["C","I","K","A"]
    print(db_ss[print_vars][print_rng])
    
    # Print equations' discrepancies for the same range
    print(db_resid_ss[print_rng])
```

### Shock simulations and exogenize/endogenize

In the following example we will use dynamic model `mParsed` created in the parsing exercise above and create its solver object based on the resulting database from the steady-state model solution.

```jl
    # Define a smaller date range for this exercise
    rng = yy(2400):yy(2600)

    # Create a new dynamic solver problem
    m_exCK = SolveModel(mParsed,db_ss,db_ss,rng)

    # Swap exogenous and endogenous variables
    endogenize!(m_exCK, ["beta","delta"])
    exogenize!( m_exCK, ["C","K"])

    # Check out underlying model
    display(m_exCK.mod)

    # Set a shock (in the middle of the range, where the initial/terminal conditions are in steady-state)
    shock!(m_exCK,"C",yy(2500),1.6)

    # Set solver display option to show iterations and solve the model
    m_exCK.display = "iter"
    solve!(m_exCK)

    # Print the result
    db_sim_exCK, db_resid_exCK = get_result(m_exCK)
    print_vars = ["I","A","beta","delta"]
    print(db_sim_exCK[print_vars])
```

### Plotting

To be able to plot in Julia one have to use [`Plots`](https://github.com/JuliaPlots/Plots.jl) package, it is a plotting metapackage which brings many different plotting packages under a single API, making it easy to swap between plotting “backends”. There are many available “backends”, but we would recommend to use [`PlotlyJS`](https://github.com/sglyon/PlotlyJS.jl) package. For your convenience it's automatically installed together with OgreTools.

Here is a simple example of how to plot a line from the results of the exercise above:

```jl
    # load main plotting metapackage Plots
    using Plots

    # specify the "backend" to be used by the Plots package
    plotly()
    
    # plot timeseries "delta" from the resulting simulation database
    plot(db_sim_exCK["delta"])
```
