using OgreTools
using Plots
using Base.Test

plotly()

# path to the tests directory
const testPath = dirname(@__FILE__)
const testModelPath = joinpath(testPath, "test_model.mod")

## todo: test exception handling

@testset "All tests" begin

  @testset "Basics" begin
    # default constructor
    @testset "- default" begin
      m0 = ParsedModel()
      @test isa(m0,ParsedModel)
      # overloaded comparison `==` for default values
      m0c = ParsedModel()
      @test m0 == m0c
    end # testset "Default"
    # from parser
    @testset "- parseFile" begin
      m1 = parseFile(testModelPath)
      @test isa(m1,ParsedModel)
      # overloaded comparison `==` for the really parsed models
      m1c = parseFile(testModelPath)
      @test m1 == m1c
    end # testset "parseFile"
    # overloaded methods for existing functions
    @testset "- overloaded" begin
      m2 = parseFile(testModelPath)
      # print()
      @test print(m2) == nothing
      # show()
      @test show(m2) == nothing
      # copy()
      m2c = copy(m2)
      @test m2 == m2c
      @test !is(m2,m2c)
      # deepcopy()
      m2c = deepcopy(m2)
      @test m2 == m2c
      @test !is(m2,m2c)
    end # testset "parseFile"
  end # testset "Basics"

  @testset "Parsing" begin
    # parse testing model file and create ParsedModel
    mParsed = parseFile(testModelPath)
    # manually create expected ParsedModel to be compared with
    mExpected = ParsedModel()
    mExpected.allnames    = ["A","C","I","K","alpha",
                             "beta","delta","gamma","rho","xi"]
    mExpected.endognames  = ["A","C","I","K"]
    mExpected.paramnames  = ["rho"]
    mExpected.exognames   = ["alpha","beta","delta","gamma","xi"]
    mExpected.eqs = [
      "C + I'n = A*K^alpha"
      "K'n = delta*K(-1) + I(-1)*(1-(1-I(-1)/I(-2))^2)"
      "C'n^(-gamma) = beta*C(+1)^(-gamma)*(A(+1)*alpha*K(+1)^(alpha-1) + delta)"
      "A'n = A(-1)^rho'p*exp(xi)"
    ]
    mExpected.eqlabs = [
      "Equilibrium"
      "Capital accumulation"
      "Euler equation"
      "Productivity process"
    ]
    mExpected.res_expr = map(parse,[
      "(db,t) -> db[t,2]+db[t,3]-(db[t,1]*db[t,4]^db[t,5])"
      "(db,t) -> db[t,4]-(db[t,7]*db[t-1,4]+db[t-1,3]*(1-(1-db[t-1,3]/db[t-2,3])^2))"
      "(db,t) -> db[t,2]^(-db[t,8])-(db[t,6]*db[t+1,2]^(-db[t,8])*(db[t,7]+db[t,5]*db[t+1,1]*db[t+1,4]^(db[t,5]-1)))"
      "(db,t) -> db[t,1]-(db[t-1,1]^db[t,9]*exp(db[t,10]))"
    ])
    mExpected.res_fun = map(eval,mExpected.res_expr)
    mExpected.jac_ind = [
        0 1 1;
        0 1 2;
        0 1 3;
        0 1 4;
        0 1 5;
       -2 2 3;
       -1 2 3;
       -1 2 4;
        0 2 4;
        0 2 7;
        1 3 1;
        0 3 2;
        1 3 2;
        1 3 4;
        0 3 5;
        0 3 6;
        0 3 7;
        0 3 8;
       -1 4 1;
        0 4 1;
        0 4 9;
        0 4 10
    ]
    mExpected.jac_expr = map(parse,[
      "(db,t) -> -(db[t,4] ^ db[t,5])"
      "(db,t) -> 1"
      "(db,t) -> 1"
      "(db,t) -> -(db[t,1]) * db[t,5] * db[t,4] ^ (db[t,5] - 1)"
      "(db,t) -> -(db[t,1]) * db[t,4] ^ db[t,5] * log(db[t,4])"
      "(db,t) -> -(db[t - 1,3]) * (-(-(-(db[t - 1,3])) / (db[t - 2,3] * db[t - 2,3])) * 2 * (1 - db[t - 1,3] / db[t - 2,3]) ^ (2 - 1))"
      "(db,t) -> -(((1 - (1 - db[t - 1,3] / db[t - 2,3]) ^ 2) + db[t - 1,3] * (-2 * (1 - db[t - 1,3] / db[t - 2,3]) ^ (2 - 1) * (-1 / db[t - 2,3]))))"
      "(db,t) -> -(db[t,7])"
      "(db,t) -> 1"
      "(db,t) -> -(db[t - 1,4])"
      "(db,t) -> -(db[t,6]) * db[t + 1,2] ^ -(db[t,8]) * db[t,5] * db[t + 1,4] ^ (db[t,5] - 1)"
      "(db,t) -> -(db[t,8]) * db[t,2] ^ (-(db[t,8]) - 1)"
      "(db,t) -> -((db[t,7] + db[t,5] * db[t + 1,1] * db[t + 1,4] ^ (db[t,5] - 1))) * db[t,6] * -(db[t,8]) * db[t + 1,2] ^ (-(db[t,8]) - 1)"
      "(db,t) -> -(db[t,6]) * db[t + 1,2] ^ -(db[t,8]) * db[t,5] * db[t + 1,1] * (db[t,5] - 1) * db[t + 1,4] ^ ((db[t,5] - 1) - 1)"
      "(db,t) -> -(db[t,6]) * db[t + 1,2] ^ -(db[t,8]) * (db[t,5] * db[t + 1,1] * db[t + 1,4] ^ (db[t,5] - 1) * log(db[t + 1,4]) + db[t + 1,1] * db[t + 1,4] ^ (db[t,5] - 1))"
      "(db,t) -> -(db[t + 1,2] ^ -(db[t,8])) * (db[t,7] + db[t,5] * db[t + 1,1] * db[t + 1,4] ^ (db[t,5] - 1))"
      "(db,t) -> -(db[t,6]) * db[t + 1,2] ^ -(db[t,8])"
      "(db,t) -> -1 * db[t,2] ^ -(db[t,8]) * log(db[t,2]) - (db[t,7] + db[t,5] * db[t + 1,1] * db[t + 1,4] ^ (db[t,5] - 1)) * db[t,6] * -1 * db[t + 1,2] ^ -(db[t,8]) * log(db[t + 1,2])"
      "(db,t) -> -(exp(db[t,10])) * db[t,9] * db[t - 1,1] ^ (db[t,9] - 1)"
      "(db,t) -> 1"
      "(db,t) -> -(exp(db[t,10])) * db[t - 1,1] ^ db[t,9] * log(db[t - 1,1])"
      "(db,t) -> -(db[t - 1,1] ^ db[t,9]) * exp(db[t,10])"
    ])
    mExpected.jac_fun = map(eval,mExpected.jac_expr)
    # assign all counters
    mExpected.nvars = 4
    mExpected.nexog = 5
    mExpected.nparams = 1
    mExpected.maxlag = 2
    mExpected.maxlead = 1

    @test mParsed == mExpected

  end # testset "Parsing"
  
  @testset "Date" begin
  
    @testset "Date-monthly" begin
  
      fd = mly(2017,11)
      ld = mly(2018,7)
      
      @test isa(fd, OgreTools.Date)
      
      @test print(fd) == nothing
      @test show(fd)  == nothing
      
      @test fd == deepcopy(fd)
      @test fd < ld
      @test ld - fd == 8
      @test ld == fd + (ld - fd)
      @test fd == ld - (ld - fd)
      
      year, per, freq = ypf(fd)
      @test year  == 2017
      @test per   == 11
      @test freq  == 12
      
      rng = fd:ld
      @test isa(rng, StepRange{OgreTools.Date,Int})
      @test (rng - fd)[1]   == 0
      @test (rng - fd)[end] == ld - fd
      
      str = dat2str(fd)
      @test str == "2017M11"
      
    end # Date-monthly
    
    @testset "Date-quarterly" begin
  
      fd = qly(2017,4)
      ld = qly(2018,3)
      
      @test isa(fd, OgreTools.Date)
      @test print(fd) == nothing
      @test show(fd)  == nothing
      
      @test fd == deepcopy(fd)
      @test fd < ld
      @test ld - fd == 3
      @test ld == fd + (ld - fd)
      @test fd == ld - (ld - fd)
      
      year, per, freq = ypf(fd)
      @test year  == 2017
      @test per   == 4
      @test freq  == 4
      
      rng = fd:ld
      @test isa(rng, StepRange{OgreTools.Date,Int})
      @test (rng - fd)[1]   == 0
      @test (rng - fd)[end] == ld - fd
      
      str = dat2str(fd)
      @test str == "2017Q4"
      
    end # Date-quarterly
    
    @testset "Date-yearly" begin
  
      fd = yly(2017)
      ld = yly(2021)
      
      @test isa(fd, OgreTools.Date)
      @test print(fd) == nothing
      @test show(fd)  == nothing
      
      @test fd == deepcopy(fd)
      @test fd < ld
      @test ld - fd == 4
      @test ld == fd + (ld - fd)
      @test fd == ld - (ld - fd)
      
      year, per, freq = ypf(fd)
      @test year  == 2017
      @test per   == 1
      @test freq  == 1
      
      rng = fd:ld
      @test isa(rng, StepRange{OgreTools.Date,Int})
      @test (rng - fd)[1]   == 0
      @test (rng - fd)[end] == ld - fd
      
      str = dat2str(fd)
      @test str == "2017"
      
    end # Date-yearly
    
    @testset "Date-errors" begin
    
      dm = mly(2017,1)
      dq = qly(2017,1)
      dy = yly(2017)
      
      @test_throws ErrorException dm == dq
      @test_throws ErrorException dq == dy
      @test_throws ErrorException dy == dm
      
      @test_throws ErrorException dm < dq
      @test_throws ErrorException dq < dy
      @test_throws ErrorException dy < dm
      
      @test_throws ErrorException dm - dq
      @test_throws ErrorException dq - dy
      @test_throws ErrorException dy - dm
    
    end  # Date-errors
  
  end # Date
  
  @testset "TimeSeries" begin
  
    fd = mly(1975,1)
    vals = randn(10)
    ts = TimeSeries(fd,vals)
  
    @testset "TimeSeries-Basic" begin
      
      @test isa(ts,TimeSeries)
      @test ts.firstdate == fd
      @test ts.values == vals
      @test ts == deepcopy(ts)
      @test print(ts) == nothing
      @test show(ts) == nothing
      @test isa(plot(ts), Plots.Plot)
    
    end # "TimeSeries-Basic"
    
    @testset "TimeSeries - set/get index" begin
      
      ind = 1
      val = 10
      ts[ind] = val
      @test all(ts.values[ind] .== val)
      ts1 = ts[ind]
      @test ts1 == TimeSeries(fd, ts.values[ind])
      
      ts[end] = 20
      @test ts.values[end] == 20
      
      ts[3:end-4] = 30
      @test all(ts.values[3:6] .== 30)
      
      d1 = mly(1975,2)
      d2 = mly(1975,5)
      
      ind = d1
      val = 40
      ts[ind] = val
      @test all(ts.values[ind-fd+1] .== val)
      ts1 = ts[ind]
      @test ts1 == TimeSeries(ind, ts.values[ind-fd+1])
      
      ind = d1:d2
      val = 50
      ts[ind] = val
      @test all(ts.values[ind-fd+1] .== val)
      ts1 = ts[ind]
      @test ts1 == TimeSeries(ind[1], ts.values[ind-fd+1])
      
    end # "TimeSeries - set/get index"
    
    @testset "TimeSeries-functions" begin
      
      # Basic function on [0,1]
      funs = [
        log,log1p,log2,log10,
        exp,expm1,
        abs,abs2,sqrt,cbrt,
        sin,cos,tan,cot,
        asin,acos,atan,acot,
        sinh,cosh,tanh,coth,
        asinh,atanh,
        erf,erfc,erfinv,erfcinv,
        gamma,lgamma,
        real,imag,fft
        ]
      
      ts = TimeSeries(fd,rand(10))
      
      for f in funs
        @test f(ts).values == f(convert(Array{Real},ts.values)) # For erf(c)inv
      end
      
      @test ifft(ts).values == ifft(convert(Array{Complex{Float64}}, ts.values)) # Here I couldn't find any other solution
      
      # Basic functions on [1,\infty]
      funs = [acosh, acoth]
      ts = TimeSeries(fd,1 + rand(10))
      for f in funs
        @test f(ts).values == f(ts.values)
      end
      
      # Statistical functions
      @test cumsum(ts).values   == cumsum(convert(Array{AbstractFloat},ts.values))
      @test cumprod(ts).values  == cumprod(ts.values) 
      
      funs = [mean, median, var, std]
      for f in funs
        @test f(ts) == f(ts.values)
      end
      
      # Function with 2 TimeSeries inputs
      ts1 = TimeSeries(fd,rand(10)) # For beta
      ts2 = TimeSeries(fd,rand(10)) # For beta
 
      @test (ts1 + ts2).values == ts1.values  + ts2.values
      @test (ts1 - ts2).values == ts1.values  - ts2.values
      @test (ts1 * ts2).values == ts1.values .* ts2.values
      @test (ts1 / ts2).values == ts1.values ./ ts2.values
      @test (ts1 ^ ts2).values == ts1.values .^ ts2.values
      
      @test cov(ts1,ts2) == cov(ts1.values,ts2.values)
      @test cor(ts1,ts2) == cor(ts1.values,ts2.values)
      @test beta(ts1,ts2).values  == beta(ts1.values,ts2.values)
      @test lbeta(ts1,ts2).values == lbeta(ts1.values,ts2.values)
      
      # Misc 1
      funs = [length, start, endof, length, isreal]
      for f in funs
        @test f(ts) == f(ts.values)
      end
      
      # Misc 2
      funs = [next, done]
      for f in funs
        @test f(ts,1) == f(ts.values,1)
      end
      
      # Misc 3
      vals = randn(100)
      ts = TimeSeries(fd,vals)
      ar =  0.9
      ma = -0.3
      @test filt([1; ma], [1; -ar], ts).values == filt([1; ma], [1; -ar], ts.values)
      
      # Misc 4
      @test diff(ts).values[2:end]  == diff(ts.values)
      @test pch(ts).values[2:end]   == 100*(ts.values[2:end] ./ ts.values[1:end-1] - 1)
      @test apch(ts).values[ts.firstdate.freq+1:end]  == 100*(ts.values[ts.firstdate.freq+1:end] ./ ts.values[1:end-ts.firstdate.freq] - 1)
      
      # Misc 5
      fd1   = mly(2017,1)
      vals1 = randn(10)
      ts1   = TimeSeries(fd1,vals1)
      
      fd2   = mly(2018,1)
      vals2 = randn(10)
      ts2   = TimeSeries(fd2,vals2)
      
      justify!(ts1,ts2)
      
      @test ts1.firstdate == ts2.firstdate == fd1
      @test length(ts1) == length(ts2)
      @test all(isnan(ts1.values[end-2:1]))
      @test all(isnan(ts2.values[1:12]))
      
      
    end # "TimeSeries-functions"
  
  end # TimeSeries
  
  @testset "DataBase" begin
  
    # Parse model
    mod_ss = parseFile(testModelPath, isSstate = true)
    
    # Parameters
    alpha = 0.4
    delta = 0.9
    beta  = 0.99
    gamma = 0.5
    rho   = 0.5

    # Steady-state
    # A = 1
    # K = ((1/beta-delta)/alpha)^(1/(alpha-1))
    # I = K*(1-delta)
    # C = K^alpha - I

    A = 1
    K = 10
    I = 1
    C = 1

    # Define the solution range and create the database
    fd  = yly(2000)
    ld  = yly(3000)

    db_init = DataBase(fd:ld)

    # Put variables/parameters into the database
    db_init["C"] = C # Constants are automatically expanded to the database range
    db_init["I"] = TimeSeries(fd:ld,I) # But it can be done manually as well
    db_init["K"] = K
    db_init["A"] = A

    db_init["alpha"] = alpha
    db_init["beta"]  = beta
    db_init["gamma"] = gamma
    db_init["delta"] = delta
    db_init["rho"]   = rho
    db_init["xi"]    = 0

    # Set the solution range (this is the longest possible)
    rng = (fd + mod_ss.maxlag) : (ld - mod_ss.maxlead)
    
    @testset "DataBase - basic" begin
    
      @test isa(db_init,DataBase)
      @test db_init             == deepcopy(db_init)
      @test show(db_init)       == nothing
      @test print(db_init)      == nothing
      @test print(db_init,0.1)  == nothing
      
    end # "DataBase - basic"
    
  end # "DataBase"

end # testset
