using OgreTools
using Base.Test

# path to the tests directory
const testPath = dirname(@__FILE__)
const testModelPath = joinpath(testPath, "test_model.mod")

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
      "K'n = delta*K(-1) + I(-1)*(1+(1-I(-1)/I(-2))^2)"
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
      "(db,t) -> db[t,4]-(db[t,7]*db[t-1,4]+db[t-1,3]*(1+(1-db[t-1,3]/db[t-2,3])^2))"
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
      "(db,t) -> -(db[t - 1,3]) * (-(-(db[t - 1,3])) / (db[t - 2,3] * db[t - 2,3])) * 2 * (1 - db[t - 1,3] / db[t - 2,3]) ^ (2 - 1)"
      "(db,t) -> -((1 + (1 - db[t - 1,3] / db[t - 2,3]) ^ 2 + db[t - 1,3] * 2 * (1 - db[t - 1,3] / db[t - 2,3]) ^ (2 - 1) * (-1 / db[t - 2,3])))"
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

    @test mParsed == mExpected

  end # testset "Parsing"

end # testset
