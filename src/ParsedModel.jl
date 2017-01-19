type ParsedModel

  ### Properties
  # A string array (vector) containing the names of the endogenous variables
  endognames::Array{String,1}
  # A string array (vector) containing the names of all exogenous variables (excluding parameters)
  exognames::Array{String,1}
  # A string array (vector) containing the names of exogenous variables explicitly declared as parameters
  paramnames::Array{String,1}
  # Returns the list of all model symbols -- endogenous, exogenous and parameters -- ordered according to their indices
  allnames::Array{String,1}
  # A string array (vector) of the model equations, as given in the model file
  eqs::Array{String,1}
  # A string array (vector) of the model equation labels, as given in the model file
  eqlabs::Array{String,1}
  # A function array (vector), where the elements are anonymous functions that calculate the residuals of the equations. Each element has two variables, the database and the time period, and the output is the residual in that period, given data in the database.
  res_fun::Array{Function,1}
  # An expression array (vector), representing the functions from the res_fun
  res_expr::Array{Expr,1}
  # An integer array (matrix), giving the indices of the nonzero elements in the blocks of the Jacobian. It has three columns. The first stores the lag/lead, the second the row index (equation), the third is the column index (variable). In the first column, a negative element is a lag, a positive element is a lead, a zero is the contemporaneous (i.e. diagonal) block.
  jac_ind::Array{Int32,2}
  # A function array (vector), storing the functions that calculate the nonzeros in the Jacobian. The length is equal to the number of rows of jac_ind. The functions have the same structure as the functions in res_fun.
  jac_fun::Array{Function,1}
  # An expression array (vector)c representing the functions from the jac_fun
  jac_expr::Array{Expr,1}
  # Returns number of endogenous variables
  nvars::Int32
  # Returns number of exogenous variables
  nexog::Int32
  # Returns number of exogenous variables declared as parameters
  nparams::Int32
  # Returns maximum lag of endogenous variables
  maxlag::Int32
  #	Returns maximum lead of endogenous variables
  maxlead::Int32

  ### Constructor
  ParsedModel() = (x = new();
    x.endognames = Array{String,1}(0);
    x.exognames = Array{String,1}(0);
    x.paramnames = Array{String,1}(0);
    x.allnames = Array{String,1}(0);
    x.eqs = Array{String,1}(0);
    x.eqlabs = Array{String,1}(0);
    x.res_fun = Array{Function,1}(0);
    x.res_expr = Array{Expr,1}(0);
    x.jac_ind = Array{Int32,2}(0,0);
    x.jac_fun = Array{Function,1}(0);
    x.jac_expr = Array{Expr,1}(0);
    x.nvars = 0; # () -> length(x.endognames);
    x.nexog = 0; # () -> length(x.exognames);
    x.nparams = 0; # () -> length(x.paramnames);
    x.maxlag = 0; # () -> (length(x.jac_ind) > 0) ? -minimum(x.jac_ind[:,1]) : 0;
    x.maxlead = 0; # () -> (length(x.jac_ind) > 0) ? maximum(x.jac_ind[:,1]) : 0;
    x)

end # type ParsedModel

# New method for the standard `show` function
# to be used with the ParsedModel type
function Base.show(io::IO, pm::ParsedModel)
  println(io, "Object of type `ParsedModel`")
  println(io, "Endogenous variables ($(pm.nvars)): $(pm.endognames)")
  println(io, "Exogenous variables ($(pm.nexog)): $(pm.exognames)")
  println(io, "Parameters ($(pm.nparams)): $(pm.paramnames)")
  println(io, "Maximum lag: $(pm.maxlag)")
  println(io, "Maximum lead: $(pm.maxlead)")
end # function Base.show()

# New method for the standard `copy` function
# to be used with the ParsedModel type
function Base.copy(pm::ParsedModel)
  pmNew = ParsedModel()
  for fldName in fieldnames(pm)
    pmField = getfield(pm,fldName)
    if !isa(pmField,Function)
      setfield!(pmNew,fldName,copy(pmField))
    end # if
  end # for
  return pmNew
end # function Base.copy()

# New method for the standard `print` function
# to be used with the ParsedModel type
function Base.print(pm::ParsedModel)
  # reporting function to be used for model symbols (endo, exo, params)
  repFun = (list) -> begin
              # padding: number of symbols per line
              npad = 4
              # the longest symbol length
              ml = maximum(vcat(0,length.(list)))
              # loop over the groups of symbols (over lines)
              for i = 1 : Int(ceil(length(list)/npad))
                map((x)->print("\t"*lpad(x,ml," ")),
                    list[((i-1)*npad+1):min(i*npad,end)])
                println("\n")
              end # for
            end # repFun
  # skip 2 lines
  println("\n")
  # report endogenous variables
  println("Endogenous variables\n")
  repFun(pm.endognames)
  # report parameters
  println("Parameters\n")
  repFun(pm.paramnames)
  # report exogenous variables
  println("Exogenous variables\n")
  repFun(pm.exognames)
  # report equations
  println("Model equations\n")
  ml = maximum(vcat(0,length.(pm.eqlabs)))
  for i = 1:pm.nvars
    println("\t" * lpad(pm.eqlabs[i],ml," ") * ":\t" * pm.eqs[i] * "\n")
  end # for
end # function Base.print()

# New method for the standard `==` function
# to be used with the ParsedModel type
function Base.:(==)(pm1::ParsedModel,pm2::ParsedModel)
  print("Comparing objects of type `ParsedModel`: ")
  # compare numbers of endogenous variables
  if (pm1.nvars != pm2.nvars)
    println("models are different\n\tnvars: $(pm1.nvars) != $(pm2.nvars)")
    return false
  end
  # compare numbers of exogenous variables
  if (pm1.nexog != pm2.nexog)
    println("models are different\n\tnexog: $(pm1.nexog) != $(pm2.nexog)")
    return false
  end
  # compare numbers of parameters
  if (pm1.nparams != pm2.nparams)
    println("models are different\n\tnparams: $(pm1.nparams) != $(pm2.nparams)")
    return false
  end
  # compare maximum lags
  if (pm1.maxlag != pm2.maxlag)
    println("models are different\n\tmaxlag: $(pm1.maxlag) != $(pm2.maxlag)")
    return false
  end
  # compare maximum leads
  if (pm1.maxlead != pm2.maxlead)
    println("models are different\n\tmaxlead: $(pm1.maxlead) != $(pm2.maxlead)")
    return false
  end
  # compare names of model symbols
  if (sort(pm1.allnames) != sort(pm2.allnames))
    println("models are different\n\tallnames: $(pm1.allnames()) != $(pm2.allnames())")
    return false
  end
  # compare equations
  if (length(pm1.eqs) != length(pm2.eqs))
    println("numbers of equations are different\n\t# eqs: $(length(pm1.eqs)) != $(length(pm2.eqs))")
    return false
  end
  pm1Sorted = sort(map((x)->replace(x," ",""),pm1.eqs))
  pm2Sorted = sort(map((x)->replace(x," ",""),pm2.eqs))
  if (pm1Sorted != pm2Sorted)
    ix = (pm1Sorted .!= pm2Sorted)
    println("models are different\n\teqs: $(pm1Sorted[ix]) != $(pm2Sorted[ix])")
    return false
  end
  # compare equations' expressions
  if (pm1.res_expr != pm2.res_expr)
    ix = (pm1.res_expr .!= pm2.res_expr)
    println("models are different\n\tresidual expr: $(pm1.res_expr[ix]) != $(pm2.res_expr[ix])")
    return false
  end
  # compare number of non-zero derivatives in Jacobian
  if (size(pm1.jac_ind) != size(pm2.jac_ind))
    println("models are different\n\tsize(Jacobian indices): $(size(pm1.jac_ind)) != $(size(pm2.jac_ind))")
    return false
  end
  # compare rows of non-zero derivatives in Jacobian
  rows1 = [ view(pm1.jac_ind,i,1:size(pm1.jac_ind,2)) for i=1:size(pm1.jac_ind,1) ]
  pm1Order = sortperm(rows1; order=Base.Order.Lexicographic)
  pm1Sorted = pm1.jac_ind[pm1Order,:]
  rows2 = [ view(pm2.jac_ind,i,1:size(pm2.jac_ind,2)) for i=1:size(pm2.jac_ind,1) ]
  pm2Order = sortperm(rows2; order=Base.Order.Lexicographic)
  pm2Sorted = pm2.jac_ind[pm2Order,:]
  if (pm1Sorted != pm2Sorted)
    ix = (pm1Sorted .!= pm2Sorted)
    println("models are different\n\tJacobian indices (rows): $(pm1Sorted[ix,:]') != $(pm2Sorted[ix,:]')")
    return false
  end
  # compare derivatives' expressions
  pm1Sorted = pm1.jac_expr[pm1Order,:]
  pm2Sorted = pm2.jac_expr[pm2Order,:]
  if (pm1Sorted != pm2Sorted)
    ix = (pm1Sorted .!= pm2Sorted)
    println("models are different\n\tJacobian expr: $(pm1Sorted[ix]) != $(pm2Sorted[ix])")
    return false
  end

  println("models are identical")
  return true
end

# Main function reading the content of the specified model file
# and calling necessary methods from the C-library to parse it
function parseFile(fileName::String; isSstate::Bool=false)
  # check if :fileName exists
  if !isfile(fileName)
    error("File `$(fileName)` not found")
  end # if

  # instantiate the model
  mod = ParsedModel()

  # get all the necessary model characteristics
  # from the parser object
  eqsRef = Ref{Ptr{Ptr{Cchar}}}(0)
  eqLabelsRef = Ref{Ptr{Ptr{Cchar}}}(0)
  eqFormulasRef = Ref{Ptr{Ptr{Cchar}}}(0)
  nEqtnRef = Ref{Cint}(0)
  varsRef = Ref{Ptr{Ptr{Cchar}}}(0)
  varTypesRef = Ref{Ptr{Ptr{Cchar}}}(0)
  nVarsRef = Ref{Cint}(0)
  derFormulasRef = Ref{Ptr{Ptr{Cchar}}}(0)
  derIndicesRef = Ref{Ptr{Ptr{Cint}}}(0)
  nDersRef = Ref{Cint}(0)
  errorMessageRef = Ref{Ptr{Cchar}}(0)
  retCode = ccall((:parse_model, _jl_libparser), Cint,
    (Cstring, Cint,
     Ref{Ptr{Ptr{Cchar}}},Ref{Ptr{Ptr{Cchar}}},Ref{Ptr{Ptr{Cchar}}},Ref{Cint},
     Ref{Ptr{Ptr{Cchar}}},Ref{Ptr{Ptr{Cchar}}},Ref{Cint},
     Ref{Ptr{Ptr{Cchar}}},Ref{Ptr{Ptr{Cint}}},Ref{Cint},
     Ref{Ptr{Cchar}}),
    fileName, isSstate,
    eqsRef, eqLabelsRef, eqFormulasRef, nEqtnRef,
    varsRef, varTypesRef, nVarsRef,
    derFormulasRef,derIndicesRef,nDersRef,
    errorMessageRef)
  if (retCode == 0)
    # assign equations and their labels
    mod.eqs = map((x)->strip(unsafe_string(x)),unsafe_wrap(Array, eqsRef[], nEqtnRef[], true))
    mod.eqlabs = map((x)->strip(unsafe_string(x)),unsafe_wrap(Array, eqLabelsRef[], nEqtnRef[], true))
    # process and assign equations' expressions
    eqStrings = map(unsafe_string,unsafe_wrap(Array, eqFormulasRef[], nEqtnRef[], true))
    # convert equation strings to expressions and anonymous functions
    mod.res_expr = map((x)->parse("(db,t)->"*x),eqStrings)
    mod.res_fun = map(eval,mod.res_expr)
    # assign variables in accordance with their types
    vars = map(unsafe_string,unsafe_wrap(Array, varsRef[], nVarsRef[], true))
    varTypes = map(unsafe_string,unsafe_wrap(Array, varTypesRef[], nVarsRef[], true))
    for i = 1 : nVarsRef[]
      if (varTypes[i] == "endogenous")
        push!(mod.endognames,vars[i])
      elseif (varTypes[i] == "exogenous")
        push!(mod.exognames,vars[i])
      elseif (varTypes[i] == "parameter")
        push!(mod.paramnames,vars[i])
      else
        error("Variable `$(vars[i])` is of unknown type `$(varTypes[i])`")
      end # if
      # push the name to the list of all model variables
      push!(mod.allnames,vars[i])
    end # for
    # process and assign derivatives' expressions
    derStrings = map(unsafe_string,unsafe_wrap(Array, derFormulasRef[], nDersRef[], true))
    # convert derivatives' strings to expressions and anonymous functions
    mod.jac_expr = map((x)->parse("(db,t)->"*x), derStrings)
    mod.jac_fun = map(eval, mod.jac_expr)
    # process and assign the array of non-zero jacobian indices
    jac_arr = map((x)->unsafe_wrap(Array,x,3,true),
                  unsafe_wrap(Array,derIndicesRef[],nDersRef[],true))
    mod.jac_ind = hcat(jac_arr...)'
    # assign all counters
    mod.nvars = length(mod.endognames)
    mod.nexog = length(mod.exognames)
    mod.nparams = length(mod.paramnames)
    if (length(mod.jac_ind) == 0)
      mod.maxlag = 0
      mod.maxlead = 0
    else
      nlags = sum(mod.jac_ind[:,1] .< 0)
      nleads = sum(mod.jac_ind[:,1] .> 0)
      mod.maxlag = (nlags > 0) ? -minimum(mod.jac_ind[:,1]) : 0
      mod.maxlead = (nleads > 0) ? maximum(mod.jac_ind[:,1]) : 0
    end # if

    return mod
  else
    # there was an error while parsing the model
    errMes = unsafe_string(errorMessageRef[])
    error("Unable to parse `$(fileName)`!\n Error code: $(retCode)\n Error message: $(errMes)")
  end # if
end # function ParseFile()
