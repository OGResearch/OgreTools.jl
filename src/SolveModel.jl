
type SolveModel

  mod::ParsedModel
  endog_ind::Array{Int,1}
  stacked_jac_ind::Array{Int,2}
  stacked_jac_fun::Array{Function,1}
  stacked_jac_expr::Array{Expr,1}
  db_init::Array{Number,2}
  db_targ::Array{Number,2}
  db_resid::Array{Number,2}
  range
  firstdate::Date
  daterange
  tolerance::Float64
  maxiter::Int64
  display::String
  resnorm::Float64
  iter::Int64
  time::DateTime

end

function SolveModel(m::ParsedModel,db_init::DataBase,db_targ::DataBase,rng)

  stacked_jac = get_stacked_jac(m,rng - db_targ.firstdate + 1)

  m = SolveModel(
      deepcopy(m), # So that end(x)ogenizing wouldn't change the parsed model in the workspace
      get_endog_ind(m),
      stacked_jac[1],
      stacked_jac[2],
      stacked_jac[3],
      db2array(db_init,m.allnames),
      db2array(db_targ,m.allnames),
      repmat([NaN],db_targ.lastdate-db_targ.firstdate+1, convert(Int,m.nvars)), # Need to get the types correct
      rng - db_targ.firstdate + 1,
      db_targ.firstdate,
      rng,
      1e-9,
      100,
      "iter",
      NaN,
      0,
      DateTime()
    )

  return m

end

function SolveModel(m::ParsedModel,db::DataBase,rng)
  return SolveModel(m, db, db, rng)
end

function ==(m1::SolveModel,m2::SolveModel)

  y = true
  fns = fieldnames(m1)
  for fn in fns
    y = y && getfield(m1,fn) == getfield(m2,fn)
  end
  return y
  
end

function Base.show(io::IO,m::SolveModel)

  w = 20

  print(io,"\n")

  print(io,lpad("Model: ",w," "))
  print(io,"$(length(m.mod.eqs)) equation(s), $(length(m.mod.paramnames)) parameter(s), $(length(m.mod.exognames)) exogenous variable(s)\n")

  print(io,lpad("Database: ",w," "))
  print(io,"first date: $(m.firstdate), last date: $(m.firstdate + size(m.db_targ,1) - 1), number of variables: $(size(m.db_targ,2))\n")

  print(io,lpad("Solution range: ",w," "))
  print(io,"$(m.daterange) ($(length(m.daterange)) periods)\n")

  print(io,lpad("Last solved: ",w," "))
  print(io,"$(m.time)\n")

  print(io,lpad("Max norm of resid: ",w," "))
  print(io,"$(m.resnorm)\n")

  print(io,lpad("Iterations: ",w," "))
  print(io,"$(m.iter)\n")

  print(io,lpad("Tolerance: ",w," "))
  print(io,"$(m.tolerance)\n")

  print(io,lpad("Maxiter: ",w," "))
  print(io,"$(m.maxiter)\n")

end

function Base.print(m::SolveModel)

  print(m.mod)

end


function get_endog_ind(m::ParsedModel)
  
  endog_ind = Array{Int}(0)
  
  for name in m.endognames
    append!(endog_ind,findin(m.allnames,[name]))
  end
  
  return endog_ind
  
end

function get_stacked_jac(m::ParsedModel, rng)

  # Initialize output
  ro    = Array{Int,1}(0)
  co    = Array{Int,1}(0)
  times = Array{Int,1}(0)
  funs  = Array{Function,1}(0)
  exprs = Array{Expr,1}(0)
  
  # Keep only derivatives w/resp to endogenous variables
  endog_ind = get_endog_ind(m)
  ind = findin(m.jac_ind[:,3], endog_ind)
  jac_ind   = m.jac_ind[ind,:]
  jac_fun   = m.jac_fun[ind]
  jac_exprs = m.jac_expr[ind]
  
  # Replace "indeces in database" with "indeces in Jacobian block"
  for i = 1:size(jac_ind,1)
    jac_ind[i,3] = findin(endog_ind,jac_ind[i,3])[1]
  end

  # de_t/dx_t

  ind = findin(jac_ind[:,1],0)
  tjac_ind    = jac_ind[ind,2:3]
  tjac_fun    = jac_fun[ind]
  tjac_exprs  = jac_exprs[ind]

  for t = rng[1]:rng[end]
    for k = 1:size(tjac_ind,1)
      push!(funs,   tjac_fun[k])
      push!(exprs,  tjac_exprs[k])
      push!(times,t)
    end
    append!(ro, tjac_ind[:,1] + (t - rng[1])*m.nvars)
    append!(co, tjac_ind[:,2] + (t - rng[1])*m.nvars)
  end

  # de_t/dx_{t-k}
  for l = 1:m.maxlag

    ind = findin(jac_ind[:,1],-l)
    tjac_ind    = jac_ind[ind,2:3]
    tjac_fun    = jac_fun[ind]
    tjac_exprs  = jac_exprs[ind]

    for t = rng[1]+l:rng[end]
      for k = 1:size(tjac_ind,1)
        push!(funs,   tjac_fun[k])
        push!(exprs,  tjac_exprs[k])
        push!(times,t)
      end
      append!(ro, tjac_ind[:,1] + (t - rng[1])*m.nvars)
      append!(co, tjac_ind[:,2] + (t - rng[1])*m.nvars - l*m.nvars)
    end

  end

  # de_t/dx_{t+k}
  for l = 1:m.maxlead

    ind = findin(jac_ind[:,1],l)
    tjac_ind    = jac_ind[ind,2:3]
    tjac_fun    = jac_fun[ind]
    tjac_exprs  = jac_exprs[ind]

    for t = rng[1]:rng[end]-l
      for k = 1:size(tjac_ind,1)
        push!(funs,   tjac_fun[k])
        push!(exprs,  tjac_exprs[k])
        push!(times,t)
      end
      append!(ro, tjac_ind[:,1] + (t - rng[1])*m.nvars)
      append!(co, tjac_ind[:,2] + (t - rng[1])*m.nvars + l*m.nvars)
    end

  end
  
  return [ro co times], funs, exprs

end

function eval_equations(m::SolveModel)

  res = Array{Float64}(0)
  for t = m.range
    for j = 1:m.mod.nvars
      push!(res, m.mod.res_fun[j](m.db_targ,t))
    end
  end
  res = reshape(res,length(res),1)
  return res

end

function eval_derivatives(m::SolveModel)

  ro    = m.stacked_jac_ind[:,1]
  co    = m.stacked_jac_ind[:,2]
  times = m.stacked_jac_ind[:,3]
  vals  = Array{Float64,1}(length(ro))
  
  # display(m.db_targ)
  
  for i = 1:length(vals)
    vals[i] = m.stacked_jac_fun[i](m.db_targ,times[i])
  end
  
  nT = m.mod.nvars*length(m.range)

  Jac = sparse(ro, co, vals, nT, nT)
  
  # display([ro co times vals])
  # display(full(Jac))

  return Jac

end

function residual(m::SolveModel,x)

  put_x_to_db!(m,x)
  res = eval_equations(m)
  return res

end

function residual!(m::SolveModel)

  res = eval_equations(m)
  m.db_resid[m.range,:] = reshape(res,convert(Int,m.mod.nvars),length(m.range))'
  return m

end

function jacobian(m::SolveModel,x)

  put_x_to_db!(m,x)
  Jac = eval_derivatives(m)
  return Jac

end

function jacobian(m::SolveModel)

  Jac = eval_derivatives(m)
  return Jac

end

function solve!(m::SolveModel)

  check_model(m)

  x       = get_x_from_db(m)
  res     = residual(m,x)
  resnorm = maxabs(res)
  if m.tolerance <= resnorm # Evaluate Jacobian at the starting value only if necessary
    Jac   = jacobian(m,x)
  end
  iter    = 1

  if m.display == "iter"
    @printf "Starting solver\n"
  end

  while m.tolerance <= resnorm && iter <= m.maxiter

    step = -Jac\res
    
    dump_factor = 1
    dump_iter   = 0
    
    x_is_not_ok = true
    while x_is_not_ok
    
      x_cand  = x + step
      
      try
        res = residual(m,x_cand)
        Jac = jacobian(m,x_cand)
        x   = x_cand
        x_is_not_ok = false
      catch
        dump_iter   = dump_iter + 1
        dump_factor = dump_factor/2
        step = dump_factor*step
        if maxabs(step) < eps()
          error("Can not move from the current iteration to a feasible point")
        end
      end
      
    end
    
    resnorm = maxabs(res)

    if m.display == "iter"
      if dump_iter == 0
        @printf "Iteration: %4.0f, resnorm: %12.4e\n" iter resnorm
      else
        @printf "Iteration: %4.0f, resnorm: %12.4e, dumping iterations: %2.0f\n" iter resnorm dump_iter
      end
    end

    iter = iter + 1

  end

  iter = iter - 1

  if m.display == "iter" || m.display == "final"
    @printf "Solver finished\n"
  end

  if m.display == "final"
    @printf "Iteration: %4.0f, resnorm: %12.4e\n" iter resnorm
  end

  put_x_to_db!(m,x)
  m.iter = iter
  m.resnorm = resnorm
  m.time = now()
  residual!(m)
  
end

function dac!(m::SolveModel,nstep)
  homotopy!(m::SolveModel,nstep)
end

function homotopy!(m::SolveModel,nstep)

  db_difi = m.db_targ - m.db_init
  
  @printf "Starting homotpy with %4.0f steps\n" nstep

  for i = 1:nstep

    if 1 < i
      m.db_init = m.db_targ
    end
    m.db_targ  = m.db_init + db_difi/nstep
    solve!(m)

    @printf "Homotopy step %4.0f, iteration: %4.0f, resnorm: %12.4e\n" i m.iter m.resnorm

  end
  
  @printf "Homotopy finished"

end

function put_x_to_db!(m::SolveModel,x)

  i = 0
  for t = m.range
    for j = 1:m.mod.nvars
      i = i+1
      m.db_targ[t,m.endog_ind[j]] = x[i]
    end
  end

end

function get_x_from_db(m::SolveModel)

  x = []
  
  for t = m.range
    for j = 1:m.mod.nvars
      push!(x,m.db_init[t,m.endog_ind[j]])
    end
  end

  return x

end

function endogenize!(m::SolveModel,varlist::Array{String,1})

  m.mod.endognames  = union(m.mod.endognames,   varlist)
  m.mod.exognames   = setdiff(m.mod.exognames,  varlist)
  m.mod.paramnames  = setdiff(m.mod.paramnames, varlist)
  
  m.mod.nvars   = length(m.mod.endognames)
  m.mod.nexog   = length(m.mod.exognames)
  m.mod.nparams = length(m.mod.paramnames)
  
  stacked_jac = get_stacked_jac(m.mod, m.range)
  
  m.endog_ind         = get_endog_ind(m.mod)
  m.stacked_jac_ind   = stacked_jac[1]
  m.stacked_jac_fun   = stacked_jac[2]
  m.stacked_jac_expr  = stacked_jac[3]
  
  return m

end

function endogenize!(m::SolveModel,varname::String)
  endogenize!(m, [varname])
end

function exogenize!(m::SolveModel,varlist::Array{String,1})

  m.mod.endognames  = setdiff(m.mod.endognames, varlist)
  m.mod.exognames   = union(m.mod.exognames,    varlist)
  m.mod.paramnames  = setdiff(m.mod.paramnames, varlist)
  
  m.mod.nvars   = length(m.mod.endognames)
  m.mod.nexog   = length(m.mod.exognames)
  m.mod.nparams = length(m.mod.paramnames)
  
  stacked_jac = get_stacked_jac(m.mod, m.range)
  
  m.endog_ind         = get_endog_ind(m.mod)
  m.stacked_jac_ind   = stacked_jac[1]
  m.stacked_jac_fun   = stacked_jac[2]
  m.stacked_jac_expr  = stacked_jac[3]
  
  return m

end

function exogenize!(m::SolveModel,varname::String)
  exogenize!(m, [varname])
end

# function shock!(m::SolveModel,varname::String,daterange::OrdinalRange{Date},vals::Array{Number,1},transf = "none")

function shock!(m::SolveModel,varname::String,daterange,vals,transf = "none")

  varind    = findin(m.mod.allnames,[varname])[1]
  timerange = daterange - m.firstdate + 1
  origvals  = m.db_targ[timerange,varind]
  
  if      transf == "none"
    shockvals = vals
  elseif  transf == "mult"
    shockvals = origvals .* vals
  elseif  transf == "add"
    shockvals = origvals + vals
  end
  
  m.db_targ[timerange,varind] = shockvals
  
  return m

end

# function shock!(m::SolveModel,varname::String,d::Date,vals::Number,transf = "none")
  # shock!(m,varname,d:d,[vals],transf)
# end

function get_result(m::SolveModel,ifall = false)

  db_targ   = DataBase(m.db_targ,   m.mod.allnames, m.firstdate)
  db_resid  = DataBase(m.db_resid,  m.mod.eqlabs,   m.firstdate)
  
  if !ifall
    db_targ   = db_targ[m.daterange]
    db_resid  = db_resid[m.daterange]
  end
  
  return db_targ, db_resid

end

function check_model(m::SolveModel)

  y = true
  
  y = y && m.mod.nvars == length(m.mod.eqs)
  if !y
    error("The number of endogenous variables ($(m.mod.nvars)) is different from the number of equations ($(length(m.mod.eqs)))")
  end
  
end