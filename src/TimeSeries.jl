
using Plots

type TimeSeries
  firstdate::Date
  values::Array{Number,1}
end

function TimeSeries(rng::OrdinalRange{Date},v::Number)
  return TimeSeries(rng[1],fill(v,length(rng)))
end

function TimeSeries(rng::Date,v::Number)
  return TimeSeries(rng,[v])
end

function ==(ts1::TimeSeries,ts2::TimeSeries)
  return ts1.firstdate == ts2.firstdate && ts1.values == ts2.values
end

function Base.show(io::IO,ts::TimeSeries)

  if isreal(ts)
    for i = 1:length(ts)
      print(io,ts.firstdate + i-1)
      print(io,@sprintf("%15.4e\n",ts.values[i]))
    end
  else
    for i = 1:length(ts)
      print(io,ts.firstdate + i-1)
      print(io,@sprintf("%15.4e + %11.4ei\n",real(ts.values[i]),imag(ts.values[i])))
    end
  end

end

function Base.setindex!(ts::TimeSeries,x,inds)
  ts.values[inds] = x
end
function Base.setindex!(ts::TimeSeries,x,inds::Date)
  extend!(ts,inds)
  ts.values[inds - ts.firstdate + 1] = x
end
function Base.setindex!(ts::TimeSeries,x,inds::OrdinalRange{Date})
  extend!(ts,inds)
  ts.values[inds - ts.firstdate + 1] = x
end

function Base.getindex(ts::TimeSeries,inds)
  return TimeSeries(ts.firstdate + inds - 1, ts.values[inds])
end
function Base.getindex(ts::TimeSeries,inds::UnitRange)
  return TimeSeries(ts.firstdate + inds[1] - 1, ts.values[inds])
end
function Base.getindex(ts::TimeSeries,inds::OrdinalRange)
  if inds.step == 1
    return TimeSeries(ts.firstdate + inds[1] - 1, ts.values[inds])
  else
    error("Range must be continuous")
  end
end
function Base.getindex(ts::TimeSeries,inds::Date)
  return TimeSeries(inds,ts.values[inds - ts.firstdate + 1])
end
function Base.getindex(ts::TimeSeries,inds::OrdinalRange{Date})
  if inds.step == 1
    return TimeSeries(inds[1],ts.values[inds - ts.firstdate + 1])
  else
    error("Range must be continuous")
  end
end
function Base.getindex(ts::TimeSeries,inds::BitArray{1}) # For logical indexing

  # The range will be the smallest containing every observation for TRUE-s in INDS
  fd  = findfirst(inds)
  ld  = findlast(inds)

  vals = ts[fd:ld].values
  vals[!inds[fd:ld]] = NaN # Set to values for FALSE-es to NaN

  return TimeSeries(ts.firstdate + fd - 1, vals)

end

function Base.endof(ts::TimeSeries)
  return length(ts.values)
end

# function Base.convert(Array,ts::TimeSeries)
  # return ts.values
# end

import Base.+
+(ts::TimeSeries,num::Number) = TimeSeries(ts.firstdate,ts.values + num)
+(num::Number,ts::TimeSeries) = TimeSeries(ts.firstdate,num + ts.values)
function +(ts1::TimeSeries,ts2::TimeSeries)
  ts11 = deepcopy(ts1)
  ts21 = deepcopy(ts2)
  justify!(ts11,ts21)
  ts = TimeSeries(ts11.firstdate,ts11.values + ts21.values)
  return ts
end

import Base.-
-(ts::TimeSeries,num::Number) = TimeSeries(ts.firstdate,ts.values - num)
-(num::Number,ts::TimeSeries) = TimeSeries(ts.firstdate,num - ts.values)
function -(ts1::TimeSeries,ts2::TimeSeries)
  ts11 = deepcopy(ts1)
  ts21 = deepcopy(ts2)
  justify!(ts11,ts21)
  ts = TimeSeries(ts11.firstdate,ts11.values - ts21.values)
  return ts
end

import Base.*
*(num::Number,ts::TimeSeries) = TimeSeries(ts.firstdate,ts.values * num)
*(ts::TimeSeries,num::Number) = TimeSeries(ts.firstdate,num * ts.values)
function *(ts1::TimeSeries,ts2::TimeSeries)
  ts11 = deepcopy(ts1)
  ts21 = deepcopy(ts2)
  justify!(ts11,ts21)
  ts = TimeSeries(ts11.firstdate,ts11.values .* ts21.values)
  return ts
end

import Base./
/(ts::TimeSeries,num::Number) = TimeSeries(ts.firstdate,ts.values / num)
/(num::Number,ts::TimeSeries) = TimeSeries(ts.firstdate,num ./ ts.values)
function /(ts1::TimeSeries,ts2::TimeSeries)
  ts11 = deepcopy(ts1)
  ts21 = deepcopy(ts2)
  justify!(ts11,ts21)
  ts = TimeSeries(ts11.firstdate,ts11.values ./ ts21.values)
  return ts
end

import Base.^
^(ts::TimeSeries,num::Number) = TimeSeries(ts.firstdate,ts.values .^ num)
^(num::Number,ts::TimeSeries) = TimeSeries(ts.firstdate,num .^ ts.values)
function ^(ts1::TimeSeries,ts2::TimeSeries)
  ts11 = deepcopy(ts1)
  ts21 = deepcopy(ts2)
  justify!(ts11,ts21)
  ts = TimeSeries(ts11.firstdate,ts11.values .^ ts21.values)
  return ts
end

# import Base.<
# function <(ts::TimeSeries,x::Real)
  # return TimeSeries(ts.firstdate, ts.values .< x)
# end
# import Base.<
# function <(x::Real,ts::TimeSeries)
  # return TimeSeries(ts.firstdate, x .< ts.values)
# end
# import Base.>
# function >(ts::TimeSeries,x::Real)
  # return TimeSeries(ts.firstdate, ts.values .> x)
# end
# import Base.>
# function >(x::Real,ts::TimeSeries)
  # return TimeSeries(ts.firstdate, x .> ts.values)
# end

import Base.<
function <(ts::TimeSeries,x::Real)
  return ts.values .< x
end
import Base.<
function <(x::Real,ts::TimeSeries)
  return x .< ts.values
end
import Base.>
function >(ts::TimeSeries,x::Real)
  return ts.values .> x
end
import Base.>
function >(x::Real,ts::TimeSeries)
  return x .> ts.values
end

function Base.log(ts::TimeSeries)
  return TimeSeries(ts.firstdate,log(ts.values))
end
function Base.log1p(ts::TimeSeries)
  return TimeSeries(ts.firstdate,log1p(ts.values)) # Accurate log(1+x) for small x
end
function Base.log2(ts::TimeSeries)
  return TimeSeries(ts.firstdate,log2(ts.values))
end
function Base.log10(ts::TimeSeries)
  return TimeSeries(ts.firstdate,log10(ts.values))
end

function Base.exp(ts::TimeSeries)
  return TimeSeries(ts.firstdate,exp(ts.values))
end
function Base.expm1(ts::TimeSeries)
  return TimeSeries(ts.firstdate,expm1(ts.values)) # Accurate exp(x)-1 for small x
end

function Base.abs(ts::TimeSeries)
  return TimeSeries(ts.firstdate,abs(ts.values))
end

function Base.abs2(ts::TimeSeries)
  return TimeSeries(ts.firstdate,abs2(ts.values))
end

function Base.sqrt(ts::TimeSeries)
  return TimeSeries(ts.firstdate,sqrt(ts.values))
end
function Base.cbrt(ts::TimeSeries)
  return TimeSeries(ts.firstdate,cbrt(ts.values)) # Cubic root
end

function Base.sin(ts::TimeSeries)
  return TimeSeries(ts.firstdate,sin(ts.values))
end
function Base.cos(ts::TimeSeries)
  return TimeSeries(ts.firstdate,cos(ts.values))
end
function Base.tan(ts::TimeSeries)
  return TimeSeries(ts.firstdate,tan(ts.values))
end
function Base.cot(ts::TimeSeries)
  return TimeSeries(ts.firstdate,cot(ts.values))
end
function Base.asin(ts::TimeSeries)
  return TimeSeries(ts.firstdate,asin(ts.values))
end
function Base.acos(ts::TimeSeries)
  return TimeSeries(ts.firstdate,acos(ts.values))
end
function Base.atan(ts::TimeSeries)
  return TimeSeries(ts.firstdate,atan(ts.values))
end
function Base.acot(ts::TimeSeries)
  return TimeSeries(ts.firstdate,acot(ts.values))
end

function Base.sinh(ts::TimeSeries)
  return TimeSeries(ts.firstdate,sinh(ts.values))
end
function Base.cosh(ts::TimeSeries)
  return TimeSeries(ts.firstdate,cosh(ts.values))
end
function Base.tanh(ts::TimeSeries)
  return TimeSeries(ts.firstdate,tanh(ts.values))
end
function Base.coth(ts::TimeSeries)
  return TimeSeries(ts.firstdate,coth(ts.values))
end
function Base.asinh(ts::TimeSeries)
  return TimeSeries(ts.firstdate,asinh(ts.values))
end
function Base.acosh(ts::TimeSeries)
  return TimeSeries(ts.firstdate,acosh(ts.values))
end
function Base.atanh(ts::TimeSeries)
  return TimeSeries(ts.firstdate,atanh(ts.values))
end
function Base.acoth(ts::TimeSeries)
  return TimeSeries(ts.firstdate,acoth(ts.values))
end

function Base.erf(ts::TimeSeries)
  return TimeSeries(ts.firstdate,erf(ts.values))
end
function Base.erfc(ts::TimeSeries)
  return TimeSeries(ts.firstdate,erfc(ts.values))
end
function Base.erfinv(ts::TimeSeries)
  return TimeSeries(ts.firstdate,erfinv(ts.values))
end
function Base.erfcinv(ts::TimeSeries)
  return TimeSeries(ts.firstdate,erfcinv(ts.values))
end

function Base.gamma(ts::TimeSeries)
  return TimeSeries(ts.firstdate,lgamma(ts.values))
end
function Base.lgamma(ts::TimeSeries)
  return TimeSeries(ts.firstdate,llgamma(ts.values))
end

function Base.beta(ts1::TimeSeries,ts2::TimeSeries)
  return TimeSeries(ts.firstdate,beta(ts1.values,ts2.values))
end
function Base.lbeta(ts1::TimeSeries,ts2::TimeSeries)
  return TimeSeries(ts.firstdate,lbeta(ts1.values,ts2.values))
end

function Base.diff(ts::TimeSeries)
  dts = TimeSeries(ts.firstdate,[NaN;diff(ts.values)])
end

function Base.cumsum(ts::TimeSeries,x...)
  return TimeSeries(ts.firstdate,cumsum(ts.values,x...))
end

function Base.cumprod(ts::TimeSeries,x...)
  return TimeSeries(ts.firstdate,cumprod(ts.values,x...))
end

function Base.mean(ts::TimeSeries,x...)
  return mean(ts.values,x...)
end

function Base.median(ts::TimeSeries,x...)
  return median(ts.values,x...)
end

function Base.var(ts::TimeSeries,x...)
  return var(ts.values,x...)
end

function Base.std(ts::TimeSeries,x...)
  return std(ts.values,x...)
end

function Base.cov(ts1::TimeSeries,ts2::TimeSeries)

  if ts1.firstdate.freq != ts2.firstdate.freq

    error("Frequencies must be equal")

  else

    ts11 = deepcopy(ts1)
    ts21 = deepcopy(ts2)
    justify!(ts11,ts21)

    return cov(ts11.values,ts21.values)

  end

end

function Base.cor(ts1::TimeSeries,ts2::TimeSeries)

  if ts1.firstdate.freq != ts2.firstdate.freq

    error("Frequencies must be equal")

  else

    ts11 = deepcopy(ts1)
    ts21 = deepcopy(ts2)
    justify!(ts11,ts21)

    return cor(ts11.values,ts21.values)

  end

end

function Base.real(ts::TimeSeries)
  return TimeSeries(ts.firstdate, real(ts.values))
end

function Base.imag(ts::TimeSeries)
  return TimeSeries(ts.firstdate, imag(ts.values))
end

function Base.isreal(ts::TimeSeries)
  return isreal(ts.values)
end

function Base.fft(ts::TimeSeries)
  return TimeSeries(ts.firstdate, fft(convert(Array{Real}, ts.values)))
end

function Base.ifft(ts::TimeSeries)
  return TimeSeries(ts.firstdate, ifft(convert(Array{Complex{Float64}}, ts.values)))
end

function Base.filt(a,b,ts::TimeSeries,init...)
  return TimeSeries(ts.firstdate, filt(a,b,ts.values,init...))
end

function Base.length(ts::TimeSeries)
  return length(ts.values)
end

function Base.start(ts::TimeSeries)
  return start(ts.values)
end

function Base.next(ts::TimeSeries,n::Int)
  return next(ts.values,n)
end

function Base.done(ts::TimeSeries,n::Int)
  return done(ts.values,n)
end

function Base.range(ts::TimeSeries)
  return range(ts.firstdate, length(ts.values))
end

function Plots.plot(ts::TimeSeries, step = [])

  fy, fp, freq = ypf(ts.firstdate)

  if length(ts) > 12
    # Use the first period of the first full year as the first xtick (assuming that the series is long enough)
    if fp == 1
      fxtick = 1
    else
      fxtick = Date(freq, fy+1, 1) - ts.firstdate + 1
    end
    if length(step) == 0
      step = freq
    end
    xticks = fxtick:step:length(ts)
  else
    xticks = 1:length(ts) # Use all ticks for short series
  end
  xticklabels = dat2str(range(ts)[xticks])

  plot(ts.values, xticks=(xticks, xticklabels), legend = false)

end

# Exported functions

function pct(ts::TimeSeries)
  freq      = ts.firstdate.freq;
  pchvalues = [NaN; 100*ts.values[2:end]./ts.values[1:end-1]]
  pchts     = TimeSeries(ts.firstdate,pchvalues)
  return pchts
end

function apct(ts::TimeSeries)
  freq      = ts.firstdate.freq;
  pchvalues = [repmat([NaN],freq); 100*ts.values[freq+1:end]./ts.values[1:end-freq]]
  pchts     = TimeSeries(ts.firstdate,pchvalues)
  return pchts
end

function hpf(ts::TimeSeries,lambda...)

  if length(lambda) == 0
    lambda = 100*ts.firstdate.freq^2
  end

  T = length(ts)

  A = spdiagm((ones(T-2),-2*ones(T-2),ones(T-2)),(0,1,2),T-2,T)
  hpmat = speye(T) + lambda*(A.'*A)

  return TimeSeries(ts.firstdate, hpmat\ts.values)

end

function extend!(ts::TimeSeries,d::Date)

  if d < ts.firstdate

    temp          = Array{Number,1}(ts.firstdate-d)
    temp[:]       = NaN

    ts.firstdate  = d
    ts.values     = [temp;ts.values]

  elseif ts.firstdate + length(ts) - 1 < d

    temp        = Array{Number,1}(d - (ts.firstdate + length(ts) - 1))
    temp[:]     = NaN

    ts.values   = [ts.values; temp]

  end

end

function extend!(ts::TimeSeries,rng::StepRange)

  extend!(ts,rng.start)
  extend!(ts,rng.stop)

end

function justify!(ts1::TimeSeries,ts2::TimeSeries)

  fd  = min(ts1.firstdate,ts2.firstdate)
  ld  = max(ts1.firstdate+length(ts1)-1, ts2.firstdate+length(ts2)-1)
  rng = fd:ld

  extend!(ts1,rng)
  extend!(ts2,rng)

end