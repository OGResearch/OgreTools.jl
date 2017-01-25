
type Date

  freq::Int
  per::Int

  function Date(freq::Int,per::Int)

    if !(freq in [1 4 12])
      error("Frequency must be 1, 4 or 12")
    elseif per < 1
      error("Period must be positive")
    else
      new(freq,per)
    end

  end

end

function Date(freq::Int,year::Int,per::Int)

    if !(freq in [1 4 12])
      error("Frequency must be 1, 4 or 12")
    elseif year < 0
      error("Year must be non-negative")
    elseif per < 1 || freq < per
      error("Period must be between 1 and frequency")
    else
      if freq == 1
        d = Date(freq,year)
      else
        d = Date(freq,freq*year + per)
      end
    end

    return d

end

function Base.show(io::IO,d::Date)
  print(io, dat2str(d))
end

function Base.print(d::Date)
  print(dat2str(d))
end

import Base.+
function +(d::Date,n::Integer)
  d1 = Date(d.freq, d.per + n)
  return d1
end
function +(n::Integer,d::Date)
  d1 = d + n
  return d1
end

import Base.-
function -(d::Date,n::Integer)
  d1 = Date(d.freq, d.per - n)
  return d1
end
function -(d1::Date,d2::Date)
  if d1.freq != d2.freq
    error("Frequencies must be equal")
  else
    d = d1.per - d2.per
    return d
  end
end
function -(rng::OrdinalRange{Date},d::Date)
  return (rng.start-d) : rng.step : (rng.stop-d)
end

function ==(d1::Date,d2::Date)
  if d1.freq != d2.freq
    error("Frequencies must be equal")
  else
    y = d1.per == d2.per
    return y
  end
end

function Base.isless(d1::Date,d2::Date)
  if d1.freq != d2.freq
    error("Frequencies must be equal")
  else
    y = d1.per < d2.per
    return y
  end
end

# Exported functions

function yy(year::Int)
  return Date(1,year)
end

function qq(year::Int,per::Int)
  return Date(4,year,per)
end

function mm(year::Int,per::Int)
  return Date(12,year,per)
end

function ypf(d::Date)

  year  = floor(Int,d.per/d.freq)
  per   = rem(d.per,d.freq)
  if d.freq != 1 && per == 0
    per   = d.freq
    year  = year - 1
  end
  if per == 0 # Annual date
    per = 1
  end
  
  return year, per, d.freq
  
end

function dat2str(d::Date)

  year, per = ypf(d)

  if d.freq == 1
    str = @sprintf("%04.0f",year)
  elseif d.freq == 4
    str = @sprintf("%04.0fQ%1.0f",year,per)
  elseif d.freq == 12
    str = @sprintf("%04.0fM%02.0f",year,per)
  end
  
  return str

end

function dat2str(rng::StepRange{Date,Int})
  return dat2str.(collect(rng))
end