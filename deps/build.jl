using BinDeps

@BinDeps.setup

deps = [
    libparser = library_dependency("libparser")
]

# add Binaries provider to the list of available providers for Linux
if is_linux()
  push!(BinDeps.defaults, Binaries)
end # if

provides(Binaries,
         URI("https://ftp.ogresearch.com/ModelParser/libparser-win-x$(Sys.WORD_SIZE).zip"),
         libparser,
         os = :Windows)

provides(Binaries,
        URI("https://ftp.ogresearch.com/ModelParser/libparser-linux-x$(Sys.WORD_SIZE).zip"),
        libparser,
        os = :Unix)

@BinDeps.install Dict([
    (:libparser, :_jl_libparser),
])

# remove Binaries provider from the list of available providers for Linux
if is_linux()
  pop!(BinDeps.defaults)
end # if
