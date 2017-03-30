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
         URI("https://nexus.ogresearch.com/repository/maven-releases/com/ogresearch/julia-model-parser/libparser-win-x$(Sys.WORD_SIZE)/0.0.1-2/libparser-win-x$(Sys.WORD_SIZE)-0.0.1-2.zip"),
         libparser,
         os = :Windows)

provides(Binaries,
        URI("https://nexus.ogresearch.com/repository/maven-releases/com/ogresearch/julia-model-parser/libparser-linux-x$(Sys.WORD_SIZE)/0.0.1-2/libparser-linux-x$(Sys.WORD_SIZE)-0.0.1-2.zip"),
        libparser,
        os = :Unix)

@BinDeps.install Dict([
    (:libparser, :_jl_libparser),
])

# remove Binaries provider from the list of available providers for Linux
if is_linux()
  pop!(BinDeps.defaults)
end # if
