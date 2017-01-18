using BinDeps

@BinDeps.setup

deps = [
    libparser = library_dependency("libparser")
]

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

if is_linux()
  pop!(BinDeps.defaults)
end # if
