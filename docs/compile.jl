function main()

    open("docs.js", "r") do fh
        for line in eachline(fh)
            captures =  match(r"ipynb :\s+\"(\S+)\"", line)
            if captures !== nothing
                name = captures[1]

                # https://github.com/JuliaLang/IJulia.jl/issues/410
                # ipython nbconvert --to=html --template=basic --execute file.ipynb

                run(`ipython nbconvert --to=html --template=basic $name.ipynb`)

                run(`mv $name.html compiled/$name.html`)
            end
       end
    end

end

main()

